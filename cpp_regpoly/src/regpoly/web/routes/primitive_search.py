"""API endpoints for primitive searches (create, list, cancel, SSE progress)."""

from __future__ import annotations

import asyncio
import json

from fastapi import APIRouter, HTTPException, Query, Request
from fastapi.responses import StreamingResponse

import regpoly._regpoly_cpp as _cpp
from regpoly.generateur import Generateur, resolve_family
from regpoly.web.database import json_dumps, json_loads
from regpoly.web.models import PrimitiveSearchCreate
from regpoly.web.tasks.primitive import run_primitive_search


router = APIRouter()


def _row_to_run(row) -> dict:
    return {
        "id": row["id"],
        "family": row["family"],
        "L": row["L"],
        "k": row["k"],
        "structural_params": json_loads(row["structural_params"]),
        "fixed_params": json_loads(row["fixed_params"]),
        "max_tries": row["max_tries"],
        "max_seconds": row["max_seconds"],
        "status": row["status"],
        "tries_done": row["tries_done"],
        "found_count": row["found_count"],
        "elapsed_seconds": row["elapsed_seconds"],
        "error_message": row["error_message"],
        "created_at": row["created_at"],
        "started_at": row["started_at"],
        "finished_at": row["finished_at"],
    }


@router.post("/primitive-searches")
async def create_primitive_search(
    request: Request, body: PrimitiveSearchCreate
) -> dict:
    db = request.app.state.db
    pool = request.app.state.pool

    family_raw = body.family
    family = resolve_family(family_raw, body.structural_params)

    # Probe once to compute k; also validates the user input.
    try:
        specs = _cpp.get_gen_param_specs(family)
    except Exception as exc:
        raise HTTPException(400, f"Unknown family '{family}': {exc}")

    try:
        fixed = dict(body.structural_params)
        for key, val in body.fixed_params.items():
            if val is not None:
                fixed[key] = val
        # Generateur.create does its own randomization and injects L so
        # samplers like tausworthe_poly see the correct output width.
        probe = Generateur.create(family_raw, body.L, **fixed)
        k = probe.k
    except Exception as exc:
        raise HTTPException(400, f"Invalid parameters: {exc}")

    cur = await db.execute(
        """
        INSERT INTO primitive_search_run
            (family, L, k, structural_params, fixed_params,
             max_tries, max_seconds, status)
        VALUES (?, ?, ?, ?, ?, ?, ?, 'pending')
        """,
        (
            family, body.L, k,
            json_dumps(body.structural_params),
            json_dumps(body.fixed_params),
            body.max_tries, body.max_seconds,
        ),
    )
    await db.commit()
    run_id = cur.lastrowid

    pool.submit("primitive", run_id, run_primitive_search)

    return {"id": run_id, "status": "pending", "k": k}


@router.get("/primitive-searches")
async def list_primitive_searches(
    request: Request,
    status: str | None = None,
    family: str | None = None,
    k: int | None = None,
    page: int = Query(1, ge=1),
    per_page: int = Query(50, ge=1, le=500),
) -> dict:
    db = request.app.state.db

    where = []
    params: list = []
    if status:
        where.append("status = ?")
        params.append(status)
    if family:
        where.append("family = ?")
        params.append(family)
    if k is not None:
        where.append("k = ?")
        params.append(k)

    where_clause = ("WHERE " + " AND ".join(where)) if where else ""

    async with db.execute(
        f"SELECT COUNT(*) FROM primitive_search_run {where_clause}", params
    ) as cur:
        row = await cur.fetchone()
        total = row[0]

    offset = (page - 1) * per_page
    async with db.execute(
        f"SELECT * FROM primitive_search_run {where_clause} "
        f"ORDER BY id DESC LIMIT ? OFFSET ?",
        [*params, per_page, offset],
    ) as cur:
        rows = await cur.fetchall()

    return {
        "items": [_row_to_run(r) for r in rows],
        "total": total,
        "page": page,
        "per_page": per_page,
    }


@router.get("/primitive-searches/{run_id}")
async def get_primitive_search(request: Request, run_id: int) -> dict:
    db = request.app.state.db
    async with db.execute(
        "SELECT * FROM primitive_search_run WHERE id = ?", (run_id,)
    ) as cur:
        row = await cur.fetchone()
    if row is None:
        raise HTTPException(404, f"Search {run_id} not found")
    return _row_to_run(row)


@router.post("/primitive-searches/{run_id}/cancel")
async def cancel_primitive_search(request: Request, run_id: int) -> dict:
    db = request.app.state.db
    status = await _current_status(db, run_id)
    if status is None:
        raise HTTPException(404, f"Search {run_id} not found")
    if status not in ("pending", "running", "paused"):
        return {"id": run_id, "status": status}
    await db.execute(
        "UPDATE primitive_search_run SET status='cancelled' WHERE id = ?",
        (run_id,),
    )
    await db.commit()
    return {"id": run_id, "status": "cancelled"}


@router.post("/primitive-searches/{run_id}/pause")
async def pause_primitive_search(request: Request, run_id: int) -> dict:
    db = request.app.state.db
    status = await _current_status(db, run_id)
    if status is None:
        raise HTTPException(404, f"Search {run_id} not found")
    if status not in ("pending", "running"):
        return {"id": run_id, "status": status}
    await db.execute(
        "UPDATE primitive_search_run SET status='paused' WHERE id = ?",
        (run_id,),
    )
    await db.commit()
    return {"id": run_id, "status": "paused"}


@router.post("/primitive-searches/{run_id}/resume")
async def resume_primitive_search(request: Request, run_id: int) -> dict:
    db = request.app.state.db
    pool = request.app.state.pool
    status = await _current_status(db, run_id)
    if status is None:
        raise HTTPException(404, f"Search {run_id} not found")
    if status not in ("paused", "cancelled", "failed"):
        return {"id": run_id, "status": status}
    await db.execute(
        "UPDATE primitive_search_run SET status='pending', "
        "error_message=NULL, finished_at=NULL WHERE id = ?",
        (run_id,),
    )
    await db.commit()
    pool.submit("primitive", run_id, run_primitive_search)
    return {"id": run_id, "status": "pending"}


@router.post("/primitive-searches/{run_id}/restart")
async def restart_primitive_search(request: Request, run_id: int) -> dict:
    """Create a NEW primitive search run with the same config as {run_id}."""
    db = request.app.state.db
    pool = request.app.state.pool
    async with db.execute(
        "SELECT family, L, k, structural_params, fixed_params, "
        "max_tries, max_seconds FROM primitive_search_run WHERE id = ?",
        (run_id,),
    ) as cur:
        row = await cur.fetchone()
    if row is None:
        raise HTTPException(404, f"Search {run_id} not found")
    cur = await db.execute(
        """
        INSERT INTO primitive_search_run
            (family, L, k, structural_params, fixed_params,
             max_tries, max_seconds, status)
        VALUES (?, ?, ?, ?, ?, ?, ?, 'pending')
        """,
        (row["family"], row["L"], row["k"],
         row["structural_params"], row["fixed_params"],
         row["max_tries"], row["max_seconds"]),
    )
    await db.commit()
    new_id = cur.lastrowid
    pool.submit("primitive", new_id, run_primitive_search)
    return {"id": new_id, "status": "pending", "cloned_from": run_id}


async def _current_status(db, run_id: int) -> str | None:
    async with db.execute(
        "SELECT status FROM primitive_search_run WHERE id = ?", (run_id,)
    ) as cur:
        row = await cur.fetchone()
    return row["status"] if row else None


@router.get("/primitive-searches/{run_id}/generators")
async def primitive_search_generators(
    request: Request, run_id: int,
    page: int = Query(1, ge=1),
    per_page: int = Query(50, ge=1, le=500),
) -> dict:
    from regpoly.web.routes.generators import _row_to_generator
    db = request.app.state.db

    async with db.execute(
        "SELECT COUNT(*) FROM primitive_generator WHERE search_run_id = ?",
        (run_id,),
    ) as cur:
        total = (await cur.fetchone())[0]

    offset = (page - 1) * per_page
    async with db.execute(
        "SELECT * FROM primitive_generator WHERE search_run_id = ? "
        "ORDER BY id DESC LIMIT ? OFFSET ?",
        (run_id, per_page, offset),
    ) as cur:
        rows = await cur.fetchall()

    return {
        "items": [_row_to_generator(r) for r in rows],
        "total": total,
        "page": page,
        "per_page": per_page,
    }


@router.get("/primitive-searches/{run_id}/progress")
async def primitive_search_progress_sse(
    request: Request, run_id: int
) -> StreamingResponse:
    settings = request.app.state.settings
    db_path = settings.db_path
    poll = settings.progress_poll_seconds

    async def event_stream():
        import aiosqlite
        async with aiosqlite.connect(db_path) as conn:
            conn.row_factory = aiosqlite.Row
            last_id = 0
            while True:
                if await request.is_disconnected():
                    break

                async with conn.execute(
                    "SELECT * FROM search_progress "
                    "WHERE search_type='primitive' AND search_run_id=? "
                    "AND id > ? ORDER BY id ASC",
                    (run_id, last_id),
                ) as cur:
                    rows = await cur.fetchall()

                for row in rows:
                    payload = {
                        "tries_done": row["tries_done"],
                        "found_count": row["found_count"],
                        "current_info": json_loads(row["current_info"]),
                        "message": row["message"],
                        "updated_at": row["updated_at"],
                    }
                    yield f"data: {json.dumps(payload)}\n\n"
                    last_id = row["id"]

                # Emit terminal event if the run has finished
                async with conn.execute(
                    "SELECT status FROM primitive_search_run WHERE id = ?",
                    (run_id,),
                ) as cur:
                    r = await cur.fetchone()
                if r and r["status"] in ("completed", "cancelled", "failed"):
                    yield f"event: end\ndata: {json.dumps({'status': r['status']})}\n\n"
                    break

                await asyncio.sleep(poll)

    return StreamingResponse(
        event_stream(),
        media_type="text/event-stream",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
    )


