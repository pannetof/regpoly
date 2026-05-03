"""API endpoints for tempering searches (create, list, cancel, SSE progress)."""

from __future__ import annotations

import asyncio
import json

from fastapi import APIRouter, HTTPException, Query, Request
from fastapi.responses import StreamingResponse

from regpoly.web.database import json_dumps, json_loads
from regpoly.web.models import TemperingSearchCreate
from regpoly.web.param_format import (
    format_gen_params, format_tempering_list,
)
from regpoly.web.tasks.tempering import run_tempering_search


router = APIRouter()


def _row_to_run(row) -> dict:
    return {
        "id": row["id"],
        "test_type": row["test_type"],
        "test_config": json_loads(row["test_config"]),
        "Lmax": row["Lmax"],
        "nb_tries": row["nb_tries"],
        "optimizer_config": json_loads(row["optimizer_config"]),
        "status": row["status"],
        "combos_total": row["combos_total"],
        "combos_done": row["combos_done"],
        "best_se": row["best_se"],
        "elapsed_seconds": row["elapsed_seconds"],
        "error_message": row["error_message"],
        "created_at": row["created_at"],
        "started_at": row["started_at"],
        "finished_at": row["finished_at"],
    }


@router.post("/tempering-searches")
async def create_tempering_search(
    request: Request, body: TemperingSearchCreate
) -> dict:
    db = request.app.state.db
    pool = request.app.state.pool

    # Validate the user-supplied generator IDs are all present.
    all_ids: set[int] = set()
    for c in body.components:
        all_ids.update(c.generator_ids)
    if not all_ids:
        raise HTTPException(400, "No generators selected")

    async with db.execute(
        f"SELECT id FROM primitive_generator WHERE id IN "
        f"({','.join('?' * len(all_ids))})", list(all_ids)
    ) as cur:
        rows = await cur.fetchall()
    found = {r[0] for r in rows}
    missing = all_ids - found
    if missing:
        raise HTTPException(400, f"Unknown generator IDs: {sorted(missing)}")

    test_type = body.test.get("type", "equidistribution")

    cur = await db.execute(
        """
        INSERT INTO tempering_search_run
            (test_type, test_config, Lmax, nb_tries, optimizer_config, status)
        VALUES (?, ?, ?, ?, ?, 'pending')
        """,
        (
            test_type,
            json_dumps(body.test),
            body.Lmax,
            body.nb_tries,
            json_dumps(body.optimizer) if body.optimizer else None,
        ),
    )
    run_id = cur.lastrowid

    for idx, c in enumerate(body.components):
        comp_cur = await db.execute(
            """
            INSERT INTO tempering_search_component
                (search_run_id, component_index, shared_with_component,
                 tempering_config)
            VALUES (?, ?, ?, ?)
            """,
            (
                run_id, idx, c.shared_with_component,
                json_dumps(c.tempering),
            ),
        )
        comp_id = comp_cur.lastrowid
        for gen_id in c.generator_ids:
            await db.execute(
                "INSERT INTO tempering_search_generator"
                "(component_id, generator_id) VALUES (?, ?)",
                (comp_id, gen_id),
            )

    await db.commit()

    pool.submit("tempering", run_id, run_tempering_search)

    return {"id": run_id, "status": "pending"}


@router.get("/tempering-searches")
async def list_tempering_searches(
    request: Request,
    status: str | None = None,
    family: str | None = None,
    k: int | None = None,
    page: int = Query(1, ge=1),
    per_page: int = Query(50, ge=1, le=500),
) -> dict:
    """List tempering-search runs.

    `family` / `k` filter by "any component in this search uses a
    primitive generator of that family / degree".
    """
    db = request.app.state.db

    where = []
    params: list = []
    joins = ""
    if status:
        where.append("tsr.status = ?")
        params.append(status)
    if family or k is not None:
        joins = (
            " JOIN tempering_search_component tsc "
            "    ON tsc.search_run_id = tsr.id "
            " JOIN tempering_search_generator tsg "
            "    ON tsg.component_id = tsc.id "
            " JOIN primitive_generator pg "
            "    ON pg.id = tsg.generator_id "
        )
        if family:
            where.append("pg.family = ?")
            params.append(family)
        if k is not None:
            where.append("pg.k = ?")
            params.append(k)

    where_clause = ("WHERE " + " AND ".join(where)) if where else ""

    count_sql = (
        f"SELECT COUNT(DISTINCT tsr.id) FROM tempering_search_run tsr"
        f"{joins}{(' ' + where_clause) if where_clause else ''}"
    )
    async with db.execute(count_sql, params) as cur:
        total = (await cur.fetchone())[0]

    offset = (page - 1) * per_page
    list_sql = (
        f"SELECT DISTINCT tsr.* FROM tempering_search_run tsr"
        f"{joins}{(' ' + where_clause) if where_clause else ''}"
        f" ORDER BY tsr.id DESC LIMIT ? OFFSET ?"
    )
    async with db.execute(list_sql, [*params, per_page, offset]) as cur:
        rows = await cur.fetchall()

    return {
        "items": [_row_to_run(r) for r in rows],
        "total": total,
        "page": page,
        "per_page": per_page,
    }


@router.get("/tempering-searches/{run_id}")
async def get_tempering_search(request: Request, run_id: int) -> dict:
    db = request.app.state.db
    async with db.execute(
        "SELECT * FROM tempering_search_run WHERE id = ?", (run_id,)
    ) as cur:
        row = await cur.fetchone()
    if row is None:
        raise HTTPException(404, f"Search {run_id} not found")

    # Attach components
    async with db.execute(
        "SELECT * FROM tempering_search_component WHERE search_run_id = ? "
        "ORDER BY component_index",
        (run_id,),
    ) as cur:
        comp_rows = await cur.fetchall()

    components = []
    for c in comp_rows:
        async with db.execute(
            "SELECT g.id, g.family, g.L, g.k, g.all_params, "
            "       g.search_run_id "
            "FROM tempering_search_generator tsg "
            "JOIN primitive_generator g ON g.id = tsg.generator_id "
            "WHERE tsg.component_id = ? "
            "ORDER BY g.id",
            (c["id"],),
        ) as cur:
            rows = await cur.fetchall()
        gen_ids = [r["id"] for r in rows]
        gen_records = [
            {
                "id": r["id"], "family": r["family"], "L": r["L"],
                "k": r["k"],
                "all_params": format_gen_params(
                    r["family"], json_loads(r["all_params"])),
                "search_run_id": r["search_run_id"],
            }
            for r in rows
        ]
        components.append({
            "component_index": c["component_index"],
            "shared_with_component": c["shared_with_component"],
            "tempering_config": format_tempering_list(
                json_loads(c["tempering_config"])),
            "generator_ids": gen_ids,
            "generators": gen_records,
        })

    result = _row_to_run(row)
    result["components"] = components
    return result


@router.post("/tempering-searches/{run_id}/cancel")
async def cancel_tempering_search(request: Request, run_id: int) -> dict:
    db = request.app.state.db
    status = await _current_status(db, run_id)
    if status is None:
        raise HTTPException(404, f"Search {run_id} not found")
    if status not in ("pending", "running", "paused"):
        return {"id": run_id, "status": status}
    await db.execute(
        "UPDATE tempering_search_run SET status='cancelled' WHERE id = ?",
        (run_id,),
    )
    await db.commit()
    return {"id": run_id, "status": "cancelled"}


@router.post("/tempering-searches/{run_id}/pause")
async def pause_tempering_search(request: Request, run_id: int) -> dict:
    db = request.app.state.db
    status = await _current_status(db, run_id)
    if status is None:
        raise HTTPException(404, f"Search {run_id} not found")
    if status not in ("pending", "running"):
        return {"id": run_id, "status": status}
    await db.execute(
        "UPDATE tempering_search_run SET status='paused' WHERE id = ?",
        (run_id,),
    )
    await db.commit()
    return {"id": run_id, "status": "paused"}


@router.post("/tempering-searches/{run_id}/resume")
async def resume_tempering_search(request: Request, run_id: int) -> dict:
    db = request.app.state.db
    pool = request.app.state.pool
    status = await _current_status(db, run_id)
    if status is None:
        raise HTTPException(404, f"Search {run_id} not found")
    if status not in ("paused", "cancelled", "failed"):
        return {"id": run_id, "status": status}
    await db.execute(
        "UPDATE tempering_search_run SET status='pending', "
        "error_message=NULL, finished_at=NULL WHERE id = ?",
        (run_id,),
    )
    await db.commit()
    pool.submit("tempering", run_id, run_tempering_search)
    return {"id": run_id, "status": "pending"}


@router.post("/tempering-searches/{run_id}/restart")
async def restart_tempering_search(request: Request, run_id: int) -> dict:
    """Clone a tempering search into a fresh run with the same config."""
    db = request.app.state.db
    pool = request.app.state.pool
    async with db.execute(
        "SELECT * FROM tempering_search_run WHERE id = ?", (run_id,)
    ) as cur:
        row = await cur.fetchone()
    if row is None:
        raise HTTPException(404, f"Search {run_id} not found")

    cur = await db.execute(
        """
        INSERT INTO tempering_search_run
            (test_type, test_config, Lmax, nb_tries, optimizer_config, status)
        VALUES (?, ?, ?, ?, ?, 'pending')
        """,
        (row["test_type"], row["test_config"], row["Lmax"],
         row["nb_tries"], row["optimizer_config"]),
    )
    new_id = cur.lastrowid

    # Clone components + generator pools
    async with db.execute(
        "SELECT * FROM tempering_search_component WHERE search_run_id = ? "
        "ORDER BY component_index",
        (run_id,),
    ) as c1:
        comp_rows = await c1.fetchall()
    for c in comp_rows:
        cur2 = await db.execute(
            """
            INSERT INTO tempering_search_component
                (search_run_id, component_index, shared_with_component,
                 tempering_config)
            VALUES (?, ?, ?, ?)
            """,
            (new_id, c["component_index"], c["shared_with_component"],
             c["tempering_config"]),
        )
        new_comp_id = cur2.lastrowid
        async with db.execute(
            "SELECT generator_id FROM tempering_search_generator "
            "WHERE component_id = ?", (c["id"],)
        ) as c3:
            for g in await c3.fetchall():
                await db.execute(
                    "INSERT INTO tempering_search_generator"
                    "(component_id, generator_id) VALUES (?, ?)",
                    (new_comp_id, g[0]),
                )
    await db.commit()

    pool.submit("tempering", new_id, run_tempering_search)
    return {"id": new_id, "status": "pending", "cloned_from": run_id}


async def _current_status(db, run_id: int) -> str | None:
    async with db.execute(
        "SELECT status FROM tempering_search_run WHERE id = ?", (run_id,)
    ) as cur:
        row = await cur.fetchone()
    return row["status"] if row else None


@router.delete("/tempering-searches/{run_id}/results")
async def delete_tempering_search_results(
    request: Request, run_id: int
) -> dict:
    """Delete every combined generator attributed to this search run
    (cascades to components and test results), then delete the search
    run row itself (which cascades to its component pool rows)."""
    db = request.app.state.db
    cur_tg = await db.execute(
        "DELETE FROM tested_generator WHERE search_run_id = ?",
        (run_id,),
    )
    tested_deleted = cur_tg.rowcount
    cur_run = await db.execute(
        "DELETE FROM tempering_search_run WHERE id = ?",
        (run_id,),
    )
    run_deleted = cur_run.rowcount
    await db.commit()
    return {
        "run_id": run_id,
        "deleted": tested_deleted,
        "run_deleted": bool(run_deleted),
    }


@router.get("/tempering-searches/{run_id}/results")
async def tempering_search_results(request: Request, run_id: int) -> dict:
    from regpoly.web.routes.tested_generators import _fetch_tested
    db = request.app.state.db

    async with db.execute(
        "SELECT id FROM tested_generator WHERE search_run_id = ? "
        "ORDER BY id DESC",
        (run_id,),
    ) as cur:
        rows = await cur.fetchall()

    items = []
    for r in rows:
        items.append(await _fetch_tested(db, r["id"]))
    return {"items": items, "total": len(items)}


@router.get("/tempering-searches/{run_id}/progress")
async def tempering_search_progress_sse(
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
                    "WHERE search_type='tempering' AND search_run_id=? "
                    "AND id > ? ORDER BY id ASC",
                    (run_id, last_id),
                ) as cur:
                    rows = await cur.fetchall()

                for row in rows:
                    payload = {
                        "current_info": json_loads(row["current_info"]),
                        "message": row["message"],
                        "updated_at": row["updated_at"],
                    }
                    yield f"data: {json.dumps(payload)}\n\n"
                    last_id = row["id"]

                async with conn.execute(
                    "SELECT status, combos_done, combos_total, best_se, "
                    "elapsed_seconds "
                    "FROM tempering_search_run WHERE id = ?", (run_id,)
                ) as cur:
                    r = await cur.fetchone()
                if r:
                    yield ("data: " + json.dumps({
                        "status": r["status"],
                        "combos_done": r["combos_done"],
                        "combos_total": r["combos_total"],
                        "best_se": r["best_se"],
                        "elapsed_seconds": r["elapsed_seconds"],
                    }) + "\n\n")
                    if r["status"] in ("completed", "cancelled", "failed"):
                        yield f"event: end\ndata: {json.dumps({'status': r['status']})}\n\n"
                        break

                await asyncio.sleep(poll)

    return StreamingResponse(
        event_stream(),
        media_type="text/event-stream",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
    )
