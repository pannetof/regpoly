"""API endpoints for primitive searches (create, list, cancel, SSE progress)."""

from __future__ import annotations

import asyncio
import json

from fastapi import APIRouter, HTTPException, Query, Request
from fastapi.responses import StreamingResponse

import regpoly._regpoly_cpp as _cpp
from regpoly.core.generator import Generator, resolve_family
from regpoly.core.parametric import NotEnumerable, build_gen_enumerator
from regpoly.web.database import json_dumps, json_loads
from regpoly.web.models import PrimitiveSearchCreate
from regpoly.web.tasks.primitive import run_primitive_search


HUGE_SPACE_THRESHOLD = 10 ** 12


router = APIRouter()


def _row_field(row, name, default=None):
    """Graceful accessor for optional columns added by migration."""
    try:
        return row[name]
    except (KeyError, IndexError):
        return default


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
        "search_mode": _row_field(row, "search_mode", "random"),
        "enum_index":  _row_field(row, "enum_index", 0) or 0,
        "enum_total":  _row_field(row, "enum_total"),
        "enum_axes":   json_loads(_row_field(row, "enum_axes")),
    }


def _resolved_for_enum(body: PrimitiveSearchCreate) -> dict:
    """Merge structural + fixed (non-None) params into one dict; strip
    ``poly`` because it is the enumerated output, never a user input."""
    resolved = dict(body.structural_params)
    for key, val in body.fixed_params.items():
        if val is not None:
            resolved[key] = val
    return resolved


def _probe_k(body: PrimitiveSearchCreate, fixed: dict) -> int:
    """Instantiate the generator with a large probe L to compute k.

    Raises HTTPException(400) on invalid inputs.  The chosen probe L
    is large enough that quicktaus-style k<=L constraints always hold.
    """
    probe_L = body.L if body.L and body.L > 0 else 65536
    try:
        probe = Generator.create(body.family, probe_L, **fixed)
        return probe.k
    except Exception as exc:
        raise HTTPException(400, f"Invalid parameters: {exc}")


def _effective_L(body: PrimitiveSearchCreate, k: int) -> int:
    """Output resolution L to use for the actual search.

    Tausworthe-family generators always run at L=64 (fixed convention);
    all other families use L=k unless the caller supplied a positive
    explicit L.
    """
    if body.L and body.L > 0:
        return body.L
    family = resolve_family(body.family, body.structural_params)
    if family == "TauswortheGen":
        return 64
    return k


def _build_enumerator_or_400(body: PrimitiveSearchCreate, L: int):
    """Return an `Enumerator` for the request or raise HTTPException(400).

    Reject ``poly`` in fixed_params when in exhaustive mode.
    """
    if "poly" in body.fixed_params and body.fixed_params["poly"] is not None:
        raise HTTPException(400, detail={
            "code": "poly_in_exhaustive",
            "message": "poly is the enumerated output; it cannot be fixed "
                       "in exhaustive mode.",
        })
    resolved = _resolved_for_enum(body)
    try:
        enumerator = build_gen_enumerator(body.family, L, resolved)
    except NotEnumerable as exc:
        raise HTTPException(400, detail={
            "code": exc.reason,
            "message": f"Exhaustive mode requires additional inputs: {exc.reason}.",
        })
    if enumerator is None:
        raise HTTPException(400, detail={
            "code": "family_not_enumerable",
            "message": f"Family '{body.family}' has no exhaustive-search enumerator.",
        })
    return enumerator


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

    fixed = dict(body.structural_params)
    for key, val in body.fixed_params.items():
        if val is not None:
            fixed[key] = val

    # L is not user-facing; it is derived automatically (L=64 for
    # Tausworthe, L=k otherwise — see _effective_L).
    k = _probe_k(body, fixed)
    effective_L = _effective_L(body, k)

    # Exhaustive-mode pre-flight: build the enumerator up front to fail
    # fast, compute the space size, and store axes metadata on the row.
    enum_total_str: str | None = None
    enum_axes_json: str | None = None
    if body.search_mode == "exhaustive":
        enumerator = _build_enumerator_or_400(body, effective_L)
        total = enumerator.total
        if total > HUGE_SPACE_THRESHOLD and not body.confirm_huge:
            raise HTTPException(400, detail={
                "code": "confirm_huge_required",
                "message": (f"Exhaustive space has {total} combinations "
                            f"(> {HUGE_SPACE_THRESHOLD}); re-submit with "
                            f"confirm_huge=true to proceed."),
                "total": str(total),
            })
        enum_total_str = str(total)
        enum_axes_json = json_dumps(enumerator.axes)

    cur = await db.execute(
        """
        INSERT INTO primitive_search_run
            (family, L, k, structural_params, fixed_params,
             max_tries, max_seconds, status,
             search_mode, enum_index, enum_total, enum_axes)
        VALUES (?, ?, ?, ?, ?, ?, ?, 'pending', ?, 0, ?, ?)
        """,
        (
            family, effective_L, k,
            json_dumps(body.structural_params),
            json_dumps(body.fixed_params),
            body.max_tries, body.max_seconds,
            body.search_mode, enum_total_str, enum_axes_json,
        ),
    )
    await db.commit()
    run_id = cur.lastrowid

    pool.submit("primitive", run_id, run_primitive_search)

    return {
        "id": run_id, "status": "pending", "k": k, "L": effective_L,
        "search_mode": body.search_mode,
        "enum_total": enum_total_str,
    }


@router.post("/primitive-searches/estimate")
async def estimate_primitive_search(body: PrimitiveSearchCreate) -> dict:
    """Dry-run the exhaustive-mode configuration and report the total
    enumeration size + per-axis metadata.  Always returns a success
    response for random mode (with total=null)."""
    if body.search_mode != "exhaustive":
        return {"search_mode": body.search_mode, "total": None, "axes": []}
    # Probe for k so the enumerator sees a self-consistent L.
    fixed = dict(body.structural_params)
    for key, val in body.fixed_params.items():
        if val is not None:
            fixed[key] = val
    k = _probe_k(body, fixed)
    enumerator = _build_enumerator_or_400(body, _effective_L(body, k))
    return {
        "search_mode": "exhaustive",
        "total": str(enumerator.total),
        "axes": enumerator.axes,
        "huge": enumerator.total > HUGE_SPACE_THRESHOLD,
    }


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
        "max_tries, max_seconds, search_mode, enum_total, enum_axes "
        "FROM primitive_search_run WHERE id = ?",
        (run_id,),
    ) as cur:
        row = await cur.fetchone()
    if row is None:
        raise HTTPException(404, f"Search {run_id} not found")
    cur = await db.execute(
        """
        INSERT INTO primitive_search_run
            (family, L, k, structural_params, fixed_params,
             max_tries, max_seconds, status,
             search_mode, enum_index, enum_total, enum_axes)
        VALUES (?, ?, ?, ?, ?, ?, ?, 'pending', ?, 0, ?, ?)
        """,
        (row["family"], row["L"], row["k"],
         row["structural_params"], row["fixed_params"],
         row["max_tries"], row["max_seconds"],
         _row_field(row, "search_mode", "random"),
         _row_field(row, "enum_total"),
         _row_field(row, "enum_axes")),
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


@router.delete("/primitive-searches/{run_id}/generators")
async def delete_primitive_search_generators(
    request: Request, run_id: int
) -> dict:
    """Delete every primitive generator attributed to this search run,
    then delete the run row itself so the search disappears from
    ``/searches`` entirely."""
    db = request.app.state.db
    cur_gen = await db.execute(
        "DELETE FROM primitive_generator WHERE search_run_id = ?",
        (run_id,),
    )
    generators_deleted = cur_gen.rowcount
    cur_run = await db.execute(
        "DELETE FROM primitive_search_run WHERE id = ?",
        (run_id,),
    )
    run_deleted = cur_run.rowcount
    await db.commit()
    return {
        "run_id": run_id,
        "deleted": generators_deleted,
        "run_deleted": bool(run_deleted),
    }


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


