"""API endpoints for tested combined generators."""

from __future__ import annotations

from fastapi import APIRouter, HTTPException, Query, Request

from regpoly_web.database import json_loads
from regpoly_web.param_format import (
    format_gen_params,
    format_tempering_list,
)
from regpoly_web.results import read_typed_results_async

router = APIRouter()


async def _fetch_tested(db, tg_id: int) -> dict | None:
    async with db.execute(
        "SELECT * FROM tested_generator WHERE id = ?", (tg_id,)
    ) as cur:
        row = await cur.fetchone()
    if row is None:
        return None

    async with db.execute(
        "SELECT * FROM tested_generator_component "
        "WHERE tested_gen_id = ? ORDER BY component_index",
        (tg_id,),
    ) as cur:
        comp_rows = await cur.fetchall()

    results = await read_typed_results_async(db, tg_id)

    keys = row.keys() if hasattr(row, "keys") else []
    return {
        "id": row["id"],
        "search_run_id": row["search_run_id"],
        "Lmax": row["Lmax"],
        "k_g": row["k_g"],
        "J": row["J"],
        "created_at": row["created_at"],
        "library_id": row["library_id"] if "library_id" in keys else None,
        "components": [
            {
                "component_index": c["component_index"],
                "generator_id": c["generator_id"],
                "family": c["family"],
                "L": c["L"],
                "k": c["k"],
                "all_params": format_gen_params(
                    c["family"], json_loads(c["all_params"])),
                "tempering_params": format_tempering_list(
                    json_loads(c["tempering_params"])),
            }
            for c in comp_rows
        ],
        "results": results,
    }


@router.get("/tested-generators")
async def list_tested_generators(
    request: Request,
    family: str | None = None,
    min_k_g: int | None = None,
    max_k_g: int | None = None,
    test_type: str | None = None,
    max_se: int | None = None,
    is_me: bool | None = None,
    search_run_id: int | None = None,
    page: int = Query(1, ge=1),
    per_page: int = Query(50, ge=1, le=500),
) -> dict:
    db = request.app.state.db

    joins = []
    where = []
    params: list = []

    if family:
        joins.append("JOIN tested_generator_component tgc "
                     "ON tgc.tested_gen_id = tg.id")
        where.append("tgc.family = ?")
        params.append(family)
    # Typed-table joins (Phase 5.4c). The legacy join on test_result was
    # generic; the typed schema requires picking the right table per
    # test_type. test_type is dropped from the legacy "any of three" mode
    # because the UI always passes a concrete value here in practice.
    join_table = None
    if test_type:
        if test_type == "equidistribution":
            join_table = "equidistribution_result"
        elif test_type == "collision_free":
            join_table = "collision_free_result"
        elif test_type == "tuplets":
            join_table = "tuplets_result"
    if join_table is None and (max_se is not None or is_me is not None):
        # No explicit test_type but user is filtering on equid-shaped
        # fields; default to equidistribution_result (most common case).
        join_table = "equidistribution_result"
    if join_table is not None:
        joins.append(f"JOIN {join_table} tr ON tr.tested_gen_id = tg.id")
        if max_se is not None:
            # Equid: se. CF: secf. Tuplets: no se filter (skip).
            if join_table == "equidistribution_result":
                where.append("tr.se <= ?")
                params.append(max_se)
            elif join_table == "collision_free_result":
                where.append("tr.secf <= ?")
                params.append(max_se)
        if is_me is not None and join_table == "equidistribution_result":
            # is_me ⇔ verified=1 AND se=0 (mirrors EquidistributionResults.is_me()).
            if is_me:
                where.append("tr.verified = 1 AND tr.se = 0")
            else:
                where.append("(tr.verified = 0 OR tr.se != 0)")
    if min_k_g is not None:
        where.append("tg.k_g >= ?")
        params.append(min_k_g)
    if max_k_g is not None:
        where.append("tg.k_g <= ?")
        params.append(max_k_g)
    if search_run_id is not None:
        where.append("tg.search_run_id = ?")
        params.append(search_run_id)

    join_clause = " ".join(joins)
    where_clause = ("WHERE " + " AND ".join(where)) if where else ""

    count_sql = (
        f"SELECT COUNT(DISTINCT tg.id) FROM tested_generator tg "
        f"{join_clause} {where_clause}"
    )
    async with db.execute(count_sql, params) as cur:
        total = (await cur.fetchone())[0]

    offset = (page - 1) * per_page
    list_sql = (
        f"SELECT DISTINCT tg.id FROM tested_generator tg "
        f"{join_clause} {where_clause} "
        f"ORDER BY tg.id DESC LIMIT ? OFFSET ?"
    )
    async with db.execute(list_sql, [*params, per_page, offset]) as cur:
        rows = await cur.fetchall()

    items = []
    for r in rows:
        item = await _fetch_tested(db, r[0])
        if item:
            items.append(item)

    return {
        "items": items,
        "total": total,
        "page": page,
        "per_page": per_page,
    }


@router.get("/tested-generators/{tg_id}")
async def get_tested_generator(request: Request, tg_id: int) -> dict:
    db = request.app.state.db
    item = await _fetch_tested(db, tg_id)
    if item is None:
        raise HTTPException(404, f"Tested generator {tg_id} not found")
    return item


@router.delete("/tested-generators/{tg_id}")
async def delete_tested_generator(request: Request, tg_id: int) -> dict:
    db = request.app.state.db
    cur = await db.execute(
        "DELETE FROM tested_generator WHERE id = ?", (tg_id,)
    )
    await db.commit()
    if cur.rowcount == 0:
        raise HTTPException(404, f"Tested generator {tg_id} not found")
    return {"deleted": tg_id}
