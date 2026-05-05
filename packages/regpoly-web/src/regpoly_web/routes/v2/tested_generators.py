"""GET /api/v2/tested-generators — chip facets + has_results filter.

`?has_results=true` joins the v2 typed-result tables
(equidistribution_result, collision_free_result, tuplets_result) so
only tested generators with at least one typed analysis row appear.
"""

from __future__ import annotations

from fastapi import APIRouter, Query, Request
from pydantic import BaseModel

router = APIRouter()


class _TestedRow(BaseModel):
    id: int
    search_run_id: int | None
    Lmax: int | None
    k_g: int | None
    J: int | None
    library_id: str | None


class TestedGeneratorList(BaseModel):
    rows: list[_TestedRow]
    total: int


@router.get("/tested-generators", response_model=TestedGeneratorList)
async def v2_list_tested_generators(
    request: Request,
    family: str | None = None,
    k_g_min: int | None = None,
    k_g_max: int | None = None,
    test_type: str | None = None,
    max_sum_delta: int | None = None,
    is_me: bool | None = None,
    has_results: bool | None = None,
    limit: int = Query(50, ge=1, le=500),
    offset: int = Query(0, ge=0),
) -> TestedGeneratorList:
    db = request.app.state.db

    where: list[str] = []
    params: list = []
    joins: list[str] = []

    if family:
        joins.append(
            "INNER JOIN tested_generator_component AS tgc "
            "ON tgc.tested_gen_id = tg.id AND tgc.component_index = 0"
        )
        where.append("tgc.family = ?")
        params.append(family)
    if k_g_min is not None:
        where.append("tg.k_g >= ?")
        params.append(k_g_min)
    if k_g_max is not None:
        where.append("tg.k_g <= ?")
        params.append(k_g_max)
    if test_type or has_results is True or max_sum_delta is not None:
        # Join the equidistribution_result table for any of these
        # filters. We accept that test_type currently only covers the
        # equidistribution result; collision_free / tuplets get added
        # in P3+.
        joins.append(
            "LEFT JOIN equidistribution_result AS er "
            "ON er.tested_gen_id = tg.id"
        )
    if test_type == "equidistribution":
        where.append("er.id IS NOT NULL")
    if max_sum_delta is not None:
        where.append("er.se <= ?")
        params.append(int(max_sum_delta))
    if has_results is True:
        # tested_generator has at least one row in any v2 typed table.
        where.append(
            "(EXISTS (SELECT 1 FROM equidistribution_result e "
            "         WHERE e.tested_gen_id = tg.id) "
            " OR EXISTS (SELECT 1 FROM collision_free_result c "
            "            WHERE c.tested_gen_id = tg.id) "
            " OR EXISTS (SELECT 1 FROM tuplets_result t "
            "            WHERE t.tested_gen_id = tg.id))"
        )
    elif has_results is False:
        where.append(
            "NOT EXISTS (SELECT 1 FROM equidistribution_result e "
            "            WHERE e.tested_gen_id = tg.id) "
            "AND NOT EXISTS (SELECT 1 FROM collision_free_result c "
            "                WHERE c.tested_gen_id = tg.id) "
            "AND NOT EXISTS (SELECT 1 FROM tuplets_result t "
            "                WHERE t.tested_gen_id = tg.id)"
        )
    if is_me is not None:
        # ME = "maximally equidistributed": Σδ == 0 across the whole
        # equidistribution profile.
        joins.append(
            "LEFT JOIN equidistribution_result AS er_me "
            "ON er_me.tested_gen_id = tg.id"
        ) if "er_me" not in " ".join(joins) else None
        where.append("er_me.se " + ("= 0" if is_me else "> 0"))

    join_clause = " ".join(j for j in joins if j)
    where_clause = ("WHERE " + " AND ".join(where)) if where else ""

    async with db.execute(
        f"SELECT COUNT(DISTINCT tg.id) FROM tested_generator AS tg "
        f"{join_clause} {where_clause}",
        params,
    ) as cur:
        total = (await cur.fetchone())[0]

    list_sql = (
        f"SELECT DISTINCT tg.id, tg.search_run_id, tg.Lmax, tg.k_g, tg.J, "
        f"       tg.library_id "
        f"FROM tested_generator AS tg {join_clause} {where_clause} "
        f"ORDER BY tg.id DESC LIMIT ? OFFSET ?"
    )
    async with db.execute(list_sql, [*params, limit, offset]) as cur:
        rows = await cur.fetchall()

    return TestedGeneratorList(
        rows=[
            _TestedRow(
                id=r["id"], search_run_id=r["search_run_id"],
                Lmax=r["Lmax"], k_g=r["k_g"], J=r["J"],
                library_id=r["library_id"] if "library_id" in r.keys() else None,
            )
            for r in rows
        ],
        total=total,
    )
