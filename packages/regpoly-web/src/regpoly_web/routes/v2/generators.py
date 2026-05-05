"""GET /api/v2/generators — server-side pagination + chip facets.

v2 list shape: {rows: [...], total: int}. Pagination via limit/offset
(0-indexed). Chip facets: family, k_min, k_max, primitive (yes/no),
has_search_results.
"""

from __future__ import annotations

from fastapi import APIRouter, Query, Request
from pydantic import BaseModel

from regpoly_web.database import json_loads
from regpoly_web.param_format import format_gen_params

router = APIRouter()


class _GeneratorRow(BaseModel):
    id: int
    search_run_id: int | None
    family: str
    L: int
    k: int
    found_at_try: int | None
    char_poly: str | None
    hamming_weight: int | None
    library_id: str | None
    structural_params: dict | None
    all_params: dict | None


class GeneratorList(BaseModel):
    rows: list[_GeneratorRow]
    total: int


def _row_to_v2(row) -> _GeneratorRow:
    keys = row.keys() if hasattr(row, "keys") else []
    family = row["family"]

    def opt(name):
        return row[name] if name in keys else None

    return _GeneratorRow(
        id=row["id"],
        search_run_id=row["search_run_id"],
        family=family,
        L=row["L"],
        k=row["k"],
        found_at_try=opt("found_at_try"),
        char_poly=opt("char_poly"),
        hamming_weight=opt("hamming_weight"),
        library_id=opt("library_id"),
        structural_params=format_gen_params(
            family, json_loads(row["structural_params"])
        ),
        all_params=format_gen_params(
            family, json_loads(row["all_params"])
        ),
    )


@router.get("/generators", response_model=GeneratorList)
async def v2_list_generators(
    request: Request,
    family: str | None = None,
    k: int | None = None,
    k_min: int | None = None,
    k_max: int | None = None,
    search_run_id: int | None = None,
    has_search_results: bool | None = None,
    primitive: str | None = Query(
        default=None, pattern="^(yes|no)$",
        description="primitive yes|no",
    ),
    limit: int = Query(50, ge=1, le=500),
    offset: int = Query(0, ge=0),
) -> GeneratorList:
    db = request.app.state.db

    where: list[str] = []
    params: list = []
    if family:
        where.append("family = ?")
        params.append(family)
    if k is not None:
        where.append("k = ?")
        params.append(k)
    if k_min is not None:
        where.append("k >= ?")
        params.append(k_min)
    if k_max is not None:
        where.append("k <= ?")
        params.append(k_max)
    if search_run_id is not None:
        where.append("search_run_id = ?")
        params.append(search_run_id)
    # has_search_results=true: only generators attached to a search run
    # (i.e. one that found at least them); =false: orphan rows.
    if has_search_results is True:
        where.append("search_run_id IS NOT NULL")
    elif has_search_results is False:
        where.append("search_run_id IS NULL")
    # primitive yes|no — primitive generators have a non-null,
    # non-error PIS analysis with a valid char poly.
    if primitive == "yes":
        where.append(
            "pis_computed_at IS NOT NULL AND pis_error IS NULL "
            "AND char_poly IS NOT NULL"
        )
    elif primitive == "no":
        where.append("(pis_error IS NOT NULL OR char_poly IS NULL)")

    where_clause = ("WHERE " + " AND ".join(where)) if where else ""

    async with db.execute(
        f"SELECT COUNT(*) FROM primitive_generator {where_clause}", params,
    ) as cur:
        total = (await cur.fetchone())[0]

    list_sql = (
        f"SELECT * FROM primitive_generator {where_clause} "
        f"ORDER BY id DESC LIMIT ? OFFSET ?"
    )
    async with db.execute(list_sql, [*params, limit, offset]) as cur:
        rows = await cur.fetchall()

    return GeneratorList(
        rows=[_row_to_v2(r) for r in rows],
        total=total,
    )


class _FamilyCounts(BaseModel):
    """Plain dict-style — Pydantic supports root-typing via __root__
    in v1 but in v2 we use a free-form dict response."""


@router.get("/generators/families/counts")
async def v2_generators_families_counts(request: Request) -> dict[str, int]:
    db = request.app.state.db
    async with db.execute(
        "SELECT family, COUNT(*) FROM primitive_generator "
        "WHERE pis_computed_at IS NOT NULL AND pis_error IS NULL "
        "AND char_poly IS NOT NULL "
        "GROUP BY family ORDER BY family"
    ) as cur:
        rows = await cur.fetchall()
    return {r[0]: int(r[1]) for r in rows}
