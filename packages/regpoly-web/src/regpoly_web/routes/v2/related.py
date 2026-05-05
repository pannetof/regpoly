"""GET /api/v2/generators/{id}/related — cross-references for the
generator detail page.

Returns:
  - primitive_search_id: the run that found this generator (None for
    orphan rows imported from YAML).
  - tempering_runs: tempering_search_run ids that pulled this generator
    into a component pool. Reverse-lookup via tempering_search_generator
    → tempering_search_component → tempering_search_run.
  - library_id: the library entry this generator points at, if any.
"""

from __future__ import annotations

from fastapi import APIRouter, HTTPException, Request
from pydantic import BaseModel

router = APIRouter()


class _RelatedResponse(BaseModel):
    primitive_search_id: int | None
    tempering_runs: list[int]
    library_id: str | None


@router.get(
    "/generators/{gen_id}/related",
    response_model=_RelatedResponse,
)
async def v2_generator_related(request: Request, gen_id: int) -> _RelatedResponse:
    db = request.app.state.db
    async with db.execute(
        "SELECT id, search_run_id, library_id FROM primitive_generator "
        "WHERE id = ?", (gen_id,),
    ) as cur:
        row = await cur.fetchone()
    if row is None:
        raise HTTPException(404, f"Generator {gen_id} not found")

    async with db.execute(
        "SELECT DISTINCT tsc.search_run_id "
        "FROM tempering_search_generator AS tsg "
        "JOIN tempering_search_component AS tsc "
        "  ON tsc.id = tsg.component_id "
        "WHERE tsg.generator_id = ?",
        (gen_id,),
    ) as cur:
        runs = [int(r["search_run_id"]) for r in await cur.fetchall()]

    return _RelatedResponse(
        primitive_search_id=row["search_run_id"],
        tempering_runs=runs,
        library_id=row["library_id"] if "library_id" in row.keys() else None,
    )
