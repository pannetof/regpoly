"""POST/DELETE /api/v2/tested-generators/{id}/publish — cataloger workflow.

Idempotent: re-publishing the same library_id returns 200; deleting an
already-empty library_id returns 200. The endpoint manipulates the
`tested_generator.library_id` column.
"""

from __future__ import annotations

from fastapi import APIRouter, HTTPException, Request
from pydantic import BaseModel

router = APIRouter()


class _PublishBody(BaseModel):
    library_id: str


class _PublishResponse(BaseModel):
    tested_gen_id: int
    library_id: str | None


async def _row_or_404(db, tg_id: int) -> dict:
    async with db.execute(
        "SELECT id, library_id FROM tested_generator WHERE id = ?",
        (tg_id,),
    ) as cur:
        row = await cur.fetchone()
    if row is None:
        raise HTTPException(404, f"Tested generator {tg_id} not found")
    return {"id": row["id"], "library_id": row["library_id"]}


@router.post(
    "/tested-generators/{tg_id}/publish", response_model=_PublishResponse,
)
async def v2_publish(
    request: Request, tg_id: int, body: _PublishBody,
) -> _PublishResponse:
    db = request.app.state.db
    await _row_or_404(db, tg_id)
    await db.execute(
        "UPDATE tested_generator SET library_id = ? WHERE id = ?",
        (body.library_id, tg_id),
    )
    await db.commit()
    return _PublishResponse(
        tested_gen_id=tg_id, library_id=body.library_id,
    )


@router.delete(
    "/tested-generators/{tg_id}/publish", response_model=_PublishResponse,
)
async def v2_unpublish(request: Request, tg_id: int) -> _PublishResponse:
    db = request.app.state.db
    await _row_or_404(db, tg_id)
    await db.execute(
        "UPDATE tested_generator SET library_id = NULL WHERE id = ?",
        (tg_id,),
    )
    await db.commit()
    return _PublishResponse(tested_gen_id=tg_id, library_id=None)
