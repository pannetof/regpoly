"""GET /api/v2/searches/histories — batch sparkline payload.

Per Researcher persona R-9: avoid the N+1 fetches on /searches by
loading every visible row's history in a single round-trip.

  ?ids=1,2,3&type=primitive&points=30
  → {histories: {id: {points: [int...], status: str, type: str}}}
"""

from __future__ import annotations

from typing import Literal

from fastapi import APIRouter, Query, Request
from pydantic import BaseModel

router = APIRouter()


class _History(BaseModel):
    points: list[int]
    status: str
    type: Literal["primitive", "tempering"]


class HistoriesEnvelope(BaseModel):
    histories: dict[str, _History]


async def _history_points_for_run(
    db, run_id: int, run_type: str, points: int,
) -> list[int]:
    """Index-bucket down-sampling — exact `min(N, points)` samples."""
    sql = (
        "SELECT id, found_count "
        "FROM search_progress "
        "WHERE search_type = ? AND search_run_id = ? "
        "ORDER BY id ASC"
    )
    async with db.execute(sql, [run_type, run_id]) as cur:
        rows = await cur.fetchall()
    if not rows:
        return []
    n = len(rows)
    if n <= points:
        return [int(r["found_count"] or 0) for r in rows]
    return [int(rows[i * n // points]["found_count"] or 0)
            for i in range(points)]


def _parse_ids(raw: str | None) -> list[int]:
    if not raw:
        return []
    out: list[int] = []
    for tok in raw.split(","):
        tok = tok.strip()
        if not tok:
            continue
        try:
            out.append(int(tok))
        except ValueError:
            continue
    return out


@router.get("/searches/histories", response_model=HistoriesEnvelope)
async def v2_searches_histories(
    request: Request,
    ids: str = Query(default="", description="comma-separated run ids"),
    type: Literal["primitive", "tempering"] = "primitive",
    points: int = Query(default=30, ge=1, le=300),
) -> HistoriesEnvelope:
    db = request.app.state.db
    id_list = _parse_ids(ids)
    histories: dict[str, _History] = {}
    if not id_list:
        return HistoriesEnvelope(histories=histories)
    table = (
        "primitive_search_run" if type == "primitive"
        else "tempering_search_run"
    )
    placeholders = ",".join("?" for _ in id_list)
    async with db.execute(
        f"SELECT id, status FROM {table} WHERE id IN ({placeholders})",
        id_list,
    ) as cur:
        meta = {r["id"]: r["status"] for r in await cur.fetchall()}
    for rid in id_list:
        if rid not in meta:
            continue
        pts = await _history_points_for_run(db, rid, type, points)
        histories[str(rid)] = _History(
            points=pts, status=meta[rid], type=type,
        )
    return HistoriesEnvelope(histories=histories)
