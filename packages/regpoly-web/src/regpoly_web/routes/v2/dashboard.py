"""GET /api/v2/dashboard/summary — single-fetch dashboard payload.

Replaces the four-fetch fan-out the legacy `dashboard()` Alpine
component did. Counts + active runs + recent runs + recent papers,
each active row carrying rate / eta_seconds / 30-point sparkline.
"""

from __future__ import annotations

from collections import OrderedDict
from typing import Literal

from fastapi import APIRouter, Request
from pydantic import BaseModel

router = APIRouter()


class _Counts(BaseModel):
    generators_tested: int
    primitive_searches_total: int
    tempering_searches_total: int
    finds_last_24h: int


class _RunRow(BaseModel):
    id: int
    type: Literal["primitive", "tempering"]
    family: str | None
    k: int | None
    status: str
    started_at: str | None
    finished_at: str | None
    rate: float | None
    eta_seconds: float | None
    sparkline: list[int]


class _PaperRow(BaseModel):
    id: str
    display: str


class DashboardSummary(BaseModel):
    counts: _Counts
    active: list[_RunRow]
    recent: list[_RunRow]
    recent_papers: list[_PaperRow]


def _format_iso(dt) -> str | None:
    if dt is None:
        return None
    if hasattr(dt, "isoformat"):
        return dt.isoformat()
    return str(dt)


async def _sparkline_for_run(db, run_id: int, run_type: str, points: int = 30) -> list[int]:
    """Down-sampled sparkline from search_progress rows for this run.

    Uses index-bucket down-sampling: pick `points` evenly-spaced rows
    from the N progress rows. Guarantees exactly `min(N, points)`
    output samples regardless of N (no modulo off-by-one).
    """
    sql = (
        "SELECT id, tries_done, found_count "
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
    return [int(rows[i * n // points]["found_count"] or 0) for i in range(points)]


@router.get("/dashboard/summary", response_model=DashboardSummary)
async def dashboard_summary(request: Request) -> DashboardSummary:
    db = request.app.state.db

    # ----- counts -----
    async with db.execute(
        "SELECT COUNT(*) FROM tested_generator"
    ) as cur:
        generators_tested = (await cur.fetchone())[0]
    async with db.execute(
        "SELECT COUNT(*) FROM primitive_search_run"
    ) as cur:
        primitive_total = (await cur.fetchone())[0]
    async with db.execute(
        "SELECT COUNT(*) FROM tempering_search_run"
    ) as cur:
        tempering_total = (await cur.fetchone())[0]
    async with db.execute(
        "SELECT COUNT(*) FROM primitive_generator "
        "WHERE created_at >= datetime('now', '-1 day')"
    ) as cur:
        finds_24h = (await cur.fetchone())[0]

    # ----- active runs (running + pending) -----
    active: list[_RunRow] = []
    async with db.execute(
        "SELECT id, family, k, status, started_at, "
        "       NULL AS finished_at "
        "FROM primitive_search_run "
        "WHERE status IN ('running', 'pending') "
        "ORDER BY started_at DESC LIMIT 20"
    ) as cur:
        for r in await cur.fetchall():
            sparkline = await _sparkline_for_run(db, r["id"], "primitive")
            active.append(_RunRow(
                id=r["id"], type="primitive",
                family=r["family"], k=r["k"], status=r["status"],
                started_at=_format_iso(r["started_at"]),
                finished_at=None,
                rate=None, eta_seconds=None,
                sparkline=sparkline,
            ))
    async with db.execute(
        "SELECT id, NULL AS family, NULL AS k, status, started_at "
        "FROM tempering_search_run "
        "WHERE status IN ('running', 'pending') "
        "ORDER BY started_at DESC LIMIT 20"
    ) as cur:
        for r in await cur.fetchall():
            sparkline = await _sparkline_for_run(db, r["id"], "tempering")
            active.append(_RunRow(
                id=r["id"], type="tempering",
                family=None, k=None, status=r["status"],
                started_at=_format_iso(r["started_at"]),
                finished_at=None,
                rate=None, eta_seconds=None,
                sparkline=sparkline,
            ))

    # ----- recent (closed) runs, capped at 10 -----
    recent: list[_RunRow] = []
    async with db.execute(
        "SELECT id, family, k, status, started_at, finished_at "
        "FROM primitive_search_run "
        "WHERE status IN ('completed', 'cancelled', 'failed') "
        "ORDER BY COALESCE(finished_at, started_at) DESC LIMIT 10"
    ) as cur:
        for r in await cur.fetchall():
            sparkline = await _sparkline_for_run(db, r["id"], "primitive")
            recent.append(_RunRow(
                id=r["id"], type="primitive",
                family=r["family"], k=r["k"], status=r["status"],
                started_at=_format_iso(r["started_at"]),
                finished_at=_format_iso(r["finished_at"]),
                rate=None, eta_seconds=None,
                sparkline=sparkline,
            ))
    async with db.execute(
        "SELECT id, NULL AS family, NULL AS k, status, started_at, "
        "       finished_at "
        "FROM tempering_search_run "
        "WHERE status IN ('completed', 'cancelled', 'failed') "
        "ORDER BY COALESCE(finished_at, started_at) DESC LIMIT 10"
    ) as cur:
        for r in await cur.fetchall():
            sparkline = await _sparkline_for_run(db, r["id"], "tempering")
            recent.append(_RunRow(
                id=r["id"], type="tempering",
                family=None, k=None, status=r["status"],
                started_at=_format_iso(r["started_at"]),
                finished_at=_format_iso(r["finished_at"]),
                rate=None, eta_seconds=None,
                sparkline=sparkline,
            ))
    # Sort merged recent by finished_at desc, then trim.
    recent.sort(
        key=lambda r: r.finished_at or r.started_at or "",
        reverse=True,
    )
    recent = recent[:10]

    # ----- recent papers (catalog) -----
    catalog = request.app.state.library
    papers = [
        _PaperRow(id=p.id, display=p.display())
        for p in catalog.papers()
    ][:6]

    return DashboardSummary(
        counts=_Counts(
            generators_tested=generators_tested,
            primitive_searches_total=primitive_total,
            tempering_searches_total=tempering_total,
            finds_last_24h=finds_24h,
        ),
        active=active,
        recent=recent,
        recent_papers=papers,
    )
