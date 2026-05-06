# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""GET /api/v2/dashboard/summary — single-fetch dashboard payload.

Replaces the four-fetch fan-out the legacy `dashboard()` Alpine
component did. Counts + active runs + recent runs + recent papers,
each active row carrying rate / eta_seconds / 30-point sparkline.
"""

from __future__ import annotations

from collections import OrderedDict, defaultdict
from typing import Literal

from fastapi import APIRouter, Request
from pydantic import BaseModel

from regpoly_web.database import json_loads
from regpoly_web.routes.v2._downsample import bucket_downsample

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


async def _sparklines_for_runs(
    db, run_ids: list[int], run_type: str, points: int = 30,
) -> dict[int, list[int]]:
    """Batched sparkline lookup — one SQL roundtrip for every active /
    recent row. Returns {run_id: downsampled tries_done series}."""
    if not run_ids:
        return {}
    placeholders = ",".join("?" for _ in run_ids)
    async with db.execute(
        f"SELECT id, search_run_id, tries_done FROM search_progress "
        f"WHERE search_type = ? AND search_run_id IN ({placeholders}) "
        f"ORDER BY id ASC",
        [run_type, *run_ids],
    ) as cur:
        rows = await cur.fetchall()
    by_run: dict[int, list[int]] = defaultdict(list)
    for r in rows:
        by_run[int(r["search_run_id"])].append(int(r["tries_done"] or 0))
    return {rid: bucket_downsample(by_run.get(rid, []), points) for rid in run_ids}


async def _live_metrics_for_run(
    db, run_id: int, run_type: str,
) -> tuple[float | None, float | None]:
    """Read rate / eta_seconds for an active run from the latest progress
    row's `current_info` JSON (and fall back to the in-memory rolling
    rate snapshot when present)."""
    rate: float | None = None
    eta_seconds: float | None = None
    async with db.execute(
        "SELECT current_info FROM search_progress "
        "WHERE search_type = ? AND search_run_id = ? "
        "ORDER BY id DESC LIMIT 1",
        [run_type, run_id],
    ) as cur:
        row = await cur.fetchone()
    if row is not None and row["current_info"]:
        info = json_loads(row["current_info"]) or {}
        # Prefer rolling rate (smoother); fall back to instantaneous.
        for key in ("rate_rolling_5s", "rate"):
            v = info.get(key)
            if v is not None:
                try:
                    rate = float(v)
                    break
                except (TypeError, ValueError):
                    pass
        v = info.get("eta_seconds")
        if v is not None:
            try:
                eta_seconds = float(v)
            except (TypeError, ValueError):
                pass

    if rate is None:
        try:
            from regpoly_web.tasks._progress_rate import snapshot
            snap = snapshot(run_id)
            if snap is not None:
                rate = snap.get("rate")
                if eta_seconds is None:
                    eta_seconds = snap.get("eta_seconds")
        except Exception:
            pass

    return rate, eta_seconds


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
        active_prim_rows = await cur.fetchall()
    async with db.execute(
        "SELECT id, NULL AS family, NULL AS k, status, started_at "
        "FROM tempering_search_run "
        "WHERE status IN ('running', 'pending') "
        "ORDER BY started_at DESC LIMIT 20"
    ) as cur:
        active_temp_rows = await cur.fetchall()

    prim_sparks = await _sparklines_for_runs(
        db, [r["id"] for r in active_prim_rows], "primitive",
    )
    temp_sparks = await _sparklines_for_runs(
        db, [r["id"] for r in active_temp_rows], "tempering",
    )
    for r in active_prim_rows:
        rate, eta = await _live_metrics_for_run(db, r["id"], "primitive")
        active.append(_RunRow(
            id=r["id"], type="primitive",
            family=r["family"], k=r["k"], status=r["status"],
            started_at=_format_iso(r["started_at"]),
            finished_at=None,
            rate=rate, eta_seconds=eta,
            sparkline=prim_sparks.get(r["id"], []),
        ))
    for r in active_temp_rows:
        rate, eta = await _live_metrics_for_run(db, r["id"], "tempering")
        active.append(_RunRow(
            id=r["id"], type="tempering",
            family=None, k=None, status=r["status"],
            started_at=_format_iso(r["started_at"]),
            finished_at=None,
            rate=rate, eta_seconds=eta,
            sparkline=temp_sparks.get(r["id"], []),
        ))

    # ----- recent (closed) runs, capped at 10 -----
    recent: list[_RunRow] = []
    async with db.execute(
        "SELECT id, family, k, status, started_at, finished_at "
        "FROM primitive_search_run "
        "WHERE status IN ('completed', 'cancelled', 'failed') "
        "ORDER BY COALESCE(finished_at, started_at) DESC LIMIT 10"
    ) as cur:
        recent_prim = await cur.fetchall()
    async with db.execute(
        "SELECT id, NULL AS family, NULL AS k, status, started_at, "
        "       finished_at "
        "FROM tempering_search_run "
        "WHERE status IN ('completed', 'cancelled', 'failed') "
        "ORDER BY COALESCE(finished_at, started_at) DESC LIMIT 10"
    ) as cur:
        recent_temp = await cur.fetchall()
    rprim_sparks = await _sparklines_for_runs(
        db, [r["id"] for r in recent_prim], "primitive",
    )
    rtemp_sparks = await _sparklines_for_runs(
        db, [r["id"] for r in recent_temp], "tempering",
    )
    for r in recent_prim:
        recent.append(_RunRow(
            id=r["id"], type="primitive",
            family=r["family"], k=r["k"], status=r["status"],
            started_at=_format_iso(r["started_at"]),
            finished_at=_format_iso(r["finished_at"]),
            rate=None, eta_seconds=None,
            sparkline=rprim_sparks.get(r["id"], []),
        ))
    for r in recent_temp:
        recent.append(_RunRow(
            id=r["id"], type="tempering",
            family=None, k=None, status=r["status"],
            started_at=_format_iso(r["started_at"]),
            finished_at=_format_iso(r["finished_at"]),
            rate=None, eta_seconds=None,
            sparkline=rtemp_sparks.get(r["id"], []),
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
