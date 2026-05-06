# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""GET /api/v2/{primitive,tempering}-searches/{id}/history — index-bucket
downsampling that returns *exactly* K points (or fewer if the run hasn't
emitted that much progress yet).

Algorithm:
  1) read all (id, current_info) rows from search_progress for the run
  2) extract the metric (`tries_done` for primitive, `best_se` for
     tempering)
  3) if len(samples) <= K: return as-is
  4) else: bucket the array into K groups by index, take last value
     per bucket
"""

from __future__ import annotations

import json
from typing import Literal

from fastapi import APIRouter, HTTPException, Query, Request
from pydantic import BaseModel

from regpoly_web.database import json_loads
from regpoly_web.routes.v2._downsample import bucket_downsample

router = APIRouter()


class HistoryResponse(BaseModel):
    points: list[float]
    type: Literal["primitive", "tempering"]
    status: str
    best_id: int | None = None


def _bucket_downsample(values: list[float], k: int) -> list[float]:
    """Backwards-compatible alias for the shared helper. Kept so
    third-party code that imports this name doesn't break."""
    return bucket_downsample(values, k)


@router.get(
    "/primitive-searches/{run_id}/history",
    response_model=HistoryResponse,
)
async def v2_primitive_history(
    request: Request, run_id: int,
    points: int = Query(50, ge=1, le=1000),
) -> HistoryResponse:
    db = request.app.state.db
    async with db.execute(
        "SELECT status FROM primitive_search_run WHERE id = ?", (run_id,),
    ) as cur:
        row = await cur.fetchone()
    if row is None:
        raise HTTPException(404, f"Search {run_id} not found")
    status = row["status"]

    async with db.execute(
        "SELECT tries_done FROM search_progress "
        "WHERE search_type='primitive' AND search_run_id = ? "
        "ORDER BY id ASC",
        (run_id,),
    ) as cur:
        rows = await cur.fetchall()
    series = [float(r["tries_done"]) for r in rows]
    return HistoryResponse(
        points=_bucket_downsample(series, points),
        type="primitive",
        status=status,
    )


@router.get(
    "/tempering-searches/{run_id}/history",
    response_model=HistoryResponse,
)
async def v2_tempering_history(
    request: Request, run_id: int,
    points: int = Query(50, ge=1, le=1000),
) -> HistoryResponse:
    db = request.app.state.db
    async with db.execute(
        "SELECT status, best_se FROM tempering_search_run WHERE id = ?",
        (run_id,),
    ) as cur:
        row = await cur.fetchone()
    if row is None:
        raise HTTPException(404, f"Search {run_id} not found")
    status = row["status"]

    async with db.execute(
        "SELECT current_info FROM search_progress "
        "WHERE search_type='tempering' AND search_run_id = ? "
        "ORDER BY id ASC",
        (run_id,),
    ) as cur:
        rows = await cur.fetchall()

    # best_se monotonically decreases over the run; for downsampling we
    # extract whichever-key-is-present from the JSON blob.
    series: list[float] = []
    for r in rows:
        info = json_loads(r["current_info"]) or {}
        v = info.get("best_overall_se")
        if v is None:
            v = info.get("se")
        if v is None:
            continue
        try:
            series.append(float(v))
        except (TypeError, ValueError):
            continue

    # Resolve best_id only once the run has finished. This signals the
    # UI's "Open in catalog" affordance.
    best_id: int | None = None
    if status == "completed":
        async with db.execute(
            "SELECT id FROM tested_generator "
            "WHERE search_run_id = ? ORDER BY id DESC LIMIT 1",
            (run_id,),
        ) as cur:
            tg = await cur.fetchone()
        if tg is not None:
            best_id = int(tg["id"])

    return HistoryResponse(
        points=_bucket_downsample(series, points),
        type="tempering",
        status=status,
        best_id=best_id,
    )


# ── Duplicate (returns prefilled form params) ──────────────────────────


class _PrefillResponse(BaseModel):
    family: str
    L: int
    k: int | None = None
    structural_params: dict
    fixed_params: dict
    max_tries: int | None = None
    max_seconds: float | None = None
    cloned_from: int


@router.get(
    "/primitive-searches/{run_id}/duplicate",
    response_model=_PrefillResponse,
)
async def v2_primitive_duplicate(request: Request, run_id: int) -> _PrefillResponse:
    db = request.app.state.db
    async with db.execute(
        "SELECT family, L, k, structural_params, fixed_params, "
        "max_tries, max_seconds FROM primitive_search_run WHERE id = ?",
        (run_id,),
    ) as cur:
        row = await cur.fetchone()
    if row is None:
        raise HTTPException(404, f"Search {run_id} not found")

    return _PrefillResponse(
        family=row["family"],
        L=int(row["L"]),
        k=int(row["k"]) if row["k"] is not None else None,
        structural_params=json.loads(row["structural_params"] or "{}"),
        fixed_params=json.loads(row["fixed_params"] or "{}"),
        max_tries=row["max_tries"],
        max_seconds=row["max_seconds"],
        cloned_from=run_id,
    )


class _TemperingPrefill(BaseModel):
    cloned_from: int
    test_config: dict
    components: list[dict]
    nb_tries: int | None = None


@router.get(
    "/tempering-searches/{run_id}/duplicate",
    response_model=_TemperingPrefill,
)
async def v2_tempering_duplicate(request: Request, run_id: int) -> _TemperingPrefill:
    db = request.app.state.db
    async with db.execute(
        "SELECT test_config, nb_tries FROM tempering_search_run WHERE id = ?",
        (run_id,),
    ) as cur:
        row = await cur.fetchone()
    if row is None:
        raise HTTPException(404, f"Search {run_id} not found")

    async with db.execute(
        "SELECT component_index, shared_with_component, tempering_config "
        "FROM tempering_search_component "
        "WHERE search_run_id = ? ORDER BY component_index ASC",
        (run_id,),
    ) as cur:
        comps = await cur.fetchall()

    return _TemperingPrefill(
        cloned_from=run_id,
        test_config=json.loads(row["test_config"] or "{}"),
        components=[
            {
                "component_index": c["component_index"],
                "shared_with_component": c["shared_with_component"],
                "tempering_config": json.loads(
                    c["tempering_config"] or "{}"),
            }
            for c in comps
        ],
        nb_tries=row["nb_tries"],
    )
