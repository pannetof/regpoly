# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Worker-status endpoint backing the Tools → Workers tab.

Combines two signals so the page works in both the host-dev and the
in-stack-compose deployments:

- DB-side "active jobs" — every row with status='running' or 'paused'
  across the three task tables (primitive_search_run,
  tempering_search_run, library_test_run). This is the source of
  truth regardless of which container's pool is running them.

- Process-side pool snapshot — `app.state.pool` /
  `app.state.analysis_pool` `_futures` mapping. Populated when the
  web container's internal pool is enabled (the host-dev mode);
  empty when `REGPOLY_WEB_DISABLE_INTERNAL_POOL=1` and the standalone
  worker container owns the dispatch.
"""

from __future__ import annotations

import time

from fastapi import APIRouter, Request

router = APIRouter()


@router.get("/workers/status")
async def workers_status(request: Request) -> dict:
    db = request.app.state.db

    # ── DB-side: active jobs across all task tables ─────────────────
    # primitive_search_run / tempering_search_run carry an explicit
    # `started_at`; library_test_run only has `created_at` so we fall
    # back to that as the elapsed-time anchor.
    active = []
    for kind, table, ts_col in (
        ("primitive",    "primitive_search_run",  "started_at"),
        ("tempering",    "tempering_search_run",  "started_at"),
        ("library_test", "library_test_run",      "created_at"),
    ):
        async with db.execute(
            f"SELECT id, status, {ts_col}, "
            f"       EXTRACT(EPOCH FROM (NOW() - {ts_col})) "
            f"FROM {table} "
            f"WHERE status IN ('running', 'paused') "
            f"ORDER BY id",
        ) as cur:
            for row in await cur.fetchall():
                started = row[2]
                active.append({
                    "kind": kind,
                    "run_id": int(row[0]),
                    "status": row[1],
                    "started_at": (started.isoformat()
                                   if hasattr(started, "isoformat")
                                   else (str(started) if started else None)),
                    "elapsed_seconds": (float(row[3])
                                        if row[3] is not None else None),
                })

    # ── Process-side: this web container's pool snapshot ────────────
    pools = []
    for name, attr in (("primitive",  "pool"),
                       ("analysis",   "analysis_pool")):
        pool_obj = getattr(request.app.state, attr, None)
        if pool_obj is None:
            pools.append({"name": name, "max_workers": None,
                          "in_flight": 0,
                          "note": "disabled in this container "
                                  "(REGPOLY_WEB_DISABLE_INTERNAL_POOL=1)"})
            continue
        executor = getattr(pool_obj, "executor", None)
        max_workers = (getattr(executor, "_max_workers", None)
                       if executor is not None else None)
        # `_futures` holds every submitted task until its done-callback
        # pops it. Split into actually-executing (running()) and queued
        # so the UI doesn't mislead users into thinking the pool is
        # stuck: with max_workers=2, only 2 can ever truly run; the
        # rest sit in the executor's submit queue waiting their turn.
        running, queued = [], []
        now = time.monotonic()
        submitted_at = getattr(pool_obj, "_submitted_at", {})
        for (k, rid), fut in pool_obj._futures.items():
            if fut.done():
                continue
            t0 = submitted_at.get((k, rid))
            elapsed = round(now - t0, 1) if t0 is not None else None
            entry = {"kind": k, "run_id": rid, "elapsed_seconds": elapsed}
            if fut.running():
                running.append(entry)
            else:
                queued.append(entry)
        pools.append({
            "name": name,
            "max_workers": max_workers,
            "running": len(running),
            "queued": len(queued),
            "in_flight": len(running) + len(queued),
            "tasks": running + queued,
            "running_tasks": running,
            "queued_tasks": queued,
        })

    return {"active_jobs": active, "pools": pools}
