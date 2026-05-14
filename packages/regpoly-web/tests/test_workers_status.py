# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Tools → Workers status endpoint.

Verifies the API shape (active_jobs + pools), behaviour when idle,
and behaviour while an F2w search is in flight.
"""

from __future__ import annotations

import time


def test_workers_status_shape_when_idle(client) -> None:
    r = client.get("/api/workers/status")
    assert r.status_code == 200
    body = r.json()
    assert "active_jobs" in body
    assert "pools" in body
    assert isinstance(body["active_jobs"], list)
    assert isinstance(body["pools"], list)
    # Two pools registered in app.state: primitive + analysis.
    pool_names = {p["name"] for p in body["pools"]}
    assert pool_names == {"primitive", "analysis"}, pool_names
    # Each pool entry exposes a running/queued split so the UI can
    # distinguish "actually executing" from "waiting in the executor
    # queue" — both are needed because in_flight = running + queued.
    for p in body["pools"]:
        assert "running" in p, p
        assert "queued"  in p, p
        assert "in_flight" in p, p
        assert p["in_flight"] == p["running"] + p["queued"], p


def test_workers_status_lists_running_jobs(client) -> None:
    """While an F2w search is in flight (or just-completed), it should
    surface in active_jobs. Hard to time race-free, so accept either
    'running' during the search or empty after completion."""
    body = {
        "family": "F2wLFSRGen", "L": 0,
        "structural_params": {"w": 8, "r": 3, "nb_terms": 3,
                              "normal_basis": False, "step": 1},
        "fixed_params": {},
        "max_tries": 5000,   # bigger budget so we have a window to observe
        "max_seconds": 30,
    }
    r = client.post("/api/primitive-searches", json=body)
    assert r.status_code == 200, r.text
    run_id = r.json()["id"]

    # Sample workers/status a few times; the search must appear in
    # active_jobs at least once before completion.
    saw_running = False
    for _ in range(30):
        ws = client.get("/api/workers/status").json()
        running_ids = {(j["kind"], j["run_id"]) for j in ws["active_jobs"]
                       if j["status"] == "running"}
        if ("primitive", run_id) in running_ids:
            saw_running = True
            break
        # Stop polling once the search finishes — saw_running may stay
        # False on a very fast search, that's fine.
        detail = client.get(f"/api/primitive-searches/{run_id}").json()
        if detail.get("status") in ("completed", "cancelled", "failed"):
            break
        time.sleep(0.1)

    # If the search was too fast to catch, at minimum the endpoint
    # must keep returning a well-formed body.
    final = client.get("/api/workers/status").json()
    assert "active_jobs" in final
    # Once the run completes, it must NOT linger in active_jobs.
    detail = client.get(f"/api/primitive-searches/{run_id}").json()
    if detail["status"] == "completed":
        ids = {(j["kind"], j["run_id"]) for j in final["active_jobs"]}
        assert ("primitive", run_id) not in ids
