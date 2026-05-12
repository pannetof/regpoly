# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Web-layer tests for the primitive-search SE post-filter.

Exercise the route, model and DB-column wiring for the ``max_se`` knob
without spinning the worker pool: POST a search, read it back, and
confirm the new column round-trips. Restart-clone preserves ``max_se``
and resets ``rejected_count``. Worker-side filter behaviour is covered
by ``test_tasks_primitive_se_filter``.
"""

from __future__ import annotations


def _tausworthe_body(**overrides) -> dict:
    body = {
        "family": "TauswortheGen",
        "L": 0,
        "structural_params": {"k": 31, "s": 3, "q": 13},
        "fixed_params": {},
        "max_tries": 1,
    }
    body.update(overrides)
    return body


def test_max_se_persists_round_trip(client) -> None:
    body = _tausworthe_body(max_se=42)
    r = client.post("/api/primitive-searches", json=body)
    assert r.status_code == 200, r.text
    run_id = r.json()["id"]

    detail = client.get(f"/api/primitive-searches/{run_id}").json()
    assert detail["max_se"] == 42
    assert detail["rejected_count"] == 0


def test_no_max_se_round_trips_as_null(client) -> None:
    """Omitting max_se leaves it NULL — the search runs unfiltered."""
    r = client.post("/api/primitive-searches", json=_tausworthe_body())
    assert r.status_code == 200, r.text
    run_id = r.json()["id"]

    detail = client.get(f"/api/primitive-searches/{run_id}").json()
    assert detail["max_se"] is None
    assert detail["rejected_count"] == 0


def test_max_se_negative_rejected(client) -> None:
    """Pydantic ge=0 must reject negative thresholds."""
    r = client.post("/api/primitive-searches",
                    json=_tausworthe_body(max_se=-1))
    assert r.status_code == 422, r.text


def test_max_se_zero_is_valid(client) -> None:
    """max_se=0 means ME-perfect only; a valid, distinct configuration
    from "no filter"."""
    r = client.post("/api/primitive-searches",
                    json=_tausworthe_body(max_se=0))
    assert r.status_code == 200, r.text
    run_id = r.json()["id"]
    detail = client.get(f"/api/primitive-searches/{run_id}").json()
    assert detail["max_se"] == 0


def test_restart_clones_max_se_and_resets_rejected(client) -> None:
    """Restart must preserve ``max_se`` but ``rejected_count`` starts
    fresh on the new run."""
    body = _tausworthe_body(max_se=7)
    r = client.post("/api/primitive-searches", json=body)
    assert r.status_code == 200, r.text
    orig = r.json()["id"]

    r = client.post(f"/api/primitive-searches/{orig}/restart")
    assert r.status_code == 200, r.text
    new_id = r.json()["id"]
    assert new_id != orig

    detail = client.get(f"/api/primitive-searches/{new_id}").json()
    assert detail["max_se"] == 7
    assert detail["rejected_count"] == 0
