# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Web-layer tests for the WELL cost-bounded primitive search.

These tests exercise the route / model / DB-column wiring without
spinning the worker pool: they POST a search, snapshot the persisted
row, and confirm the validation rules and the new ``max_cost`` column
round-trip end to end. Worker-side cost enforcement is covered by the
C++ test ``WellMaxCostBoundsAllHits`` and by the Python sampler tests.
"""

from __future__ import annotations


def _well_body(**overrides) -> dict:
    body = {
        "family": "WELLRNG",
        "L": 32,
        "structural_params": {"w": 32, "r": 16, "p": 0},
        "fixed_params": {"m1": 13, "m2": 9, "m3": 5},
        "max_tries": 50,
    }
    body.update(overrides)
    return body


def test_max_cost_persists_round_trip(client) -> None:
    body = _well_body(max_cost=12)
    r = client.post("/api/primitive-searches", json=body)
    assert r.status_code == 200, r.text
    run_id = r.json()["id"]

    r = client.get(f"/api/primitive-searches/{run_id}")
    assert r.status_code == 200
    assert r.json()["max_cost"] == 12


def test_well_search_pins_L_to_w(client) -> None:
    """WELL outputs one full w-bit word per step, so the search must
    always run at L=w regardless of the (sometimes huge) state size k.
    Previously the route defaulted to L=k for non-Tausworthe families."""
    body = _well_body(max_cost=12)
    body["L"] = 0   # no explicit override → server picks
    r = client.post("/api/primitive-searches", json=body)
    assert r.status_code == 200, r.text
    run_id = r.json()["id"]

    detail = client.get(f"/api/primitive-searches/{run_id}").json()
    assert detail["L"] == 32, detail   # body.structural_params['w'] == 32
    # k for (w=32, r=16, p=0) is 32*16 = 512; L must NOT be k anymore.
    assert detail["k"] == 512, detail
    assert detail["L"] != detail["k"]


def test_max_cost_out_of_range_rejected(client) -> None:
    r = client.post("/api/primitive-searches", json=_well_body(max_cost=99))
    assert r.status_code == 400
    detail = r.json()["detail"]
    assert detail["code"] == "max_cost_out_of_range"


def test_max_cost_non_well_family_rejected(client) -> None:
    body = {
        "family": "TGFSRGen",
        "L": 32,
        "structural_params": {"w": 32, "r": 3},
        "fixed_params": {"m": 1},
        "max_tries": 5,
        "max_cost": 12,
    }
    r = client.post("/api/primitive-searches", json=body)
    assert r.status_code == 400
    assert r.json()["detail"]["code"] == "max_cost_well_only"


def test_pinned_matrices_with_max_cost_rejected(client) -> None:
    pinned = {
        "T0": {"M": 1}, "T1": {"M": 1}, "T2": {"M": 1}, "T3": {"M": 1},
        "T4": {"M": 1}, "T5": {"M": 1}, "T6": {"M": 1}, "T7": {"M": 1},
    }
    body = _well_body(
        max_cost=12,
        fixed_params={"m1": 13, "m2": 9, "m3": 5, "matrices": pinned},
    )
    r = client.post("/api/primitive-searches", json=body)
    assert r.status_code == 400
    assert r.json()["detail"]["code"] == "max_cost_with_pinned_matrices"


def test_no_max_cost_is_allowed(client) -> None:
    """Backwards-compat: omitting max_cost yields a row with NULL.
    Without a cap, WELL needs `matrices` pinned for the probe."""
    pinned = {f"T{i}": {"M": 0} for i in range(8)}
    body = _well_body(
        fixed_params={"m1": 13, "m2": 9, "m3": 5, "matrices": pinned},
    )
    r = client.post("/api/primitive-searches", json=body)
    assert r.status_code == 200, r.text
    run_id = r.json()["id"]
    detail = client.get(f"/api/primitive-searches/{run_id}").json()
    assert detail["max_cost"] is None
