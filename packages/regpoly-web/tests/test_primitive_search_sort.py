# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Server-side sort on /primitive-searches/{id}/generators.

Verifies (a) allow-list rejects unknown sort keys, (b) the order
returned matches the requested sort_by + sort_dir, (c) the default
behaviour (no params) stays `id DESC` for backwards compatibility
with clients that haven't migrated.
"""

from __future__ import annotations

import time


def _kick_off_f2w_search(client) -> int:
    """Launch a small F2w search and wait for completion, return run id."""
    body = {
        "family": "F2wLFSRGen", "L": 0,
        "structural_params": {"w": 8, "r": 3, "nb_terms": 3,
                              "normal_basis": False, "step": 1},
        "fixed_params": {},
        "max_tries": 300,
        "max_seconds": 30,
    }
    r = client.post("/api/primitive-searches", json=body)
    assert r.status_code == 200, r.text
    run_id = r.json()["id"]
    for _ in range(60):
        detail = client.get(f"/api/primitive-searches/{run_id}").json()
        if detail.get("status") in ("completed", "cancelled", "failed"):
            break
        time.sleep(0.2)
    assert detail.get("status") == "completed", detail
    assert detail.get("found_count", 0) > 0
    return run_id


def test_sort_default_is_id_desc(client) -> None:
    run_id = _kick_off_f2w_search(client)
    items = client.get(
        f"/api/primitive-searches/{run_id}/generators").json()["items"]
    assert len(items) >= 2
    ids = [g["id"] for g in items]
    assert ids == sorted(ids, reverse=True), ids


def test_sort_by_id_asc(client) -> None:
    run_id = _kick_off_f2w_search(client)
    items = client.get(
        f"/api/primitive-searches/{run_id}/generators"
        "?sort_by=id&sort_dir=asc").json()["items"]
    ids = [g["id"] for g in items]
    assert ids == sorted(ids), ids


def test_sort_by_k_works(client) -> None:
    """k is the same for every hit in this search (24), but the query
    must still complete successfully — exercises the order-by branch."""
    run_id = _kick_off_f2w_search(client)
    r = client.get(
        f"/api/primitive-searches/{run_id}/generators"
        "?sort_by=k&sort_dir=desc")
    assert r.status_code == 200


def test_sort_by_unknown_key_rejected(client) -> None:
    run_id = _kick_off_f2w_search(client)
    r = client.get(
        f"/api/primitive-searches/{run_id}/generators"
        "?sort_by=hamming_weight")
    assert r.status_code == 400


def test_sort_dir_validated(client) -> None:
    run_id = _kick_off_f2w_search(client)
    r = client.get(
        f"/api/primitive-searches/{run_id}/generators"
        "?sort_dir=upward")
    assert r.status_code == 400
