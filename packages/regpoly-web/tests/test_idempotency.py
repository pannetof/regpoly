"""Phase 3 red — POST lifecycle endpoints are idempotent.

Calling cancel/pause/resume/restart twice never returns 409. Contract
test for both v1 and v2 endpoints.
"""

from __future__ import annotations


def _create_psr(client) -> int:
    r = client.post(
        "/api/primitive-searches",
        json={
            "family": "MTGen",
            "L": 64,
            "structural_params": {"w": 32, "n": 2, "m": 1, "r": 31, "u": 11},
            "fixed_params": {},
            "max_tries": 1,
            "max_seconds": None,
        },
    )
    assert r.status_code in (200, 201), r.text
    return r.json()["id"]


def test_cancel_twice_returns_200(client) -> None:
    pid = _create_psr(client)
    a = client.post(f"/api/primitive-searches/{pid}/cancel")
    b = client.post(f"/api/primitive-searches/{pid}/cancel")
    assert a.status_code == 200, a.text
    assert b.status_code == 200, b.text


def test_pause_already_paused_returns_200(client) -> None:
    pid = _create_psr(client)
    a = client.post(f"/api/primitive-searches/{pid}/pause")
    b = client.post(f"/api/primitive-searches/{pid}/pause")
    assert a.status_code == 200, a.text
    assert b.status_code == 200, b.text
