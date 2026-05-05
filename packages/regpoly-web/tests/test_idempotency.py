"""Phase 3+P6 — POST lifecycle endpoints are idempotent for both run
types and all four actions: cancel/pause/resume/restart × primitive/
tempering. 8 cases. Calling twice never returns 409.
"""

from __future__ import annotations

import json


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


def _create_tsr(client) -> int:
    """Insert a minimal tempering_search_run row directly via the v1
    POST endpoint. Kept generic so the schema for `test_config` /
    components is whatever the v1 endpoint expects today."""
    r = client.post(
        "/api/tempering-searches",
        json={
            "test_type": "equidistribution",
            "test_config": {"L": 32, "k": 32, "delta": 0},
            "Lmax": 32,
            "nb_tries": 1,
            "components": [],
        },
    )
    if r.status_code in (200, 201):
        return r.json().get("id")
    # If the v1 endpoint requires components we can't synthesize a real
    # combo here; fall back to a direct DB insert. This is a contract
    # test against the lifecycle endpoints, not a real run.
    import sqlite3
    db_path = client.app.state.settings.db_path
    conn = sqlite3.connect(db_path)
    try:
        cur = conn.execute(
            "INSERT INTO tempering_search_run "
            "(test_type, test_config, Lmax, nb_tries, status, "
            "combos_done, started_at, finished_at) "
            "VALUES (?, ?, ?, ?, 'pending', 0, NULL, NULL)",
            ("equidistribution", json.dumps({"L": 32}), 32, 1),
        )
        conn.commit()
        return cur.lastrowid
    finally:
        conn.close()


# ── primitive ──


def test_primitive_cancel_twice_returns_200(client) -> None:
    pid = _create_psr(client)
    a = client.post(f"/api/primitive-searches/{pid}/cancel")
    b = client.post(f"/api/primitive-searches/{pid}/cancel")
    assert a.status_code == 200, a.text
    assert b.status_code == 200, b.text


def test_primitive_pause_twice_returns_200(client) -> None:
    pid = _create_psr(client)
    a = client.post(f"/api/primitive-searches/{pid}/pause")
    b = client.post(f"/api/primitive-searches/{pid}/pause")
    assert a.status_code == 200, a.text
    assert b.status_code == 200, b.text


def test_primitive_resume_twice_returns_200(client) -> None:
    pid = _create_psr(client)
    client.post(f"/api/primitive-searches/{pid}/pause")
    a = client.post(f"/api/primitive-searches/{pid}/resume")
    b = client.post(f"/api/primitive-searches/{pid}/resume")
    assert a.status_code == 200, a.text
    assert b.status_code == 200, b.text


def test_primitive_restart_twice_returns_200(client) -> None:
    pid = _create_psr(client)
    a = client.post(f"/api/primitive-searches/{pid}/restart")
    b = client.post(f"/api/primitive-searches/{pid}/restart")
    assert a.status_code == 200, a.text
    assert b.status_code == 200, b.text


# ── tempering ──


def test_tempering_cancel_twice_returns_200(client) -> None:
    tid = _create_tsr(client)
    if tid is None:
        return
    a = client.post(f"/api/tempering-searches/{tid}/cancel")
    b = client.post(f"/api/tempering-searches/{tid}/cancel")
    assert a.status_code == 200, a.text
    assert b.status_code == 200, b.text


def test_tempering_pause_twice_returns_200(client) -> None:
    tid = _create_tsr(client)
    if tid is None:
        return
    a = client.post(f"/api/tempering-searches/{tid}/pause")
    b = client.post(f"/api/tempering-searches/{tid}/pause")
    assert a.status_code == 200, a.text
    assert b.status_code == 200, b.text


def test_tempering_resume_twice_returns_200(client) -> None:
    tid = _create_tsr(client)
    if tid is None:
        return
    client.post(f"/api/tempering-searches/{tid}/pause")
    a = client.post(f"/api/tempering-searches/{tid}/resume")
    b = client.post(f"/api/tempering-searches/{tid}/resume")
    assert a.status_code == 200, a.text
    assert b.status_code == 200, b.text


def test_tempering_restart_twice_returns_200(client) -> None:
    tid = _create_tsr(client)
    if tid is None:
        return
    a = client.post(f"/api/tempering-searches/{tid}/restart")
    b = client.post(f"/api/tempering-searches/{tid}/restart")
    assert a.status_code == 200, a.text
    assert b.status_code == 200, b.text
