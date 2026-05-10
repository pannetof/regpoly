# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Tests for the unauthenticated /healthz probe.

Required by the dockerize plan Phase 2: Caddy exempts /healthz from
basic auth; the compose `web` service uses ``curl /healthz`` as its
healthcheck; UptimeRobot's external monitor keys on the ``"db":"up"``
substring. The endpoint must return 200 with a stable shape on
success and 503 with a degraded payload when PG is unreachable.
"""

from __future__ import annotations


def test_healthz_returns_ok_when_pg_up(client) -> None:
    r = client.get("/healthz")
    assert r.status_code == 200, r.text
    body = r.json()
    assert body == {"status": "ok", "db": "up"}


def test_healthz_unauthenticated(client) -> None:
    """The probe must not require basic auth.

    The TestClient doesn't enforce Caddy basic auth (that runs in
    front of the app), but we can still assert the route is mounted
    OUTSIDE the /api prefix — Caddy exempts the exact path /healthz,
    so a route at /api/healthz would be auth-gated by mistake.
    """
    spec = client.get("/openapi.json").json()
    assert "/healthz" in spec["paths"], (
        "/healthz must be at the root, not under /api/"
    )
    assert "/api/healthz" not in spec["paths"]


def test_healthz_returns_503_when_pool_unreachable(client) -> None:
    """If the pool's connection acquire times out / errors, return 503.

    Replace the pool with one whose ``connection()`` raises an
    OperationalError so the route's exception path runs.
    """
    import psycopg

    class _BrokenPool:
        def connection(self):
            raise psycopg.OperationalError("simulated outage")

    saved = client.app.state.dbpool
    client.app.state.dbpool = _BrokenPool()
    try:
        r = client.get("/healthz")
        assert r.status_code == 503, r.text
        body = r.json()
        assert body["status"] == "degraded"
        assert body["db"] == "down"
    finally:
        client.app.state.dbpool = saved


def test_healthz_skipped_when_no_pool(client) -> None:
    """When ``state.dbpool`` is None the route returns ok with
    ``db: "skipped"`` — this is the shape the route uses if a future
    test app is built without a pool. Confirms the null-check path."""
    saved = client.app.state.dbpool
    client.app.state.dbpool = None
    try:
        r = client.get("/healthz")
        assert r.status_code == 200, r.text
        assert r.json() == {"status": "ok", "db": "skipped"}
    finally:
        client.app.state.dbpool = saved
