# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Production-mode shape: web container has no in-process worker pool.

The dockerize plan splits compute off the web container. With
``REGPOLY_WEB_DISABLE_INTERNAL_POOL=1`` the web app must:

- Open the asyncpg pool (so routes can read/write).
- Skip creating the search/analysis :class:`TaskPool` instances.
- Skip starting the in-process schedulers.
- Still accept ``POST /api/primitive-searches`` — the row lands as
  ``status='pending'`` and waits for the worker container to pick
  it up.

A failure mode the plan explicitly calls out: a route that still
called ``app.state.pool.submit(...)`` would raise ``AttributeError``
the first time someone POSTed a search. This test ensures that
contract: routes never touch ``app.state.pool`` after the cutover.
"""

from __future__ import annotations

import json
from pathlib import Path

from fastapi.testclient import TestClient


def test_disable_internal_pool_skips_in_process_pools(tmp_db_url: str) -> None:
    """Lifespan must not allocate ProcessPoolExecutors in prod mode."""
    from regpoly_web.app import create_app
    from regpoly_web.config import Settings

    settings = Settings(
        db_url=tmp_db_url,
        reload=False,
        pool_size=1,
        disable_internal_pool=True,
        import_root=Path("/tmp").resolve(),
    )
    app = create_app(settings)
    with TestClient(app):
        assert app.state.pool is None, (
            "Production-mode app must not create the search TaskPool"
        )
        assert app.state.analysis_pool is None, (
            "Production-mode app must not create the analysis TaskPool"
        )
        # The asyncpg-style pool MUST still be present — routes need it.
        assert app.state.dbpool is not None
        assert app.state.db is not None


def test_post_primitive_search_in_prod_mode_yields_pending(
    tmp_db_url: str,
) -> None:
    """A submit against a worker-less web app lands as `pending`.

    Reproduces the contract the plan codifies in `routes` after the
    Phase 2 worker split: every route is "DB write only", never a
    direct dispatch. Without a worker, the row stays `pending`
    indefinitely; with one, the worker container's scheduler picks
    it up.
    """
    import psycopg
    from regpoly_web.app import create_app
    from regpoly_web.config import Settings

    settings = Settings(
        db_url=tmp_db_url,
        reload=False,
        pool_size=1,
        disable_internal_pool=True,
        import_root=Path("/tmp").resolve(),
    )
    app = create_app(settings)
    with TestClient(app) as client:
        r = client.post(
            "/api/primitive-searches",
            json={
                "family": "MTGen",
                "L": 64,
                "structural_params": {"w": 32, "n": 2, "m": 1, "r": 31, "u": 11},
                "fixed_params": {},
                "max_tries": 1,
            },
        )
        assert r.status_code in (200, 201), r.text
        body = r.json()
        run_id = body["id"]
        assert body["status"] == "pending"

    # After the lifespan exits, the row must STILL be pending — no
    # in-process worker should have touched it. We re-open a sync
    # connection to confirm.
    conn = psycopg.connect(tmp_db_url, autocommit=True)
    try:
        row = conn.execute(
            "SELECT status FROM primitive_search_run WHERE id = %s",
            (run_id,),
        ).fetchone()
    finally:
        conn.close()
    assert row is not None
    # The lifespan teardown reaps pending+running → cancelled in dev
    # mode, but that path only runs when disable_internal_pool=False.
    # In prod mode the row should be untouched.
    assert row[0] == "pending", (
        f"prod-mode app must not mutate pending rows on shutdown; got "
        f"{row[0]!r}"
    )
