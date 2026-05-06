# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""pytest fixtures for regpoly-web.

Phase 5.3 introduces a FastAPI TestClient fixture backed by an
isolated SQLite file in a tmpdir. The v2 redesign (P2 prep) promotes
the e2e-only `seeded_db` fixture to this root conftest so non-e2e
contract tests (test_dashboard_api.py, test_history_endpoint.py,
test_filter_endpoints.py, …) can consume it. The e2e conftest
re-imports it.
"""

from __future__ import annotations

import json
import sqlite3
from collections.abc import Iterator
from pathlib import Path

import pytest
from fastapi.testclient import TestClient


@pytest.fixture
def tmp_db_path(tmp_path: Path) -> str:
    """An ephemeral SQLite file used by a single web app instance.

    SQLite + aiosqlite require a real file path (not :memory:) because
    the worker pool processes open their own connections to the same
    DB. tmp_path gives each test its own isolated DB.
    """
    return str(tmp_path / "test_regpoly.db")


@pytest.fixture
def web_settings(tmp_db_path: str):
    """Settings for an app rooted at the tmp DB.

    Disables hot-reload so the catalog is loaded once at startup
    (matches the production lifespan).
    """
    from regpoly_web.config import Settings

    from pathlib import Path
    return Settings(
        db_path=tmp_db_path,
        reload=False,
        pool_size=1,
        # Confine import-dir to the tmp tree so the path-traversal
        # security test can prove the guard is active.
        import_root=Path(tmp_db_path).parent.resolve(),
    )


@pytest.fixture
def client(web_settings) -> Iterator[TestClient]:
    """A FastAPI TestClient with the app's lifespan executed
    (DB opened, catalog loaded, pool started, scheduler running)."""
    from regpoly_web.app import create_app

    app = create_app(web_settings)
    # `with TestClient(app)` runs the lifespan; without the context
    # manager startup events don't fire and app.state is empty.
    with TestClient(app) as c:
        yield c


# --------------------------------------------------------------------
# Seeded DB shared across the v2 redesign's contract tests.
# Promoted from tests/e2e/conftest.py during P2 prep so non-e2e tests
# can consume it. Session-scoped: rows live for the entire pytest run.
# --------------------------------------------------------------------


@pytest.fixture(scope="session")
def seeded_db_path(tmp_path_factory: pytest.TempPathFactory) -> str:
    """One isolated SQLite file shared across the seeded session."""
    p = tmp_path_factory.mktemp("seeded_db") / "regpoly_seeded.db"
    return str(p)


@pytest.fixture(scope="session")
def seeded_db(seeded_db_path: str) -> str:
    """Initialize the schema and seed enough rows to drive the v2
    redesign's contract tests (and the legacy e2e golden paths).

    Seeded rows:
      - 1 completed primitive_search_run for MTGen (psr_id=auto)
      - 1 primitive_generator row attached to that run
      - 1 tested_generator id=4242 with one component
      - 1 equidistribution_result (Lmax=32, sparse δ at ℓ=5,12; se=3)
    """
    from regpoly_web.database import init_sync

    init_sync(seeded_db_path)

    conn = sqlite3.connect(seeded_db_path)
    conn.row_factory = sqlite3.Row
    try:
        conn.execute(
            "INSERT INTO primitive_search_run"
            "(family, L, k, structural_params, fixed_params, "
            " status, tries_done, found_count, elapsed_seconds) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            ("MTGen", 19937, 32, json.dumps({}), json.dumps({}),
             "completed", 100, 1, 0.5),
        )
        psr_id = conn.execute(
            "SELECT id FROM primitive_search_run ORDER BY id DESC LIMIT 1"
        ).fetchone()[0]

        conn.execute(
            "INSERT INTO primitive_generator"
            "(search_run_id, family, L, k, structural_params, "
            " search_params, all_params, found_at_try, char_poly, "
            " hamming_weight, pis_se, pis_computed_at) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, datetime('now'))",
            (psr_id, "MTGen", 19937, 32,
             json.dumps({"L": 19937}),
             json.dumps({"a": "0x9908b0df"}),
             json.dumps({"L": 19937, "a": "0x9908b0df"}),
             1, "0xdeadbeef", 135, 0),
        )

        conn.execute(
            "INSERT INTO tested_generator"
            "(id, search_run_id, Lmax, k_g, J) VALUES (?, ?, ?, ?, ?)",
            (4242, None, 32, 19937, 1),
        )
        conn.execute(
            "INSERT INTO tested_generator_component"
            "(tested_gen_id, component_index, family, L, k, "
            " all_params, tempering_params) "
            "VALUES (?, ?, ?, ?, ?, ?, ?)",
            (4242, 0, "MTGen", 19937, 32,
             json.dumps({"L": 19937}), json.dumps([])),
        )
        ecart = [0] * 33
        ecart[5] = 1
        ecart[12] = 2
        conn.execute(
            "INSERT INTO equidistribution_result"
            "(tested_gen_id, test_config, kg, L, Lmax, ecart_json, se, "
            " verified, elapsed_seconds) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (4242, json.dumps({"type": "equidistribution", "L": 32}),
             19937, 32, 32, json.dumps(ecart), 3, 1, 0.1),
        )

        conn.commit()
    finally:
        conn.close()
    return seeded_db_path


@pytest.fixture
def seeded_client(seeded_db: str) -> Iterator[TestClient]:
    """A FastAPI TestClient over the seeded DB (session-shared rows).

    Use this for non-e2e contract tests that need the seeded primitive
    search run / generator / tested-generator + result. Each test
    still gets a fresh app instance (lifespan rerun) to avoid state
    leakage across tests.
    """
    from regpoly_web.app import create_app
    from regpoly_web.config import Settings

    settings = Settings(db_path=seeded_db, reload=False, pool_size=1)
    app = create_app(settings)
    with TestClient(app) as c:
        yield c
