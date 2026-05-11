# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""pytest fixtures for regpoly-web (PostgreSQL backend).

The dockerize-and-deploy migration switched the backing store from
SQLite to PostgreSQL 16. This conftest spins up an in-process PG
instance (via the ``pgserver`` package, which vendors PG binaries —
no system PG required) and hands out per-test databases.

Two fixture chains preserved from the SQLite era:

- **Per-test ephemeral**: ``tmp_db_url`` / ``web_settings`` / ``client``
  — each test gets a fresh database; rows do not leak across tests.
- **Session-shared seeded**: ``seeded_db_url`` / ``seeded_db`` /
  ``seeded_client`` — a single seeded database used by the v2 redesign
  contract tests (test_dashboard_api.py, test_history_endpoint.py,
  test_filter_endpoints.py, …).

Backward-compat fixture names: the old ``tmp_db_path`` /
``seeded_db_path`` are still emitted (returning DSNs, not file paths)
for tests that haven't been ported yet.
"""

from __future__ import annotations

import json
import threading
import uuid
from collections.abc import Iterator
from pathlib import Path

import psycopg
import pytest
from fastapi.testclient import TestClient

# ─── Session-wide pgserver instance ─────────────────────────────────

_pg_lock = threading.Lock()
_pg_server = None


@pytest.fixture(scope="session")
def _pg_uri(tmp_path_factory: pytest.TempPathFactory) -> str:
    """Start a single in-process PG server for the whole pytest session.

    pgserver vendors the PG binaries; no system installation is needed.
    The server lives in a session tmp dir; it is stopped at session
    teardown via the trailing yield.
    """
    global _pg_server
    import pgserver

    with _pg_lock:
        if _pg_server is None:
            data_dir = tmp_path_factory.mktemp("pgserver_data")
            _pg_server = pgserver.get_server(str(data_dir))

    return _pg_server.get_uri()


def _admin_dsn(uri: str) -> str:
    """Strip ``/<dbname>`` so we can connect to the default ``postgres``
    database for CREATE DATABASE / DROP DATABASE."""
    # pgserver hands out URIs of the form
    # `postgresql://postgres:@/postgres?host=/tmp/...`. We need the
    # trailing `?host=...` preserved; just rewrite the database name.
    if "?" in uri:
        prefix, rest = uri.split("?", 1)
        params = "?" + rest
    else:
        prefix, params = uri, ""
    if "/" in prefix.replace("://", "", 1):
        prefix = prefix.rsplit("/", 1)[0]
    return f"{prefix}/postgres{params}"


def _make_url(admin_uri: str, dbname: str) -> str:
    if "?" in admin_uri:
        prefix, rest = admin_uri.split("?", 1)
        params = "?" + rest
    else:
        prefix, params = admin_uri, ""
    base = prefix.rsplit("/", 1)[0]
    return f"{base}/{dbname}{params}"


def _create_db(admin_dsn: str, dbname: str) -> None:
    """CREATE DATABASE outside any transaction (PG quirk)."""
    conn = psycopg.connect(admin_dsn, autocommit=True)
    try:
        conn.execute(f'CREATE DATABASE "{dbname}"')
    finally:
        conn.close()


def _drop_db(admin_dsn: str, dbname: str) -> None:
    conn = psycopg.connect(admin_dsn, autocommit=True)
    try:
        conn.execute(
            f'SELECT pg_terminate_backend(pid) FROM pg_stat_activity '
            f"WHERE datname = '{dbname}'"
        )
        conn.execute(f'DROP DATABASE IF EXISTS "{dbname}"')
    finally:
        conn.close()


# ─── Per-test ephemeral DB ───────────────────────────────────────────

@pytest.fixture
def tmp_db_url(_pg_uri: str) -> Iterator[str]:
    """A fresh PG database for one test. Created + migrated, dropped
    at teardown."""
    from regpoly_web.database import init_db_sync

    admin_dsn = _admin_dsn(_pg_uri)
    dbname = f"test_{uuid.uuid4().hex[:12]}"
    _create_db(admin_dsn, dbname)
    url = _make_url(admin_dsn, dbname)
    init_db_sync(url)
    try:
        yield url
    finally:
        _drop_db(admin_dsn, dbname)


# Backward-compat alias for tests that still ask for `tmp_db_path`.
@pytest.fixture
def tmp_db_path(tmp_db_url: str) -> str:
    """Legacy name; returns a PG DSN, NOT a filesystem path."""
    return tmp_db_url


@pytest.fixture
def web_settings(tmp_db_url: str):
    """Settings for an app rooted at the per-test DB."""
    from regpoly_web.config import Settings

    return Settings(
        db_url=tmp_db_url,
        reload=False,
        pool_size=1,
        # Confine import-dir to a tmp tree so the path-traversal
        # security test can prove the guard is active.
        import_root=Path("/tmp").resolve(),
    )


@pytest.fixture
def client(web_settings) -> Iterator[TestClient]:
    """A FastAPI TestClient with the app's lifespan executed
    (DB pool opened, catalog loaded, scheduler running)."""
    from regpoly_web.app import create_app

    app = create_app(web_settings)
    with TestClient(app) as c:
        yield c


# ─── Session-shared seeded DB ────────────────────────────────────────

@pytest.fixture(scope="session")
def seeded_db_url(_pg_uri: str) -> Iterator[str]:
    """One PG database shared across the seeded session. Rows seeded
    by the ``seeded_db`` fixture below.
    """
    from regpoly_web.database import init_db_sync

    admin_dsn = _admin_dsn(_pg_uri)
    dbname = "test_seeded"
    _drop_db(admin_dsn, dbname)
    _create_db(admin_dsn, dbname)
    url = _make_url(admin_dsn, dbname)
    init_db_sync(url)
    try:
        yield url
    finally:
        _drop_db(admin_dsn, dbname)


# Backward-compat alias.
@pytest.fixture(scope="session")
def seeded_db_path(seeded_db_url: str) -> str:
    """Legacy name; returns a PG DSN, NOT a filesystem path."""
    return seeded_db_url


@pytest.fixture(scope="session")
def seeded_db(seeded_db_url: str) -> str:
    """Initialize seed rows for the v2 redesign contract tests.

    Seeded rows:
      - 1 completed primitive_search_run for MTGen
      - 1 primitive_generator row attached to that run
      - 1 tested_generator id=4242 with one component
      - 1 equidistribution_result (Lmax=32, sparse δ at ℓ=5,12; se=3)
    """
    conn = psycopg.connect(seeded_db_url, autocommit=False)
    try:
        with conn.cursor() as cur:
            cur.execute(
                'INSERT INTO primitive_search_run'
                '(family, l, k, structural_params, fixed_params, '
                ' status, tries_done, found_count, elapsed_seconds) '
                "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s) RETURNING id",
                ("MTGen", 19937, 32, json.dumps({}), json.dumps({}),
                 "completed", 100, 1, 0.5),
            )
            psr_id = cur.fetchone()[0]

            # A 32-entry gap profile seeded so the /generators/{id}
            # detail page's equidist_chart has data to render in tests.
            pis_gaps = [0] * 32
            pis_gaps[5] = 1
            pis_gaps[12] = 2
            cur.execute(
                'INSERT INTO primitive_generator'
                '(search_run_id, family, l, k, structural_params, '
                ' search_params, all_params, found_at_try, char_poly, '
                ' hamming_weight, pis_se, pis_gaps, pis_computed_at) '
                "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, NOW())",
                (psr_id, "MTGen", 19937, 32,
                 json.dumps({"L": 19937}),
                 json.dumps({"a": "0x9908b0df"}),
                 json.dumps({"L": 19937, "a": "0x9908b0df"}),
                 1, "0xdeadbeef", 135, 0, json.dumps(pis_gaps)),
            )

            cur.execute(
                'INSERT INTO tested_generator'
                '(id, search_run_id, lmax, k_g, j) '
                "VALUES (%s, %s, %s, %s, %s)",
                (4242, None, 32, 19937, 1),
            )
            cur.execute(
                'INSERT INTO tested_generator_component'
                '(tested_gen_id, component_index, family, l, k, '
                ' all_params, tempering_params) '
                "VALUES (%s, %s, %s, %s, %s, %s, %s)",
                (4242, 0, "MTGen", 19937, 32,
                 json.dumps({"L": 19937}), json.dumps([])),
            )
            ecart = [0] * 33
            ecart[5] = 1
            ecart[12] = 2
            cur.execute(
                'INSERT INTO equidistribution_result'
                '(tested_gen_id, test_config, kg, l, lmax, ecart_json, se, '
                ' verified, elapsed_seconds) '
                "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)",
                (4242, json.dumps({"type": "equidistribution", "L": 32}),
                 19937, 32, 32, json.dumps(ecart), 3, True, 0.1),
            )
        conn.commit()
    finally:
        conn.close()
    return seeded_db_url


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

    settings = Settings(db_url=seeded_db, reload=False, pool_size=1)
    app = create_app(settings)
    with TestClient(app) as c:
        yield c
