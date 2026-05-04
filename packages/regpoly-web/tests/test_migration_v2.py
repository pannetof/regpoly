"""Phase 5.4 — schema v1 → v2 migration tests.

Builds a hand-crafted v1 DB (the Phase-0 schema, without the typed
result tables and with schema_version=1), then drives the migration
two ways:

  1. Direct: call migrate_v2_inplace(conn) and check counters + rows.
  2. Via init_sync: opening the older DB through the normal startup
     path should trigger the in-place upgrade automatically.

Also covers idempotency (running twice is a no-op) and the fresh-DB
case (init_sync on a brand new file lands at v2 with no work to do).
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import pytest


# ── v1 schema (subset reproducing the structure we need to migrate) ──

V1_SCHEMA = """
CREATE TABLE schema_version (
    version    INTEGER NOT NULL,
    applied_at TEXT    NOT NULL DEFAULT (datetime('now'))
);

CREATE TABLE primitive_search_run (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    family            TEXT NOT NULL,
    L                 INTEGER NOT NULL,
    k                 INTEGER NOT NULL,
    structural_params TEXT NOT NULL,
    fixed_params      TEXT NOT NULL,
    max_tries         INTEGER,
    max_seconds       REAL,
    status            TEXT NOT NULL DEFAULT 'pending',
    tries_done        INTEGER NOT NULL DEFAULT 0,
    found_count       INTEGER NOT NULL DEFAULT 0,
    elapsed_seconds   REAL,
    error_message     TEXT,
    created_at        TEXT NOT NULL DEFAULT (datetime('now')),
    started_at        TEXT,
    finished_at       TEXT
);

CREATE TABLE primitive_generator (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    search_run_id     INTEGER,
    family            TEXT NOT NULL,
    L                 INTEGER NOT NULL,
    k                 INTEGER NOT NULL,
    structural_params TEXT NOT NULL,
    search_params     TEXT NOT NULL,
    all_params        TEXT NOT NULL,
    found_at_try      INTEGER,
    created_at        TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE TABLE tested_generator (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    search_run_id     INTEGER,
    Lmax              INTEGER NOT NULL,
    k_g               INTEGER NOT NULL,
    J                 INTEGER NOT NULL,
    created_at        TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE TABLE tested_generator_component (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    tested_gen_id     INTEGER NOT NULL,
    component_index   INTEGER NOT NULL,
    generator_id      INTEGER,
    family            TEXT NOT NULL,
    L                 INTEGER NOT NULL,
    k                 INTEGER NOT NULL,
    all_params        TEXT NOT NULL,
    tempering_params  TEXT NOT NULL,
    UNIQUE(tested_gen_id, component_index)
);

CREATE TABLE tempering_search_component (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    search_run_id     INTEGER NOT NULL,
    component_index   INTEGER NOT NULL,
    shared_with_component INTEGER,
    tempering_config  TEXT NOT NULL,
    UNIQUE(search_run_id, component_index)
);

CREATE TABLE tempering_search_generator (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    component_id      INTEGER NOT NULL,
    generator_id      INTEGER NOT NULL
);

CREATE TABLE search_progress (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    search_type       TEXT NOT NULL,
    search_run_id     INTEGER NOT NULL,
    tries_done        INTEGER NOT NULL DEFAULT 0,
    found_count       INTEGER NOT NULL DEFAULT 0,
    current_info      TEXT,
    message           TEXT,
    updated_at        TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE TABLE yaml_import (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    file_path         TEXT NOT NULL UNIQUE,
    import_type       TEXT NOT NULL,
    row_count         INTEGER NOT NULL,
    imported_at       TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE TABLE test_result (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    tested_gen_id     INTEGER NOT NULL,
    test_type         TEXT NOT NULL,
    test_config       TEXT NOT NULL,
    se                INTEGER,
    is_me             INTEGER,
    secf              INTEGER,
    is_cf             INTEGER,
    score             REAL,
    detail            TEXT NOT NULL,
    elapsed_seconds   REAL,
    created_at        TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE TABLE tempering_search_run (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    test_type         TEXT NOT NULL,
    test_config       TEXT NOT NULL,
    Lmax              INTEGER NOT NULL,
    nb_tries          INTEGER NOT NULL,
    optimizer_config  TEXT,
    status            TEXT NOT NULL DEFAULT 'pending',
    combos_total      INTEGER,
    combos_done       INTEGER NOT NULL DEFAULT 0,
    best_se           INTEGER,
    elapsed_seconds   REAL,
    error_message     TEXT,
    created_at        TEXT NOT NULL DEFAULT (datetime('now')),
    started_at        TEXT,
    finished_at       TEXT
);
"""


def _build_v1_db(path: Path) -> None:
    conn = sqlite3.connect(str(path))
    try:
        conn.executescript(V1_SCHEMA)
        conn.execute("INSERT INTO schema_version(version) VALUES (1)")
        conn.execute(
            "INSERT INTO tested_generator(id, search_run_id, Lmax, k_g, J) "
            "VALUES (?, ?, ?, ?, ?)",
            (10, None, 32, 19937, 1),
        )
        conn.execute(
            "INSERT INTO tested_generator(id, search_run_id, Lmax, k_g, J) "
            "VALUES (?, ?, ?, ?, ?)",
            (20, None, 32, 521, 1),
        )
        conn.execute(
            "INSERT INTO tested_generator(id, search_run_id, Lmax, k_g, J) "
            "VALUES (?, ?, ?, ?, ?)",
            (30, None, 32, 1024, 1),
        )
        # Equidistribution row with a sparse `ecart` map.
        conn.execute(
            "INSERT INTO test_result"
            "(tested_gen_id, test_type, test_config, se, is_me, score, detail) "
            "VALUES (?, ?, ?, ?, ?, ?, ?)",
            (10, "equidistribution",
             json.dumps({"L": 32, "Lmax": 32, "kg": 19937}),
             7, 0, 7.0,
             json.dumps({"se": 7, "ecart": {"5": 1, "12": 2, "30": 4}})),
        )
        # Collision-free row.
        conn.execute(
            "INSERT INTO test_result"
            "(tested_gen_id, test_type, test_config, se, is_me, secf, is_cf, score, detail) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (20, "collision_free",
             json.dumps({"L": 32, "kg": 521}),
             None, None, 3, 1, 3.0,
             json.dumps({"ecart_cf": {"100": 1, "300": 2}})),
        )
        # Tuplets row.
        conn.execute(
            "INSERT INTO test_result"
            "(tested_gen_id, test_type, test_config, se, is_me, score, detail) "
            "VALUES (?, ?, ?, ?, ?, ?, ?)",
            (30, "tuplets",
             json.dumps({"L": 32, "kg": 1024, "tupd": 4,
                         "tuph": [0, 8, 8, 8, 8],
                         "testtype": 1, "threshold": 1.0}),
             None, None, None,
             json.dumps({})),
        )
        conn.commit()
    finally:
        conn.close()


def test_migrate_v2_inplace_backfills_typed_tables(tmp_path: Path) -> None:
    db_path = tmp_path / "v1.db"
    _build_v1_db(db_path)

    # The web app's schema.sql is what defines the v2 typed tables.
    from regpoly_web.config import SCHEMA_PATH
    from regpoly_web.migrations.v2 import migrate_v2_inplace

    conn = sqlite3.connect(str(db_path))
    try:
        conn.executescript(Path(SCHEMA_PATH).read_text())
        counters = migrate_v2_inplace(conn)
        conn.commit()

        assert counters["equidistribution_inserted"] == 1
        assert counters["collision_free_inserted"] == 1
        assert counters["tuplets_inserted"] == 1
        assert counters["skipped_unknown_type"] == 0
        assert counters["skipped_already_present"] == 0

        eq = conn.execute(
            "SELECT tested_gen_id, kg, L, Lmax, ecart_json, se, verified "
            "FROM equidistribution_result"
        ).fetchone()
        assert eq[0] == 10
        assert eq[1] == 19937
        assert eq[2] == 32
        assert eq[3] == 32
        ecart = json.loads(eq[4])
        assert len(ecart) == 33  # Lmax+1
        assert ecart[5] == 1
        assert ecart[12] == 2
        assert ecart[30] == 4
        assert eq[5] == 7
        assert eq[6] == 1

        cf = conn.execute(
            "SELECT tested_gen_id, kg, L, ecart_cf_json, secf, verified "
            "FROM collision_free_result"
        ).fetchone()
        assert cf[0] == 20
        assert cf[1] == 521
        ecart_cf = json.loads(cf[3])
        assert len(ecart_cf) == 522  # kg+1
        assert ecart_cf[100] == 1
        assert ecart_cf[300] == 2
        assert cf[4] == 3

        tp = conn.execute(
            "SELECT tested_gen_id, kg, L, tupd, testtype, threshold, "
            " tuph_json, gap_json, DELTA_json, pourcentage_json "
            "FROM tuplets_result"
        ).fetchone()
        assert tp[0] == 30
        assert tp[1] == 1024
        assert tp[3] == 4
        assert tp[4] == 1
        assert tp[5] == pytest.approx(1.0)
        assert json.loads(tp[6]) == [0, 8, 8, 8, 8]
    finally:
        conn.close()


def test_migrate_v2_is_idempotent(tmp_path: Path) -> None:
    db_path = tmp_path / "v1.db"
    _build_v1_db(db_path)

    from regpoly_web.config import SCHEMA_PATH
    from regpoly_web.migrations.v2 import migrate_v2_inplace

    conn = sqlite3.connect(str(db_path))
    try:
        conn.executescript(Path(SCHEMA_PATH).read_text())
        first = migrate_v2_inplace(conn)
        conn.commit()

        second = migrate_v2_inplace(conn)
        conn.commit()

        assert first["equidistribution_inserted"] == 1
        assert second["equidistribution_inserted"] == 0
        assert second["skipped_already_present"] == 3

        # Each typed table should have exactly one row.
        for tbl in ("equidistribution_result", "collision_free_result", "tuplets_result"):
            n = conn.execute(f"SELECT COUNT(*) FROM {tbl}").fetchone()[0]
            assert n == 1
    finally:
        conn.close()


def test_init_sync_upgrades_v1_db_to_v2(tmp_path: Path) -> None:
    """Opening an older DB through the normal startup path should
    trigger the in-place upgrade automatically."""
    db_path = tmp_path / "v1.db"
    _build_v1_db(db_path)

    from regpoly_web.database import SCHEMA_VERSION, init_sync

    init_sync(str(db_path))

    conn = sqlite3.connect(str(db_path))
    try:
        version = conn.execute(
            "SELECT MAX(version) FROM schema_version"
        ).fetchone()[0]
        assert version == SCHEMA_VERSION

        # Backfilled rows should be present.
        assert conn.execute(
            "SELECT COUNT(*) FROM equidistribution_result"
        ).fetchone()[0] == 1
        assert conn.execute(
            "SELECT COUNT(*) FROM collision_free_result"
        ).fetchone()[0] == 1
        assert conn.execute(
            "SELECT COUNT(*) FROM tuplets_result"
        ).fetchone()[0] == 1
    finally:
        conn.close()


def test_init_sync_fresh_db_lands_at_v2(tmp_path: Path) -> None:
    db_path = tmp_path / "fresh.db"

    from regpoly_web.database import SCHEMA_VERSION, init_sync

    init_sync(str(db_path))

    conn = sqlite3.connect(str(db_path))
    try:
        version = conn.execute(
            "SELECT MAX(version) FROM schema_version"
        ).fetchone()[0]
        assert version == SCHEMA_VERSION

        # Typed tables exist but are empty.
        for tbl in ("equidistribution_result", "collision_free_result", "tuplets_result"):
            n = conn.execute(f"SELECT COUNT(*) FROM {tbl}").fetchone()[0]
            assert n == 0
    finally:
        conn.close()
