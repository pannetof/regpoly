# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Tests for the WELL mat_types renumbering migration."""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import pytest

from regpoly_web.config import SCHEMA_PATH
from regpoly_web.migrations.well_mat_types import (
    OLD_TO_NEW,
    MARKER_PATH,
    already_migrated,
    migrate,
)


@pytest.fixture
def fresh_db(tmp_path: Path) -> sqlite3.Connection:
    db_path = tmp_path / "regpoly.db"
    conn = sqlite3.connect(str(db_path))
    conn.execute("PRAGMA foreign_keys=ON")
    conn.executescript(Path(SCHEMA_PATH).read_text())
    return conn


def _ins_well_gen(
    conn: sqlite3.Connection, mat_types: list[int], extras: dict | None = None
) -> int:
    extras = extras or {}
    payload = {"mat_types": mat_types, **extras}
    cur = conn.execute(
        """INSERT INTO primitive_generator
           (family, L, k, structural_params, search_params, all_params)
           VALUES ('WELLGen', 32, 19937, ?, ?, ?)""",
        (json.dumps(payload), json.dumps(payload), json.dumps(payload)),
    )
    return cur.lastrowid  # type: ignore[return-value]


def test_remaps_well_row_in_place(fresh_db: sqlite3.Connection) -> None:
    # Legacy WELL19937a row 1: types {0,0,3,0,1,0,0,0}.
    rid = _ins_well_gen(fresh_db, [0, 0, 3, 0, 1, 0, 0, 0])

    counters = migrate(fresh_db)
    fresh_db.commit()

    assert counters.rows_touched >= 1
    row = fresh_db.execute(
        "SELECT all_params FROM primitive_generator WHERE id=?", (rid,)
    ).fetchone()
    payload = json.loads(row[0])
    # Expected post-remap: {3,3,2,3,1,3,3,3} (paper Table II WELL19937a).
    assert payload["mat_types"] == [3, 3, 2, 3, 1, 3, 3, 3]


def test_skips_obsolete_old_type6(fresh_db: sqlite3.Connection) -> None:
    rid = _ins_well_gen(fresh_db, [0, 6, 0, 0, 0, 0, 0, 0])

    counters = migrate(fresh_db)

    assert counters.rows_touched == 0
    assert counters.rows_skipped == 1
    assert ("primitive_generator", rid) in counters.skipped_old_type6
    # The row's payload must be unchanged.
    fresh_db.commit()
    row = fresh_db.execute(
        "SELECT all_params FROM primitive_generator WHERE id=?", (rid,)
    ).fetchone()
    assert json.loads(row[0])["mat_types"] == [0, 6, 0, 0, 0, 0, 0, 0]


def test_dry_run_does_not_persist(fresh_db: sqlite3.Connection) -> None:
    rid = _ins_well_gen(fresh_db, [0, 0, 3, 0, 1, 0, 0, 0])
    fresh_db.commit()  # persist the seed row before the dry-run rollback.

    migrate(fresh_db, dry_run=True)
    fresh_db.rollback()

    row = fresh_db.execute(
        "SELECT all_params FROM primitive_generator WHERE id=?", (rid,)
    ).fetchone()
    assert json.loads(row[0])["mat_types"] == [0, 0, 3, 0, 1, 0, 0, 0]
    assert not already_migrated(fresh_db)


def test_idempotency_marker_blocks_rerun(fresh_db: sqlite3.Connection) -> None:
    _ins_well_gen(fresh_db, [0, 0, 3, 0, 1, 0, 0, 0])
    migrate(fresh_db)
    fresh_db.commit()

    with pytest.raises(RuntimeError, match=MARKER_PATH):
        migrate(fresh_db)


def test_force_overrides_marker(fresh_db: sqlite3.Connection) -> None:
    _ins_well_gen(fresh_db, [0, 0, 3, 0, 1, 0, 0, 0])
    migrate(fresh_db)
    fresh_db.commit()

    # A second pass with --force should not crash. The values are
    # already in 0..6 range so OLD_TO_NEW would shuffle them again,
    # which is the expected (and footgun-y) behaviour of --force.
    counters = migrate(fresh_db, force=True)
    assert counters.rows_scanned >= 1


def test_non_well_rows_are_not_touched(fresh_db: sqlite3.Connection) -> None:
    payload = {"a": 0xDEADBEEF, "w": 32}
    fresh_db.execute(
        """INSERT INTO primitive_generator
           (family, L, k, structural_params, search_params, all_params)
           VALUES ('MTGen', 32, 19937, ?, ?, ?)""",
        (json.dumps(payload), json.dumps(payload), json.dumps(payload)),
    )

    counters = migrate(fresh_db)
    assert counters.rows_touched == 0


def test_old_to_new_table_matches_paper() -> None:
    # Sanity-check the OLD_TO_NEW table mirrors the cpp source.
    assert OLD_TO_NEW == [3, 1, 4, 2, 5, 6, -1, 0]
