# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Tests for the WELL matrices payload migration (flat triple →
structured `matrices: {T0..T7}` map)."""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import pytest
from regpoly_web.config import SCHEMA_PATH
from regpoly_web.migrations.well_matrices import (
    NEW_MARKER,
    OLD_MARKER,
    OLD_TO_NEW,
    _get_state,
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


# WELL19937a row 1 (paper Table II): legacy mat_types {0,0,3,0,1,0,0,0}
# (legacy 0..7 numbering) decodes to paper-form {3,3,2,3,1,3,3,3} via
# OLD_TO_NEW. mat_pi is [-25,0,0, 27,0,0, 9,0,0, 1,0,0, 0,0,0, -9,0,0,
# -21,0,0, 21,0,0]. mat_pu is all zeros.
LEGACY_WELL19937A = {
    "mat_types": [0, 0, 3, 0, 1, 0, 0, 0],
    "mat_pi": [-25, 0, 0, 27, 0, 0, 9, 0, 0, 1, 0, 0,
                 0, 0, 0, -9, 0, 0, -21, 0, 0, 21, 0, 0],
    "mat_pu": [0] * 24,
}

# Same generator after the 45cf403 integer remap but pre-shape-conversion
# (the "transitional" state).
TRANSITIONAL_WELL19937A = {
    "mat_types": [3, 3, 2, 3, 1, 3, 3, 3],
    "mat_pi": [-25, 0, 0, 27, 0, 0, 9, 0, 0, 1, 0, 0,
                 0, 0, 0, -9, 0, 0, -21, 0, 0, 21, 0, 0],
    "mat_pu": [0] * 24,
}

EXPECTED_STRUCTURED = {
    "matrices": {
        "T0": {"M": 3, "t": -25},
        "T1": {"M": 3, "t":  27},
        "T2": {"M": 2, "t":   9},
        "T3": {"M": 3, "t":   1},
        "T4": {"M": 1},
        "T5": {"M": 3, "t":  -9},
        "T6": {"M": 3, "t": -21},
        "T7": {"M": 3, "t":  21},
    },
}


def _ins_well_gen(conn: sqlite3.Connection, payload: dict) -> int:
    cur = conn.execute(
        """INSERT INTO primitive_generator
           (family, L, k, structural_params, search_params, all_params)
           VALUES ('WELLGen', 32, 19937, ?, ?, ?)""",
        (json.dumps(payload), json.dumps(payload), json.dumps(payload)),
    )
    return cur.lastrowid  # type: ignore[return-value]


def _set_marker(conn: sqlite3.Connection, marker: str) -> None:
    conn.execute(
        "INSERT OR REPLACE INTO yaml_import"
        " (file_path, import_type, row_count) VALUES (?, ?, ?)",
        (marker, marker, 0),
    )


# ─── State machine ───────────────────────────────────────────────────────


def test_state_fresh_db(fresh_db: sqlite3.Connection) -> None:
    assert _get_state(fresh_db) == "fresh"
    assert not already_migrated(fresh_db)


def test_state_with_old_marker_only(fresh_db: sqlite3.Connection) -> None:
    _set_marker(fresh_db, OLD_MARKER)
    fresh_db.commit()
    assert _get_state(fresh_db) == "transitional"
    assert not already_migrated(fresh_db)


def test_state_with_new_marker(fresh_db: sqlite3.Connection) -> None:
    _set_marker(fresh_db, NEW_MARKER)
    fresh_db.commit()
    assert _get_state(fresh_db) == "done"
    assert already_migrated(fresh_db)


def test_state_with_both_markers_treated_as_done(fresh_db: sqlite3.Connection) -> None:
    _set_marker(fresh_db, OLD_MARKER)
    _set_marker(fresh_db, NEW_MARKER)
    fresh_db.commit()
    assert _get_state(fresh_db) == "done"


# ─── Fresh-DB conversion (integer remap + shape) ─────────────────────────


def test_fresh_db_full_conversion(fresh_db: sqlite3.Connection) -> None:
    rid = _ins_well_gen(fresh_db, LEGACY_WELL19937A)

    counters = migrate(fresh_db)
    fresh_db.commit()

    assert counters.state == "fresh"
    assert counters.rows_touched >= 1
    row = fresh_db.execute(
        "SELECT all_params FROM primitive_generator WHERE id=?", (rid,)
    ).fetchone()
    payload = json.loads(row[0])
    # Legacy keys are gone.
    assert "mat_types" not in payload
    assert "mat_pi" not in payload
    assert "mat_pu" not in payload
    # Structured matrices match paper Table II for WELL19937a.
    assert payload["matrices"] == EXPECTED_STRUCTURED["matrices"]
    # New marker is written; old marker is absent.
    assert _get_state(fresh_db) == "done"


# ─── Transitional-DB conversion (shape only) ─────────────────────────────


def test_transitional_db_shape_only(fresh_db: sqlite3.Connection) -> None:
    """A DB that ran 45cf403's CLI has integers already in 0..6 range
    plus the OLD marker. The new CLI should NOT reapply OLD_TO_NEW
    (which would shuffle the integers); only shape-convert."""
    rid = _ins_well_gen(fresh_db, TRANSITIONAL_WELL19937A)
    _set_marker(fresh_db, OLD_MARKER)
    fresh_db.commit()
    assert _get_state(fresh_db) == "transitional"

    counters = migrate(fresh_db)
    fresh_db.commit()

    assert counters.state == "transitional"
    payload = json.loads(fresh_db.execute(
        "SELECT all_params FROM primitive_generator WHERE id=?", (rid,)
    ).fetchone()[0])
    assert payload["matrices"] == EXPECTED_STRUCTURED["matrices"]
    # OLD marker removed, NEW marker set.
    assert _get_state(fresh_db) == "done"
    assert fresh_db.execute(
        "SELECT 1 FROM yaml_import WHERE file_path = ?", (OLD_MARKER,)
    ).fetchone() is None


# ─── Edge cases ──────────────────────────────────────────────────────────


def test_skip_obsolete_old_type6(fresh_db: sqlite3.Connection) -> None:
    bad = {**LEGACY_WELL19937A,
           "mat_types": [0, 6, 0, 0, 0, 0, 0, 0]}    # legacy type 6 → -1
    rid = _ins_well_gen(fresh_db, bad)

    counters = migrate(fresh_db)

    assert counters.rows_skipped == 1
    assert any(t == "primitive_generator" and r == rid
               for t, r in counters.skipped_old_type6)
    # Payload unchanged.
    fresh_db.commit()
    payload = json.loads(fresh_db.execute(
        "SELECT all_params FROM primitive_generator WHERE id=?", (rid,)
    ).fetchone()[0])
    assert payload["mat_types"][1] == 6   # untouched


def test_dry_run_does_not_persist(fresh_db: sqlite3.Connection) -> None:
    rid = _ins_well_gen(fresh_db, LEGACY_WELL19937A)
    fresh_db.commit()

    migrate(fresh_db, dry_run=True)
    fresh_db.rollback()

    payload = json.loads(fresh_db.execute(
        "SELECT all_params FROM primitive_generator WHERE id=?", (rid,)
    ).fetchone()[0])
    assert "mat_types" in payload
    assert "matrices" not in payload
    assert not already_migrated(fresh_db)


def test_idempotency_blocks_rerun(fresh_db: sqlite3.Connection) -> None:
    _ins_well_gen(fresh_db, LEGACY_WELL19937A)
    migrate(fresh_db)
    fresh_db.commit()
    assert already_migrated(fresh_db)

    with pytest.raises(RuntimeError, match=NEW_MARKER):
        migrate(fresh_db)


def test_force_overrides_marker(fresh_db: sqlite3.Connection) -> None:
    _ins_well_gen(fresh_db, LEGACY_WELL19937A)
    migrate(fresh_db)
    fresh_db.commit()
    # Re-run with --force; rows are now structured, so nothing should
    # actually be touched (idempotent on structured form). No exception.
    counters = migrate(fresh_db, force=True)
    assert counters.rows_touched == 0


def test_non_well_rows_untouched(fresh_db: sqlite3.Connection) -> None:
    payload = {"a": 0xDEADBEEF, "w": 32}
    fresh_db.execute(
        """INSERT INTO primitive_generator
           (family, L, k, structural_params, search_params, all_params)
           VALUES ('MTGen', 32, 19937, ?, ?, ?)""",
        (json.dumps(payload), json.dumps(payload), json.dumps(payload)),
    )
    counters = migrate(fresh_db)
    assert counters.rows_touched == 0


def test_old_to_new_table() -> None:
    assert OLD_TO_NEW == [3, 1, 4, 2, 5, 6, -1, 0]
