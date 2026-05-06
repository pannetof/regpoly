# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""WELL mat_types renumbering migration: legacy 0..7 → paper M0..M6.

Walks every JSON column that may contain a WELL ``mat_types`` array
and remaps each element via OLD_TO_NEW. The on-disk schema is
untouched — only payload bytes change. Old type 6 (multi-shift XOR)
has no paper equivalent and is rejected per row.

Tables and columns walked:

  primitive_search_run.{structural_params, fixed_params}
  primitive_generator.{structural_params, search_params, all_params}
  tempering_search_component.tempering_config
  tested_generator_component.all_params
  test_result.detail

Filter: rows where ``family = 'WELLGen'`` (the canonical name written
by ``routes/import_export.py:_import_one_generators_file`` after
``resolve_family()`` normalisation). For ``tempering_search_component``
and ``test_result``, ``family`` is not on the row directly; we walk
the JSON conservatively and skip any payload that doesn't contain a
``mat_types`` key.

Idempotency: re-runs are gated by a marker row in ``yaml_import``
(``import_type = 'well-mat-types-migrated'``). ``--force`` overrides
the gate.
"""

from __future__ import annotations

import json
import sqlite3
from dataclasses import dataclass, field
from typing import Any

OLD_TO_NEW: list[int] = [3, 1, 4, 2, 5, 6, -1, 0]
MARKER_PATH = "well-mat-types-migrated"


@dataclass
class MigrationCounters:
    rows_scanned: int = 0
    rows_touched: int = 0
    rows_skipped: int = 0
    skipped_old_type6: list[tuple[str, int]] = field(default_factory=list)
    skipped_other: list[tuple[str, int, str]] = field(default_factory=list)
    per_table: dict[str, dict[str, int]] = field(default_factory=dict)


def _remap_mat_types(payload: Any) -> tuple[Any, bool, str | None]:
    """Apply OLD_TO_NEW to ``payload['mat_types']`` if present.

    Returns ``(new_payload, touched, error)`` where:
      - ``touched`` is True iff ``mat_types`` was found and remapped,
      - ``error`` is a string when an old type 6 (or out-of-range int)
        was found and the row should be skipped.
    """
    if not isinstance(payload, dict):
        return payload, False, None
    mt = payload.get("mat_types")
    if not isinstance(mt, list):
        return payload, False, None
    if not all(isinstance(t, int) for t in mt):
        return payload, False, None

    # Heuristic: if every value is already in 0..6 AND the *count* of
    # legal-but-rare new values (0 = M0 = ZERO) plus any value that's
    # impossible under legacy (none) suggests the payload is already
    # post-rename, leave it alone. We can't *prove* it from numbers
    # alone — but the caller gates the entire migration via the
    # MARKER_PATH idempotency check, so this function is only called
    # for never-migrated DBs.
    new_mt: list[int] = []
    for t in mt:
        if t < 0 or t >= len(OLD_TO_NEW):
            return payload, False, f"out-of-range type {t}"
        nt = OLD_TO_NEW[t]
        if nt < 0:
            return payload, False, "obsolete old type 6 (multi-shift)"
        new_mt.append(nt)
    if new_mt == mt:
        return payload, False, None
    new_payload = dict(payload)
    new_payload["mat_types"] = new_mt
    return new_payload, True, None


def _bump(counters: MigrationCounters, table: str, key: str) -> None:
    counters.per_table.setdefault(table, {})
    counters.per_table[table][key] = counters.per_table[table].get(key, 0) + 1


def _walk_well_rows(
    conn: sqlite3.Connection,
    table: str,
    json_columns: list[str],
    where: str,
    counters: MigrationCounters,
    *,
    dry_run: bool,
) -> None:
    """Walk ``SELECT id, <json_columns> FROM <table> WHERE <where>``.

    For each row, decode each JSON column, apply ``_remap_mat_types``,
    and ``UPDATE`` the row in place if any column was touched. Rows
    with old-type-6 payloads are logged and skipped wholesale.
    """
    cols = ", ".join(json_columns)
    rows = conn.execute(
        f"SELECT id, {cols} FROM {table} WHERE {where}"
    ).fetchall()
    for row in rows:
        counters.rows_scanned += 1
        _bump(counters, table, "scanned")
        rid = row[0]
        new_values: list[Any] = []
        any_error: str | None = None
        any_touched = False
        for j, col in enumerate(json_columns):
            raw = row[j + 1]
            if raw is None or raw == "":
                new_values.append(raw)
                continue
            try:
                payload = json.loads(raw)
            except (json.JSONDecodeError, TypeError):
                new_values.append(raw)
                continue
            new_payload, touched, err = _remap_mat_types(payload)
            if err is not None:
                any_error = err
                break
            if touched:
                any_touched = True
                new_values.append(json.dumps(new_payload, sort_keys=True))
            else:
                new_values.append(raw)
        if any_error is not None:
            counters.rows_skipped += 1
            _bump(counters, table, "skipped")
            if "obsolete" in any_error:
                counters.skipped_old_type6.append((table, rid))
            else:
                counters.skipped_other.append((table, rid, any_error))
            continue
        if not any_touched:
            continue
        counters.rows_touched += 1
        _bump(counters, table, "touched")
        if not dry_run:
            placeholders = ", ".join(f"{c} = ?" for c in json_columns)
            conn.execute(
                f"UPDATE {table} SET {placeholders} WHERE id = ?",
                (*new_values, rid),
            )


def already_migrated(conn: sqlite3.Connection) -> bool:
    """Return True iff the marker row exists in yaml_import."""
    cur = conn.execute(
        "SELECT 1 FROM yaml_import WHERE file_path = ?", (MARKER_PATH,)
    )
    return cur.fetchone() is not None


def write_marker(conn: sqlite3.Connection, rows_touched: int) -> None:
    conn.execute(
        "INSERT OR REPLACE INTO yaml_import"
        " (file_path, import_type, row_count) VALUES (?, ?, ?)",
        (MARKER_PATH, "well-mat-types-migrated", rows_touched),
    )


def migrate(
    conn: sqlite3.Connection,
    *,
    dry_run: bool = False,
    force: bool = False,
) -> MigrationCounters:
    """Apply the WELL mat_types renumbering to all relevant rows.

    Returns a populated ``MigrationCounters``. The caller is
    responsible for committing or rolling back. We do NOT open our
    own transaction — Python's sqlite3 already runs DML inside an
    implicit transaction, and stacking a ``BEGIN`` raises
    "cannot start a transaction within a transaction". For
    concurrent-writer protection use ``conn.execute('BEGIN
    IMMEDIATE')`` *before* calling ``migrate()`` if needed.
    """
    counters = MigrationCounters()

    if not force and already_migrated(conn):
        raise RuntimeError(
            "DB already migrated (yaml_import marker '"
            + MARKER_PATH
            + "' present). Use --force to re-run."
        )

    # Family-tagged tables: filter to WELLGen rows only.
    well_filter = "family = 'WELLGen'"
    _walk_well_rows(
        conn,
        "primitive_search_run",
        ["structural_params", "fixed_params"],
        well_filter,
        counters,
        dry_run=dry_run,
    )
    _walk_well_rows(
        conn,
        "primitive_generator",
        ["structural_params", "search_params", "all_params"],
        well_filter,
        counters,
        dry_run=dry_run,
    )
    _walk_well_rows(
        conn,
        "tested_generator_component",
        ["all_params"],
        well_filter,
        counters,
        dry_run=dry_run,
    )

    # Untagged JSON columns: scan everything, remap only if mat_types
    # is present. ``_remap_mat_types`` is a no-op on payloads without
    # the key. (test_result.detail and tempering_search_component
    # rarely contain mat_types, but we walk for completeness.)
    _walk_well_rows(
        conn,
        "tempering_search_component",
        ["tempering_config"],
        "1=1",
        counters,
        dry_run=dry_run,
    )
    _walk_well_rows(
        conn,
        "test_result",
        ["detail"],
        "1=1",
        counters,
        dry_run=dry_run,
    )

    if not dry_run:
        write_marker(conn, counters.rows_touched)

    return counters
