# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""WELL matrices payload migration: flat triple → structured map.

The paper-aligned redesign replaced `mat_types`/`mat_pi`/`mat_pu`
(three flat lists per WELL row in the DB JSON) with a single
`matrices: {T0..T7}` map. This migration walks every JSON payload
that may contain a WELL row, converts the legacy flat triple to the
new structured form, and removes the legacy keys.

For WELL ints that pre-date the M0..M6 rename (the "well-mat-types-
migrated" marker is absent), it ALSO applies the integer remap from
the legacy 0..7 numbering to the paper 0..6 numbering as part of the
shape conversion. So one CLI invocation handles every state.

State machine driven by `yaml_import` marker rows:

  - Neither marker present: full conversion (integer remap + shape
    conversion + drop legacy keys). Used for fresh DBs that never ran
    the old migration.
  - Only `well-mat-types-migrated` present (45cf403's marker): shape
    conversion only — integer remap was done; just restructure.
    The old marker is deleted; the new marker is written.
  - Only `well-matrices-restructured` present (this commit's marker):
    skip; nothing to do.
  - Both present: skip; nothing to do; print warning.

Tables and JSON columns walked:

  primitive_search_run.{structural_params, fixed_params}
  primitive_generator.{structural_params, search_params, all_params}
  tested_generator_component.{all_params, tempering_params}
  tempering_search_component.tempering_config
  test_result.detail

Filter: rows where `family = 'WELLGen'` (the canonical name written
by `routes/import_export.py:_import_one_generators_file` after
`resolve_family()` normalisation). Untagged JSON columns are walked
for any payload that contains `mat_types`/`mat_pi`/`mat_pu` keys.

Idempotency: re-runs are gated by the new marker. `--force` overrides.

The migration unconditionally **deletes** the legacy keys after
conversion — single source of truth.
"""

from __future__ import annotations

import json
import sqlite3
from dataclasses import dataclass, field
from typing import Any

# Legacy → paper Mi numbering (from the M0..M6 rename, commit 3e83478).
OLD_TO_NEW: list[int] = [3, 1, 4, 2, 5, 6, -1, 0]

OLD_MARKER = "well-mat-types-migrated"            # 45cf403
NEW_MARKER = "well-matrices-restructured"         # this commit


@dataclass
class MigrationCounters:
    rows_scanned: int = 0
    rows_touched: int = 0
    rows_skipped: int = 0
    skipped_old_type6: list[tuple[str, int]] = field(default_factory=list)
    skipped_other: list[tuple[str, int, str]] = field(default_factory=list)
    per_table: dict[str, dict[str, int]] = field(default_factory=dict)
    state: str = "fresh"  # 'fresh' | 'transitional' | 'done' | 'inconsistent'


# ─── d_s and test-mask decoders (mirror C++ well_legacy_decode) ──────────

M32 = 0xFFFFFFFF


def _decode_d_s_mask(mask: int, ctx: str) -> int:
    inv = (~int(mask)) & M32
    if inv == 0 or inv & (inv - 1):
        raise ValueError(
            f"{ctx}: cannot decode WELL paper M6 d_s mask 0x{int(mask):08x}: "
            f"expected exactly one zero bit (got {bin(inv).count('1')})"
        )
    lsb = (inv & -inv).bit_length() - 1
    return 31 - lsb


def _decode_test_mask(mask: int, ctx: str) -> int:
    m = int(mask) & M32
    if m == 0 or m & (m - 1):
        raise ValueError(
            f"{ctx}: cannot decode WELL paper M6 test mask 0x{int(mask):08x}: "
            f"expected exactly one set bit"
        )
    lsb = (m & -m).bit_length() - 1
    return 31 - lsb


# ─── Payload conversion ──────────────────────────────────────────────────


def _flat_to_structured(
    payload: dict[str, Any],
    *,
    apply_integer_remap: bool,
    ctx: str,
) -> tuple[dict[str, Any], bool, str | None]:
    """Convert a payload with `mat_types`/`mat_pi`/`mat_pu` to
    `matrices: {T0..T7}` form. Returns ``(new_payload, touched, error)``.
    Already-structured payloads are passed through untouched.
    """
    if not isinstance(payload, dict):
        return payload, False, None
    mat_types = payload.get("mat_types")
    if mat_types is None:
        return payload, False, None
    mat_pi = payload.get("mat_pi") or [0] * 24
    mat_pu = payload.get("mat_pu") or [0] * 24
    if not isinstance(mat_types, list) or len(mat_types) < 8:
        return payload, False, "mat_types is not an 8-int list"

    matrices: dict[str, dict] = {}
    for j in range(8):
        raw = int(mat_types[j])
        if apply_integer_remap:
            if raw < 0 or raw >= 8 or OLD_TO_NEW[raw] < 0:
                return payload, False, (
                    f"obsolete legacy WELL transformation type {raw}"
                )
            Mi = OLD_TO_NEW[raw]
        else:
            if raw < 0 or raw > 6:
                return payload, False, (
                    f"out-of-range Mi {raw} (post-rename payload should be 0..6)"
                )
            Mi = raw
        pi0 = int(mat_pi[j * 3]) if len(mat_pi) > j * 3 else 0
        pu0 = int(mat_pu[j * 3]) if len(mat_pu) > j * 3 else 0
        pu1 = int(mat_pu[j * 3 + 1]) if len(mat_pu) > j * 3 + 1 else 0
        pu2 = int(mat_pu[j * 3 + 2]) if len(mat_pu) > j * 3 + 2 else 0
        slot_ctx = f"{ctx} slot T{j}"
        entry: dict[str, int] = {"M": Mi}
        if Mi in (0, 1):
            pass
        elif Mi in (2, 3):
            entry["t"] = pi0
        elif Mi == 4:
            entry["a"] = pu0 & M32
        elif Mi == 5:
            entry["t"] = pi0
            entry["b"] = pu0 & M32
        elif Mi == 6:
            try:
                entry["q"] = pi0
                entry["t"] = _decode_test_mask(pu2, slot_ctx)
                entry["s"] = _decode_d_s_mask(pu1, slot_ctx)
                entry["a"] = pu0 & M32
            except ValueError as exc:
                return payload, False, str(exc)
        matrices[f"T{j}"] = entry

    new_payload = {k: v for k, v in payload.items()
                   if k not in ("mat_types", "mat_pi", "mat_pu")}
    new_payload["matrices"] = matrices
    return new_payload, True, None


# ─── Marker state machine ────────────────────────────────────────────────


def _get_state(conn: sqlite3.Connection) -> str:
    has_old = conn.execute(
        "SELECT 1 FROM yaml_import WHERE file_path = ?", (OLD_MARKER,)
    ).fetchone() is not None
    has_new = conn.execute(
        "SELECT 1 FROM yaml_import WHERE file_path = ?", (NEW_MARKER,)
    ).fetchone() is not None
    if has_new and not has_old: return "done"
    if has_old and has_new:     return "done"  # treat both-present as done
    if has_old:                 return "transitional"
    return "fresh"


def already_migrated(conn: sqlite3.Connection) -> bool:
    return _get_state(conn) == "done"


# ─── Walking ─────────────────────────────────────────────────────────────


def _bump(counters: MigrationCounters, table: str, key: str) -> None:
    counters.per_table.setdefault(table, {})
    counters.per_table[table][key] = counters.per_table[table].get(key, 0) + 1


def _walk_rows(
    conn: sqlite3.Connection,
    table: str,
    json_columns: list[str],
    where: str,
    counters: MigrationCounters,
    *,
    apply_integer_remap: bool,
    dry_run: bool,
) -> None:
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
            new_payload, touched, err = _flat_to_structured(
                payload,
                apply_integer_remap=apply_integer_remap,
                ctx=f"{table}.{col}#{rid}",
            )
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


def migrate(
    conn: sqlite3.Connection,
    *,
    dry_run: bool = False,
    force: bool = False,
) -> MigrationCounters:
    """Apply the WELL matrices migration. Returns counters and final
    marker state. Caller must commit (or rollback) the transaction.
    """
    counters = MigrationCounters()
    state = _get_state(conn)
    counters.state = state

    if state == "done" and not force:
        raise RuntimeError(
            f"DB already migrated (yaml_import marker '{NEW_MARKER}' present)."
            " Use --force to re-run."
        )

    apply_integer_remap = (state == "fresh")

    well_filter = "family = 'WELLGen'"
    _walk_rows(conn, "primitive_search_run",
               ["structural_params", "fixed_params"],
               well_filter, counters,
               apply_integer_remap=apply_integer_remap, dry_run=dry_run)
    _walk_rows(conn, "primitive_generator",
               ["structural_params", "search_params", "all_params"],
               well_filter, counters,
               apply_integer_remap=apply_integer_remap, dry_run=dry_run)
    _walk_rows(conn, "tested_generator_component",
               ["all_params", "tempering_params"],
               well_filter, counters,
               apply_integer_remap=apply_integer_remap, dry_run=dry_run)
    _walk_rows(conn, "tempering_search_component",
               ["tempering_config"],
               "1=1", counters,
               apply_integer_remap=apply_integer_remap, dry_run=dry_run)
    _walk_rows(conn, "test_result",
               ["detail"],
               "1=1", counters,
               apply_integer_remap=apply_integer_remap, dry_run=dry_run)

    if not dry_run:
        # Replace markers: delete old (if present), insert new.
        conn.execute(
            "DELETE FROM yaml_import WHERE file_path = ?", (OLD_MARKER,)
        )
        conn.execute(
            "INSERT OR REPLACE INTO yaml_import"
            " (file_path, import_type, row_count) VALUES (?, ?, ?)",
            (NEW_MARKER, NEW_MARKER, counters.rows_touched),
        )

    return counters
