# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 5.4 — schema v1 → v2 in-place migration.

The v2 schema (see schema.sql) introduces three typed result tables
that mirror the C++ analysis structs exactly:

  equidistribution_result  ⇄  MatricialEquidResult { ecart, se, verified }
  collision_free_result    ⇄  CollisionFreeResult  { ecart_cf, secf, verified }
  tuplets_result           ⇄  TupletsRunResult     { tupd, tuph, gap, DELTA,
                                                     pourcentage, firstpart_*,
                                                     secondpart_* }

This module backfills those tables from the legacy `test_result` rows
that v1 wrote. Each `test_result` row is dispatched by `test_type`:

  - "equidistribution" / "matricial" → equidistribution_result
  - "collision_free" / "collisions" → collision_free_result
  - "tuplets" → tuplets_result

The legacy `test_result` table is left in place so currently-deployed
routes keep working. Phase 5.4b/c will flip routes/workers off it and
drop it in Phase 8.

The migration is **idempotent**: re-running it skips any tested_gen_id
that already has a typed row of the matching kind.
"""

from __future__ import annotations

import json
import sqlite3
from typing import Any

_EQUID_TYPES = frozenset({"equidistribution", "matricial", "matricial_equidistribution"})
_CF_TYPES = frozenset({"collision_free", "collisions", "cf"})
_TUPLETS_TYPES = frozenset({"tuplets", "tuplet"})


def migrate_v2_inplace(conn: sqlite3.Connection) -> dict[str, int]:
    """Backfill v2 typed result tables from legacy `test_result` rows.

    Returns a counters dict with keys:
      - equidistribution_inserted
      - collision_free_inserted
      - tuplets_inserted
      - skipped_unknown_type
      - skipped_already_present
      - skipped_bad_detail

    Caller is responsible for committing the connection and bumping
    `schema_version`. A single transaction wraps every insert so a
    failure leaves the DB at the previous version.
    """
    counters = {
        "equidistribution_inserted": 0,
        "collision_free_inserted": 0,
        "tuplets_inserted": 0,
        "skipped_unknown_type": 0,
        "skipped_already_present": 0,
        "skipped_bad_detail": 0,
    }

    rows = conn.execute(
        "SELECT id, tested_gen_id, test_type, test_config, "
        "       se, secf, score, detail "
        "FROM test_result ORDER BY id"
    ).fetchall()

    for row in rows:
        tr_id = row[0]
        tested_gen_id = row[1]
        test_type = (row[2] or "").lower()
        test_config = row[3] or "{}"
        se = row[4]
        secf = row[5]
        # `score` (row[6]) was duplicative of `se`; not preserved separately.
        try:
            detail = json.loads(row[7]) if row[7] else {}
        except (json.JSONDecodeError, TypeError):
            counters["skipped_bad_detail"] += 1
            continue

        if test_type in _EQUID_TYPES:
            if _has_typed_row(conn, "equidistribution_result", tested_gen_id, test_config):
                counters["skipped_already_present"] += 1
                continue
            _insert_equid(conn, tested_gen_id, test_config, detail, se)
            counters["equidistribution_inserted"] += 1
        elif test_type in _CF_TYPES:
            if _has_typed_row(conn, "collision_free_result", tested_gen_id, test_config):
                counters["skipped_already_present"] += 1
                continue
            _insert_cf(conn, tested_gen_id, test_config, detail, secf if secf is not None else 0)
            counters["collision_free_inserted"] += 1
        elif test_type in _TUPLETS_TYPES:
            if _has_typed_row(conn, "tuplets_result", tested_gen_id, test_config):
                counters["skipped_already_present"] += 1
                continue
            _insert_tuplets(conn, tested_gen_id, test_config, detail)
            counters["tuplets_inserted"] += 1
        else:
            counters["skipped_unknown_type"] += 1

    return counters


def _has_typed_row(
    conn: sqlite3.Connection, table: str, tested_gen_id: int, test_config: str
) -> bool:
    cur = conn.execute(
        f"SELECT 1 FROM {table} WHERE tested_gen_id = ? AND test_config = ? LIMIT 1",
        (tested_gen_id, test_config),
    )
    return cur.fetchone() is not None


def _insert_equid(
    conn: sqlite3.Connection,
    tested_gen_id: int,
    test_config_json: str,
    detail: dict[str, Any],
    se: int | None,
) -> None:
    """Reconstruct an equidistribution_result row from a legacy detail blob.

    Legacy `detail` for equidistribution carries:
      {"se": int, "ecart": {"l": gap, ...}}  (sparse, sentinel-filtered)

    The new schema stores the full ecart array. Any missing positions
    are filled with 0 (plan documents this convention: 0 = verified,
    INT_MAX = unverified beyond maxl). For migrated rows we cannot
    distinguish those cases, so we use 0 throughout — this is lossy
    but safe for downstream routes that mostly read `se`.
    """
    cfg = json.loads(test_config_json) if test_config_json else {}
    Lmax = int(cfg.get("Lmax", cfg.get("L", _max_l_from_ecart(detail))))
    L = int(cfg.get("L", Lmax))
    kg = int(cfg.get("kg", _infer_kg(conn, tested_gen_id) or 0))
    ecart_arr = [0] * (Lmax + 1)
    for k, v in (detail.get("ecart") or {}).items():
        try:
            idx = int(k)
        except (TypeError, ValueError):
            continue
        if 0 <= idx <= Lmax:
            ecart_arr[idx] = int(v)
    se_val = int(se) if se is not None else int(detail.get("se", 0))
    conn.execute(
        "INSERT INTO equidistribution_result"
        "(tested_gen_id, test_config, kg, L, Lmax, ecart_json, se, verified) "
        "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
        (tested_gen_id, test_config_json, kg, L, Lmax,
         json.dumps(ecart_arr), se_val, 1),
    )


def _insert_cf(
    conn: sqlite3.Connection,
    tested_gen_id: int,
    test_config_json: str,
    detail: dict[str, Any],
    secf: int,
) -> None:
    cfg = json.loads(test_config_json) if test_config_json else {}
    L = int(cfg.get("L", 0))
    kg = int(cfg.get("kg", _infer_kg(conn, tested_gen_id) or 0))
    raw = detail.get("ecart_cf") or detail.get("ecart") or {}
    ecart_arr = [0] * (kg + 1)
    if isinstance(raw, dict):
        for k, v in raw.items():
            try:
                idx = int(k)
            except (TypeError, ValueError):
                continue
            if 0 <= idx <= kg:
                ecart_arr[idx] = int(v)
    elif isinstance(raw, list):
        for i, v in enumerate(raw):
            if i <= kg:
                ecart_arr[i] = int(v)
    conn.execute(
        "INSERT INTO collision_free_result"
        "(tested_gen_id, test_config, kg, L, ecart_cf_json, secf, verified) "
        "VALUES (?, ?, ?, ?, ?, ?, ?)",
        (tested_gen_id, test_config_json, kg, L,
         json.dumps(ecart_arr), int(secf), 1),
    )


def _insert_tuplets(
    conn: sqlite3.Connection,
    tested_gen_id: int,
    test_config_json: str,
    detail: dict[str, Any],
) -> None:
    """The v1 tuplets path never wrote a typed `detail` blob, so we
    insert a placeholder row that preserves the test_config for later
    re-computation. firstpart_*/secondpart_* default to 0.0 sentinels.
    """
    cfg = json.loads(test_config_json) if test_config_json else {}
    L = int(cfg.get("L", 0))
    kg = int(cfg.get("kg", _infer_kg(conn, tested_gen_id) or 0))
    tuph_in = cfg.get("tuph") or detail.get("tuph") or [0]
    tupd = int(cfg.get("tupd", len(tuph_in) - 1 if tuph_in else 0))
    testtype = int(cfg.get("testtype", 1))
    threshold = float(cfg.get("threshold", 0.0))
    gap_list = detail.get("gap") or [0.0] * ((tuph_in[1] if len(tuph_in) > 1 else 0) + 1)
    DELTA = detail.get("DELTA") or [0.0] * (tupd + 1)
    pourcentage = detail.get("pourcentage") or [0.0] * (tupd + 1)
    conn.execute(
        "INSERT INTO tuplets_result"
        "(tested_gen_id, test_config, kg, L, tupd, testtype, threshold, "
        " tuph_json, gap_json, DELTA_json, pourcentage_json, "
        " firstpart_max, firstpart_sum, secondpart_max, secondpart_sum) "
        "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        (
            tested_gen_id, test_config_json, kg, L, tupd, testtype, threshold,
            json.dumps(list(tuph_in)),
            json.dumps(list(gap_list)),
            json.dumps(list(DELTA)),
            json.dumps(list(pourcentage)),
            float(detail.get("firstpart_max", 0.0)),
            float(detail.get("firstpart_sum", 0.0)),
            float(detail.get("secondpart_max", 0.0)),
            float(detail.get("secondpart_sum", 0.0)),
        ),
    )


def _max_l_from_ecart(detail: dict[str, Any]) -> int:
    raw = detail.get("ecart") or {}
    if not isinstance(raw, dict) or not raw:
        return 0
    try:
        return max(int(k) for k in raw.keys())
    except ValueError:
        return 0


def _infer_kg(conn: sqlite3.Connection, tested_gen_id: int) -> int | None:
    row = conn.execute(
        "SELECT k_g FROM tested_generator WHERE id = ?", (tested_gen_id,)
    ).fetchone()
    return int(row[0]) if row else None
