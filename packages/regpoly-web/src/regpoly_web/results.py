# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Persistence helpers for typed analysis results (Phase 5.4b/c).

`save_typed_result` takes a Python analysis-result object
(EquidistributionResults / CollisionFreeResults / TupletsResults from
regpoly.analyses), inspects its class, and writes the matching row to
the v2 typed table introduced in Phase 5.4a.

`read_typed_results_*` produce the legacy test_result-shaped dict
list from those typed tables, so the existing UI consumers (which
expect `{test_type, se, is_me, secf, is_cf, score, detail, ...}`)
keep working without template changes.
"""

from __future__ import annotations

import json
import sqlite3
import time
from typing import Any

import aiosqlite
from regpoly.analyses.collision_free_results import CollisionFreeResults
from regpoly.analyses.equidistribution_results import EquidistributionResults
from regpoly.analyses.tuplets_results import TupletsResults


def save_typed_result(
    conn: sqlite3.Connection,
    tested_gen_id: int,
    test_config: dict[str, Any],
    result: Any,
    *,
    kg: int,
    L: int,
    elapsed_seconds: float | None = None,
) -> int | None:
    """Insert one row into the typed table that matches `result`'s class.

    Returns the new row id, or None if the result class is not recognised
    (caller should still write to legacy test_result in that case).

    `kg` is the combined-generator k_g; `L` is the analysis-side L (for
    equidistribution this is typically the active Lmax). The caller
    has these from the Combination it just exercised.
    """
    test_config_json = json.dumps(test_config, sort_keys=True)

    if isinstance(result, EquidistributionResults):
        return _insert_equid(conn, tested_gen_id, test_config_json,
                             result, kg=kg, L=L, elapsed=elapsed_seconds)
    if isinstance(result, CollisionFreeResults):
        return _insert_cf(conn, tested_gen_id, test_config_json,
                          result, kg=kg, L=L, elapsed=elapsed_seconds)
    if isinstance(result, TupletsResults):
        return _insert_tuplets(conn, tested_gen_id, test_config_json,
                               result, kg=kg, L=L, elapsed=elapsed_seconds)
    return None


def _insert_equid(
    conn: sqlite3.Connection,
    tested_gen_id: int,
    test_config_json: str,
    result: EquidistributionResults,
    *,
    kg: int,
    L: int,
    elapsed: float | None,
) -> int:
    Lmax = len(result.ecart) - 1  # ecart is 1-indexed, length L+1
    ecart_arr = [int(x) for x in result.ecart]
    cur = conn.execute(
        "INSERT INTO equidistribution_result"
        "(tested_gen_id, test_config, kg, L, Lmax, ecart_json, se, "
        " verified, elapsed_seconds) "
        "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        (tested_gen_id, test_config_json, int(kg), int(L), int(Lmax),
         json.dumps(ecart_arr), int(result.se),
         1 if getattr(result, "_verified", False) else 0,
         elapsed),
    )
    return int(cur.lastrowid)


def _insert_cf(
    conn: sqlite3.Connection,
    tested_gen_id: int,
    test_config_json: str,
    result: CollisionFreeResults,
    *,
    kg: int,
    L: int,
    elapsed: float | None,
) -> int:
    ecart_cf_arr = [int(x) for x in result.ecart_cf]
    cur = conn.execute(
        "INSERT INTO collision_free_result"
        "(tested_gen_id, test_config, kg, L, ecart_cf_json, secf, "
        " verified, elapsed_seconds) "
        "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
        (tested_gen_id, test_config_json, int(kg), int(L),
         json.dumps(ecart_cf_arr), int(result.secf),
         1 if getattr(result, "_verified", False) else 0,
         elapsed),
    )
    return int(cur.lastrowid)


def _insert_tuplets(
    conn: sqlite3.Connection,
    tested_gen_id: int,
    test_config_json: str,
    result: TupletsResults,
    *,
    kg: int,
    L: int,
    elapsed: float | None,
) -> int:
    cur = conn.execute(
        "INSERT INTO tuplets_result"
        "(tested_gen_id, test_config, kg, L, tupd, testtype, threshold, "
        " tuph_json, gap_json, DELTA_json, pourcentage_json, "
        " firstpart_max, firstpart_sum, secondpart_max, secondpart_sum, "
        " elapsed_seconds) "
        "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        (
            tested_gen_id, test_config_json, int(kg), int(L),
            int(result.tupd), int(result.testtype), float(result.treshold),
            json.dumps([int(x) for x in result.tuph]),
            json.dumps([float(x) for x in result.gap]),
            json.dumps([float(x) for x in result.DELTA]),
            json.dumps([float(x) for x in result.pourcentage]),
            float(result.firstpart_max), float(result.firstpart_sum),
            float(result.secondpart_max), float(result.secondpart_sum),
            elapsed,
        ),
    )
    return int(cur.lastrowid)


def _safe_json(v):
    """Tolerant decoder: PG JSONB columns return ``dict``/``list``
    natively via psycopg3, but legacy SQLite-text payloads round-trip
    as strings. Accept both shapes."""
    if v is None or v == "":
        return None
    if isinstance(v, (dict, list)):
        return v
    return json.loads(v)


def now_seconds() -> float:
    """Tiny wrapper so callers can record `elapsed = now() - t0` without
    importing time themselves."""
    return time.time()


# ── Read path ─────────────────────────────────────────────────────────
#
# The route layer used to call `SELECT * FROM test_result`, building
# UI-shaped dicts directly from those columns. After 5.4c the typed
# tables are the source of truth, so we union the three SELECTs into
# the same legacy shape the templates already consume:
#
#   {test_type, test_config, se, is_me, secf, is_cf, score, detail,
#    elapsed_seconds, created_at}
#
# The synthetic fields (is_me / is_cf / score / detail.ecart) are
# derived from the typed columns:
#
#   - is_me  = (table == equid)   AND verified == 1 AND se == 0
#   - is_cf  = (table == cf)      AND verified == 1 AND secf == 0
#   - score  = se (equid) / secf (cf) / firstpart_max+secondpart_max (tuplets)
#   - detail = {ecart: {l → gap}}  (equid; sentinel-filtered as before)
#              {ecart_cf: {l → gap}}  (cf)
#              {gap, DELTA, pourcentage, ...}  (tuplets — verbatim)


_EQUID_SELECT = (
    "SELECT 'equidistribution' AS kind, id, tested_gen_id, test_config, "
    "       kg, L, Lmax, ecart_json, NULL::jsonb AS ecart_cf_json, "
    "       se, NULL::integer AS secf, verified, "
    "       NULL::integer AS tupd, NULL::integer AS testtype, "
    "       NULL::double precision AS threshold, "
    "       NULL::jsonb AS tuph_json, NULL::jsonb AS gap_json, "
    "       NULL::jsonb AS delta_json, "
    "       NULL::jsonb AS pourcentage_json, "
    "       NULL::double precision AS firstpart_max, "
    "       NULL::double precision AS firstpart_sum, "
    "       NULL::double precision AS secondpart_max, "
    "       NULL::double precision AS secondpart_sum, "
    "       elapsed_seconds, created_at "
    "FROM equidistribution_result WHERE tested_gen_id = %s"
)
_CF_SELECT = (
    "SELECT 'collision_free' AS kind, id, tested_gen_id, test_config, "
    "       kg, L, NULL::integer AS Lmax, NULL::jsonb AS ecart_json, "
    "       ecart_cf_json, "
    "       NULL::integer AS se, secf, verified, "
    "       NULL::integer AS tupd, NULL::integer AS testtype, "
    "       NULL::double precision AS threshold, "
    "       NULL::jsonb AS tuph_json, NULL::jsonb AS gap_json, "
    "       NULL::jsonb AS delta_json, "
    "       NULL::jsonb AS pourcentage_json, "
    "       NULL::double precision AS firstpart_max, "
    "       NULL::double precision AS firstpart_sum, "
    "       NULL::double precision AS secondpart_max, "
    "       NULL::double precision AS secondpart_sum, "
    "       elapsed_seconds, created_at "
    "FROM collision_free_result WHERE tested_gen_id = %s"
)
_TUPLETS_SELECT = (
    "SELECT 'tuplets' AS kind, id, tested_gen_id, test_config, "
    "       kg, L, NULL::integer AS Lmax, NULL::jsonb AS ecart_json, "
    "       NULL::jsonb AS ecart_cf_json, "
    "       NULL::integer AS se, NULL::integer AS secf, "
    "       NULL::boolean AS verified, "
    "       tupd, testtype, threshold, "
    "       tuph_json, gap_json, delta_json, pourcentage_json, "
    "       firstpart_max, firstpart_sum, secondpart_max, secondpart_sum, "
    "       elapsed_seconds, created_at "
    "FROM tuplets_result WHERE tested_gen_id = %s"
)

_UNION_SELECT = (
    f"{_EQUID_SELECT} UNION ALL {_CF_SELECT} UNION ALL {_TUPLETS_SELECT} "
    f"ORDER BY created_at, id"
)

# Sentinel: any single ecart entry larger than this is an INT_MAX
# placeholder for "uncomputed beyond maxl" — filter it out so it
# never reaches the UI. Mirrors the rule baked into the legacy
# tasks/tempering.py::_build_result_detail.
_ECART_SENTINEL = 1 << 40


def _row_to_legacy_dict(row: Any) -> dict[str, Any]:
    """Convert a unioned row (dict-like) into the legacy shape."""
    kind = row["kind"]
    test_config = row["test_config"]
    try:
        test_config = _safe_json(test_config) if test_config else None
    except (TypeError, json.JSONDecodeError):
        pass
    base: dict[str, Any] = {
        "test_type": kind,
        "test_config": test_config,
        "se": row["se"],
        "is_me": None,
        "secf": row["secf"],
        "is_cf": None,
        "score": None,
        "detail": {},
        "elapsed_seconds": row["elapsed_seconds"],
        "created_at": row["created_at"],
    }
    if kind == "equidistribution":
        ecart_arr = _safe_json(row["ecart_json"]) if row["ecart_json"] else []
        sparse = {
            str(i): int(v)
            for i, v in enumerate(ecart_arr)
            if v and 0 < int(v) < _ECART_SENTINEL
        }
        base["detail"] = {"ecart": sparse} if sparse else {}
        base["is_me"] = (
            bool(row["verified"]) and row["se"] is not None and int(row["se"]) == 0
        )
        if row["se"] is not None:
            base["score"] = float(row["se"])
    elif kind == "collision_free":
        ecart_cf_arr = _safe_json(row["ecart_cf_json"]) if row["ecart_cf_json"] else []
        sparse = {
            str(i): int(v)
            for i, v in enumerate(ecart_cf_arr)
            if v and 0 < int(v) < _ECART_SENTINEL
        }
        base["detail"] = {"ecart_cf": sparse} if sparse else {}
        base["is_cf"] = (
            bool(row["verified"]) and row["secf"] is not None and int(row["secf"]) == 0
        )
        if row["secf"] is not None:
            base["score"] = float(row["secf"])
    elif kind == "tuplets":
        gap = _safe_json(row["gap_json"]) if row["gap_json"] else []
        DELTA = _safe_json(row["DELTA_json"]) if row["DELTA_json"] else []
        pourcentage = _safe_json(row["pourcentage_json"]) if row["pourcentage_json"] else []
        base["detail"] = {
            "tupd": row["tupd"],
            "tuph": _safe_json(row["tuph_json"]) if row["tuph_json"] else [],
            "gap": gap,
            "DELTA": DELTA,
            "pourcentage": pourcentage,
            "firstpart_max": row["firstpart_max"],
            "firstpart_sum": row["firstpart_sum"],
            "secondpart_max": row["secondpart_max"],
            "secondpart_sum": row["secondpart_sum"],
        }
        # Score for ranking: max-style sum of the two firstpart/secondpart
        # max columns; matches the Python TupletsResults.is_ok() ordering.
        if row["firstpart_max"] is not None and row["secondpart_max"] is not None:
            base["score"] = float(row["firstpart_max"]) + float(row["secondpart_max"])
    return base


async def read_typed_results_async(
    db: aiosqlite.Connection, tested_gen_id: int
) -> list[dict[str, Any]]:
    """Async read for the FastAPI route layer."""
    out: list[dict[str, Any]] = []
    async with db.execute(
        _UNION_SELECT, (tested_gen_id, tested_gen_id, tested_gen_id)
    ) as cur:
        rows = await cur.fetchall()
    for row in rows:
        out.append(_row_to_legacy_dict(row))
    return out


def read_typed_results_sync(
    conn: sqlite3.Connection, tested_gen_id: int
) -> list[dict[str, Any]]:
    """Sync read for worker / CLI use."""
    rows = conn.execute(
        _UNION_SELECT, (tested_gen_id, tested_gen_id, tested_gen_id)
    ).fetchall()
    return [_row_to_legacy_dict(row) for row in rows]


def find_existing_typed_result(
    conn: sqlite3.Connection,
    tested_gen_id: int,
    test_type: str,
    test_config_json: str,
) -> dict[str, Any] | None:
    """Look up a previously-saved typed result for the same
    (tested_gen_id, test_type, test_config). Returns the legacy-shaped
    dict + a synthetic id, or None.

    Used by tasks/library_test.py to short-circuit when the user
    re-runs a test that's already been computed.
    """
    table = {
        "equidistribution": ("equidistribution_result", _EQUID_SELECT),
        "matricial": ("equidistribution_result", _EQUID_SELECT),
        "collision_free": ("collision_free_result", _CF_SELECT),
        "tuplets": ("tuplets_result", _TUPLETS_SELECT),
    }.get(test_type.lower())
    if table is None:
        return None
    _, select_sql = table
    sql = f"{select_sql} AND test_config = %s ORDER BY id DESC LIMIT 1"
    row = conn.execute(sql, (tested_gen_id, test_config_json)).fetchone()
    if row is None:
        return None
    base = _row_to_legacy_dict(row)
    base["id"] = row["id"]
    return base
