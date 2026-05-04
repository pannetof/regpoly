"""Persistence helpers for typed analysis results (Phase 5.4b).

`save_typed_result` takes a Python analysis-result object
(EquidistributionResults / CollisionFreeResults / TupletsResults from
regpoly.analyses), inspects its class, and writes the matching row to
the v2 typed table introduced in Phase 5.4a.

The legacy `test_result` insert remains in the calling code (dual
write) until Phase 5.4c flips the route reads onto the typed tables.
That keeps the test_result.detail JSON consumers (tempered route, UI)
working unchanged through the transition.
"""

from __future__ import annotations

import json
import sqlite3
import time
from typing import Any

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


def now_seconds() -> float:
    """Tiny wrapper so callers can record `elapsed = now() - t0` without
    importing time themselves."""
    return time.time()
