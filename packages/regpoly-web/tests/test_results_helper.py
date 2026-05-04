"""Phase 5.4b — typed-result write helper tests.

Exercises `regpoly_web.results.save_typed_result` against the v2 typed
tables introduced in 5.4a, ensuring each result class lands in the
correct table with the correct typed columns.
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import pytest


def _open_v2_db(tmp_path: Path) -> sqlite3.Connection:
    """Apply the live schema.sql to a fresh DB and insert one
    tested_generator row to satisfy the FK."""
    from regpoly_web.config import SCHEMA_PATH

    db = tmp_path / "results.db"
    conn = sqlite3.connect(str(db))
    conn.executescript(Path(SCHEMA_PATH).read_text())
    conn.execute(
        "INSERT INTO tested_generator(id, search_run_id, Lmax, k_g, J) "
        "VALUES (?, ?, ?, ?, ?)",
        (42, None, 32, 19937, 1),
    )
    conn.commit()
    return conn


def test_save_typed_result_equidistribution(tmp_path: Path) -> None:
    from regpoly.analyses.equidistribution_results import EquidistributionResults

    from regpoly_web.results import save_typed_result

    L = 8
    ecart = [0] * (L + 1)
    ecart[1] = 0
    ecart[3] = 5
    ecart[7] = 12
    res = EquidistributionResults(
        L=L,
        ecart=ecart,
        psi12=[True] * (L + 1),
        se=17,
        verified=True,
        mse=0,
        meverif=True,
        delta=[0] * (L + 1),
    )
    conn = _open_v2_db(tmp_path)
    try:
        rid = save_typed_result(
            conn, 42, {"type": "equidistribution", "L": L},
            res, kg=19937, L=L, elapsed_seconds=0.5,
        )
        conn.commit()
        assert rid is not None

        row = conn.execute(
            "SELECT tested_gen_id, kg, L, Lmax, ecart_json, se, verified, "
            "       elapsed_seconds "
            "FROM equidistribution_result WHERE id = ?", (rid,)
        ).fetchone()
        assert row[0] == 42
        assert row[1] == 19937
        assert row[2] == L
        assert row[3] == L
        assert json.loads(row[4]) == ecart
        assert row[5] == 17
        assert row[6] == 1
        assert row[7] == pytest.approx(0.5)
    finally:
        conn.close()


def test_save_typed_result_collision_free(tmp_path: Path) -> None:
    from regpoly.analyses.collision_free_results import CollisionFreeResults

    from regpoly_web.results import save_typed_result

    kg = 8
    ecart_cf = [0] * (kg + 1)
    ecart_cf[3] = 1
    ecart_cf[7] = 2
    res = CollisionFreeResults(
        ecart_cf=ecart_cf, secf=3, verified=True, msecf=0,
    )
    conn = _open_v2_db(tmp_path)
    try:
        rid = save_typed_result(
            conn, 42, {"type": "collision_free", "L": 32},
            res, kg=kg, L=32,
        )
        conn.commit()
        assert rid is not None

        row = conn.execute(
            "SELECT tested_gen_id, kg, L, ecart_cf_json, secf, verified "
            "FROM collision_free_result WHERE id = ?", (rid,)
        ).fetchone()
        assert row[0] == 42
        assert row[1] == kg
        assert row[2] == 32
        assert json.loads(row[3]) == ecart_cf
        assert row[4] == 3
        assert row[5] == 1
    finally:
        conn.close()


def test_save_typed_result_tuplets(tmp_path: Path) -> None:
    from regpoly.analyses.tuplets_results import TupletsResults

    from regpoly_web.results import save_typed_result

    tupd = 3
    tuph = [0, 4, 4, 4]                                 # 1-indexed
    gap = [0.0, 0.0, 1.0, 2.0, 3.0]                     # length tuph[1]+1 = 5
    DELTA = [0.0, 0.5, 1.0, 1.5]                        # length tupd+1
    pourcentage = [0.0, 0.1, 0.2, 0.3]
    res = TupletsResults(
        tupletsverif=True,
        tupd=tupd,
        tuph=tuph,
        gap=gap,
        DELTA=DELTA,
        pourcentage=pourcentage,
        firstpart_max=2.5,
        firstpart_sum=4.0,
        secondpart_max=1.0,
        secondpart_sum=1.5,
        treshold=3.0,
        testtype=1,
        verified=True,
    )
    conn = _open_v2_db(tmp_path)
    try:
        rid = save_typed_result(
            conn, 42, {"type": "tuplets", "L": 32},
            res, kg=19937, L=32,
        )
        conn.commit()
        assert rid is not None

        row = conn.execute(
            "SELECT tupd, testtype, threshold, "
            "       tuph_json, gap_json, DELTA_json, pourcentage_json, "
            "       firstpart_max, firstpart_sum, secondpart_max, secondpart_sum "
            "FROM tuplets_result WHERE id = ?", (rid,)
        ).fetchone()
        assert row[0] == tupd
        assert row[1] == 1
        assert row[2] == pytest.approx(3.0)
        assert json.loads(row[3]) == tuph
        assert json.loads(row[4]) == gap
        assert json.loads(row[5]) == DELTA
        assert json.loads(row[6]) == pourcentage
        assert row[7] == pytest.approx(2.5)
        assert row[8] == pytest.approx(4.0)
        assert row[9] == pytest.approx(1.0)
        assert row[10] == pytest.approx(1.5)
    finally:
        conn.close()


def test_save_typed_result_unknown_type_returns_none(tmp_path: Path) -> None:
    """Result objects that aren't one of the three recognised classes
    return None so the caller can still write the legacy test_result row."""
    from regpoly_web.results import save_typed_result

    class Unknown:
        se = 0

    conn = _open_v2_db(tmp_path)
    try:
        rid = save_typed_result(
            conn, 42, {"type": "weird"}, Unknown(), kg=1, L=1,
        )
        assert rid is None
    finally:
        conn.close()
