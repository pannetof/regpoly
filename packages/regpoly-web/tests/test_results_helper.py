# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 5.4b — typed-result write helper tests.

Exercises `regpoly_web.results.save_typed_result` against the v2 typed
tables introduced in 5.4a, ensuring each result class lands in the
correct table with the correct typed columns.

NOTE (dockerize-plan Phase 1.2): this suite was written against
SQLite directly (sqlite3.connect on the schema.sql file). After the
PG cutover the helpers expect a psycopg connection and the SQL uses
``%s`` placeholders. The suite is skipped pending a PG-aware port
that uses the conftest's pgserver fixture.
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import pytest

pytestmark = pytest.mark.skip(
    reason="SQLite-shaped helper tests; port to pgserver fixture in a follow-up "
    "(dockerize-plan Phase 1.2)"
)


def _open_v2_db(tmp_path: Path) -> sqlite3.Connection:
    """Apply the live schema.sql to a fresh DB and insert one
    tested_generator row to satisfy the FK."""
    from regpoly_web.config import SCHEMA_PATH

    db = tmp_path / "results.db"
    conn = sqlite3.connect(str(db))
    conn.row_factory = sqlite3.Row
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


# ── Phase 5.4c — read path ────────────────────────────────────────────


def _save_one_of_each(conn: sqlite3.Connection, tested_id: int) -> None:
    from regpoly.analyses.collision_free_results import CollisionFreeResults
    from regpoly.analyses.equidistribution_results import EquidistributionResults
    from regpoly.analyses.tuplets_results import TupletsResults
    from regpoly_web.results import save_typed_result

    L = 4
    ecart = [0, 0, 1, 0, 2]
    save_typed_result(
        conn, tested_id, {"type": "equidistribution", "L": L},
        EquidistributionResults(
            L=L, ecart=ecart, psi12=[True] * (L + 1), se=3,
            verified=True, mse=0, meverif=True, delta=[0] * (L + 1),
        ),
        kg=19937, L=L,
    )
    ecart_cf = [0, 0, 0, 1, 0]
    save_typed_result(
        conn, tested_id, {"type": "collision_free", "L": 32},
        CollisionFreeResults(ecart_cf=ecart_cf, secf=1, verified=True, msecf=0),
        kg=521, L=32,
    )
    save_typed_result(
        conn, tested_id, {"type": "tuplets", "L": 32},
        TupletsResults(
            tupletsverif=True, tupd=2, tuph=[0, 4, 4],
            gap=[0.0, 0.0, 1.0, 2.0, 3.0],
            DELTA=[0.0, 0.5, 1.0],
            pourcentage=[0.0, 0.1, 0.2],
            firstpart_max=2.5, firstpart_sum=4.0,
            secondpart_max=1.0, secondpart_sum=1.5,
            treshold=3.0, testtype=1, verified=True,
        ),
        kg=19937, L=32,
    )


def test_read_typed_results_sync_unions_all_three_tables(tmp_path: Path) -> None:
    from regpoly_web.results import read_typed_results_sync

    conn = _open_v2_db(tmp_path)
    try:
        _save_one_of_each(conn, 42)
        conn.commit()

        rows = read_typed_results_sync(conn, 42)
        assert len(rows) == 3
        kinds = {r["test_type"] for r in rows}
        assert kinds == {"equidistribution", "collision_free", "tuplets"}

        eq = next(r for r in rows if r["test_type"] == "equidistribution")
        assert eq["se"] == 3
        assert eq["is_me"] is False  # se != 0
        assert eq["detail"] == {"ecart": {"2": 1, "4": 2}}
        assert eq["score"] == pytest.approx(3.0)

        cf = next(r for r in rows if r["test_type"] == "collision_free")
        assert cf["secf"] == 1
        assert cf["is_cf"] is False  # secf != 0
        assert cf["detail"] == {"ecart_cf": {"3": 1}}

        tp = next(r for r in rows if r["test_type"] == "tuplets")
        assert tp["detail"]["tupd"] == 2
        assert tp["detail"]["firstpart_max"] == pytest.approx(2.5)
    finally:
        conn.close()


def test_find_existing_typed_result_dedupes_by_config(tmp_path: Path) -> None:
    from regpoly_web.results import find_existing_typed_result

    conn = _open_v2_db(tmp_path)
    try:
        _save_one_of_each(conn, 42)
        conn.commit()

        prior = find_existing_typed_result(
            conn, 42, "equidistribution",
            json.dumps({"type": "equidistribution", "L": 4}, sort_keys=True),
        )
        assert prior is not None
        assert prior["se"] == 3
        assert "id" in prior

        # Different config → no hit.
        miss = find_existing_typed_result(
            conn, 42, "equidistribution",
            json.dumps({"type": "equidistribution", "L": 999}, sort_keys=True),
        )
        assert miss is None
    finally:
        conn.close()


def test_route_returns_typed_results_for_tested_generator(tmp_path: Path) -> None:
    """End-to-end: save typed rows, hit /api/tested-generators/{id}, and
    verify the legacy-shaped JSON is what the UI gets."""
    import sqlite3 as _sqlite3

    from fastapi.testclient import TestClient
    from regpoly_web.app import create_app
    from regpoly_web.config import Settings

    db_path = tmp_path / "route.db"
    settings = Settings(db_url=str(db_path), reload=False, pool_size=1)
    app = create_app(settings)
    with TestClient(app) as client:
        # Insert tested_generator + a component, then save typed results.
        conn = _sqlite3.connect(str(db_path))
        try:
            conn.execute(
                "INSERT INTO tested_generator(id, search_run_id, Lmax, k_g, J) "
                "VALUES (?, ?, ?, ?, ?)",
                (777, None, 32, 19937, 1),
            )
            conn.execute(
                "INSERT INTO tested_generator_component"
                "(tested_gen_id, component_index, family, L, k, "
                " all_params, tempering_params) "
                "VALUES (?, ?, ?, ?, ?, ?, ?)",
                (777, 0, "MTGen", 19937, 32,
                 json.dumps({}), json.dumps([])),
            )
            conn.commit()
            _save_one_of_each(conn, 777)
            conn.commit()
        finally:
            conn.close()

        r = client.get("/api/tested-generators/777")
        assert r.status_code == 200, r.text
        body = r.json()
        assert body["id"] == 777
        assert body["k_g"] == 19937
        assert len(body["results"]) == 3
        kinds = {res["test_type"] for res in body["results"]}
        assert kinds == {"equidistribution", "collision_free", "tuplets"}
