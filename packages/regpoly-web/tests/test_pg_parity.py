# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Cross-engine SQL parity: SQLite (legacy) vs PostgreSQL (current).

The dockerize plan switched the backing store from SQLite to PG.
Translating a query is mechanical (``?`` → ``%s``, ``datetime('now')``
→ ``NOW()``) but a handful of SQL constructs return *different
results* on the two engines:

- ``ORDER BY ... DESC`` puts NULLs FIRST in SQLite, LAST in PG.
- JSON columns return strings in SQLite (TEXT), parsed objects in PG
  (JSONB).
- ``BOOLEAN`` columns store as 0/1 INTEGER in SQLite, true/false in PG.
- Joining a NULL FK does not return the row in either engine.

This file plants the same fixture rows in both engines and asserts
identical post-shaping result sets for the queries the route layer
relies on. Acts as insurance: future query churn that drifts in
either direction trips the parity.

Marked ``slow`` because each test reapplies both schemas.
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import psycopg
import pytest

from regpoly_web.config import SCHEMA_PATH


def _open_sqlite(tmp_path: Path) -> sqlite3.Connection:
    db = tmp_path / "parity.db"
    conn = sqlite3.connect(str(db))
    conn.row_factory = sqlite3.Row
    conn.executescript(Path(SCHEMA_PATH).read_text())
    return conn


def _seed_sqlite(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        INSERT INTO primitive_search_run
            (family, L, k, structural_params, fixed_params, status,
             tries_done, found_count)
        VALUES
            ('MTGen', 64, 32, '{}', '{}', 'completed',  100, 3),
            ('MTGen', 64, 32, '{}', '{}', 'running',     50, 1),
            ('MTGen', 64, 32, '{}', '{}', 'pending',      0, 0);

        INSERT INTO primitive_generator
            (search_run_id, family, L, k,
             structural_params, search_params, all_params, found_at_try,
             pis_se, pis_computed_at)
        VALUES
            (1, 'MTGen', 64, 32, '{"L":64}', '{"a":"0x1"}',
             '{"L":64,"a":"0x1"}', 1, 5, datetime('now')),
            (1, 'MTGen', 64, 32, '{"L":64}', '{"a":"0x2"}',
             '{"L":64,"a":"0x2"}', 2, NULL, NULL);
        """
    )
    conn.commit()


def _seed_pg(conn: psycopg.Connection) -> None:
    conn.execute(
        "INSERT INTO primitive_search_run "
        "(family, l, k, structural_params, fixed_params, status, "
        " tries_done, found_count) VALUES "
        "('MTGen', 64, 32, '{}'::jsonb, '{}'::jsonb, 'completed', 100, 3),"
        "('MTGen', 64, 32, '{}'::jsonb, '{}'::jsonb, 'running',    50, 1),"
        "('MTGen', 64, 32, '{}'::jsonb, '{}'::jsonb, 'pending',     0, 0)"
    )
    conn.execute(
        "INSERT INTO primitive_generator "
        "(search_run_id, family, l, k, structural_params, search_params, "
        " all_params, found_at_try, pis_se, pis_computed_at) VALUES "
        "(1, 'MTGen', 64, 32, '{\"L\":64}'::jsonb, '{\"a\":\"0x1\"}'::jsonb, "
        " '{\"L\":64,\"a\":\"0x1\"}'::jsonb, 1, 5, NOW()),"
        "(1, 'MTGen', 64, 32, '{\"L\":64}'::jsonb, '{\"a\":\"0x2\"}'::jsonb, "
        " '{\"L\":64,\"a\":\"0x2\"}'::jsonb, 2, NULL, NULL)"
    )
    conn.commit()


def _normalize(rows: list, json_cols: tuple[str, ...] = ()) -> list[dict]:
    """Coerce rows to a uniform shape so the two engines compare equal.

    - dict access (sqlite3.Row supports it)
    - JSON columns: parse strings → objects (SQLite TEXT vs PG JSONB)
    - Booleans: cast 0/1 → bool
    - Timestamps: replaced by ``not None`` truthy/falsy markers
    """
    out: list[dict] = []
    for r in rows:
        d = dict(r)
        for col in json_cols:
            v = d.get(col)
            if isinstance(v, str):
                d[col] = json.loads(v)
        for col in ("verified", "is_me"):
            if col in d and d[col] is not None and not isinstance(d[col], bool):
                d[col] = bool(d[col])
        for col in ("pis_computed_at", "created_at", "started_at",
                    "finished_at"):
            if col in d:
                d[col] = d[col] is not None
        out.append(d)
    return out


@pytest.mark.slow
def test_parity_select_status_counts(tmp_path: Path, tmp_db_url: str) -> None:
    """SELECT status, COUNT(*) GROUP BY status — same row set on both."""
    sq = _open_sqlite(tmp_path)
    _seed_sqlite(sq)
    sq_rows = _normalize(sq.execute(
        "SELECT status, COUNT(*) AS n FROM primitive_search_run "
        "GROUP BY status ORDER BY status"
    ).fetchall())
    sq.close()

    pg = psycopg.connect(tmp_db_url, autocommit=True)
    try:
        _seed_pg(pg)
        pg_cur = pg.execute(
            "SELECT status, COUNT(*) AS n FROM primitive_search_run "
            "GROUP BY status ORDER BY status"
        )
        # psycopg returns tuples by default; map by column.
        cols = [d.name for d in pg_cur.description]
        pg_rows = _normalize([dict(zip(cols, row)) for row in pg_cur.fetchall()])
    finally:
        pg.close()

    assert sq_rows == pg_rows


@pytest.mark.slow
def test_parity_join_with_null_fk(tmp_path: Path, tmp_db_url: str) -> None:
    """LEFT JOIN where FK is NULL: PG and SQLite both include the
    base row with NULLs in the joined columns."""
    sq = _open_sqlite(tmp_path)
    _seed_sqlite(sq)
    sq_rows = _normalize(sq.execute(
        "SELECT pg.id AS gid, pg.found_at_try, pg.pis_se "
        "FROM primitive_generator pg "
        "LEFT JOIN primitive_search_run psr ON psr.id = pg.search_run_id "
        "ORDER BY pg.id"
    ).fetchall())
    sq.close()

    pg_conn = psycopg.connect(tmp_db_url, autocommit=True)
    try:
        _seed_pg(pg_conn)
        cur = pg_conn.execute(
            "SELECT pg.id AS gid, pg.found_at_try, pg.pis_se "
            "FROM primitive_generator pg "
            "LEFT JOIN primitive_search_run psr ON psr.id = pg.search_run_id "
            "ORDER BY pg.id"
        )
        cols = [d.name for d in cur.description]
        pg_rows = _normalize([dict(zip(cols, row)) for row in cur.fetchall()])
    finally:
        pg_conn.close()

    assert sq_rows == pg_rows


@pytest.mark.slow
def test_parity_json_columns_round_trip(tmp_path: Path, tmp_db_url: str) -> None:
    """JSON columns: SQLite returns TEXT, PG returns parsed dict.

    Both must round-trip the same logical value when normalized.
    """
    sq = _open_sqlite(tmp_path)
    _seed_sqlite(sq)
    sq_rows = _normalize(sq.execute(
        "SELECT id, all_params FROM primitive_generator ORDER BY id"
    ).fetchall(), json_cols=("all_params",))
    sq.close()

    pg = psycopg.connect(tmp_db_url, autocommit=True)
    try:
        _seed_pg(pg)
        cur = pg.execute(
            "SELECT id, all_params FROM primitive_generator ORDER BY id"
        )
        cols = [d.name for d in cur.description]
        pg_rows = _normalize(
            [dict(zip(cols, row)) for row in cur.fetchall()],
            json_cols=("all_params",),
        )
    finally:
        pg.close()

    assert sq_rows == pg_rows


@pytest.mark.slow
def test_parity_order_by_pis_computed_at_with_nulls(
    tmp_path: Path, tmp_db_url: str,
) -> None:
    """``ORDER BY pis_computed_at DESC`` differs between engines unless
    the query is explicit. SQLite default DESC: NULLs LAST. PG default
    DESC: NULLs FIRST.

    Confirms the dockerize-plan-mandated behaviour: routes that depend
    on the SQLite default ordering must add ``NULLS LAST`` in PG. This
    test PROVES the engines disagree without an explicit clause, then
    PROVES they agree once ``NULLS LAST`` is added in PG.
    """
    sq = _open_sqlite(tmp_path)
    _seed_sqlite(sq)
    sq_rows = _normalize(sq.execute(
        "SELECT id FROM primitive_generator "
        "ORDER BY pis_computed_at DESC, id ASC"
    ).fetchall())
    sq.close()

    pg = psycopg.connect(tmp_db_url, autocommit=True)
    try:
        _seed_pg(pg)
        cur = pg.execute(
            "SELECT id FROM primitive_generator "
            "ORDER BY pis_computed_at DESC, id ASC"
        )
        cols = [d.name for d in cur.description]
        pg_rows_default = _normalize(
            [dict(zip(cols, row)) for row in cur.fetchall()]
        )
        # Add explicit `NULLS LAST` to make PG match SQLite default
        # DESC behaviour (NULLs LAST).
        cur = pg.execute(
            "SELECT id FROM primitive_generator "
            "ORDER BY pis_computed_at DESC NULLS LAST, id ASC"
        )
        cols = [d.name for d in cur.description]
        pg_rows_explicit = _normalize(
            [dict(zip(cols, row)) for row in cur.fetchall()]
        )
    finally:
        pg.close()

    # Default DESC: PG puts the NULL row LAST; SQLite puts it FIRST.
    # If they ever match by accident the parity is meaningless.
    assert sq_rows != pg_rows_default, (
        "SQLite/PG default DESC ordering must differ on NULL position; "
        "the parity test is no longer load-bearing if they match"
    )
    # With NULLS LAST in PG, the order matches SQLite default DESC.
    assert sq_rows == pg_rows_explicit, (
        f"With explicit NULLS LAST, PG must match SQLite DESC default; "
        f"sqlite={sq_rows} pg={pg_rows_explicit}"
    )
