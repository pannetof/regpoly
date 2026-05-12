# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Worker-side tests for the primitive-search SE post-filter.

Drive the real ``run_primitive_search`` worker function against an
ephemeral PG, but monkey-patch ``Generator.create``,
``is_full_period``, and ``analyze_single_generator`` so the test
doesn't depend on the C++ search loop. Verifies:

- a full-period candidate with SE > max_se increments ``rejected_count``
  and does NOT insert a row;
- a full-period candidate with SE ≤ max_se inserts a row with
  ``pis_*`` columns populated;
- ``rejected_count`` survives pause/resume.
"""

from __future__ import annotations

import json

import psycopg
import pytest
from regpoly_web.tasks import primitive as primitive_task


class _FakeGen:
    """Minimal stand-in for ``regpoly.core.generator.Generator``.

    Carries the ``params`` and ``k`` attributes the worker reads, plus a
    canned ``is_full_period`` answer.
    """

    def __init__(self, full_period: bool, params: dict, k: int = 31):
        self._full_period = full_period
        self.params = params
        self.k = k

    def is_full_period(self) -> bool:
        return self._full_period


def _insert_se_run(db_url: str, *, max_se: int,
                   rejected_count: int = 0,
                   status: str = "pending") -> int:
    """Insert a minimal primitive_search_run with max_se set."""
    conn = psycopg.connect(db_url, autocommit=True)
    try:
        cur = conn.execute(
            """
            INSERT INTO primitive_search_run
                (family, l, k, structural_params, fixed_params,
                 max_tries, max_se, rejected_count, status)
            VALUES ('TauswortheGen', 64, 31, %s, %s,
                    %s, %s, %s, %s)
            RETURNING id
            """,
            (json.dumps({"k": 31, "s": 3, "q": 13}),
             json.dumps({}), 1, max_se, rejected_count, status),
        )
        return int(cur.fetchone()[0])
    finally:
        conn.close()


def _read_row(db_url: str, run_id: int) -> dict:
    conn = psycopg.connect(db_url, autocommit=True)
    try:
        cur = conn.execute(
            "SELECT status, tries_done, found_count, rejected_count "
            "FROM primitive_search_run WHERE id = %s",
            (run_id,),
        )
        cols = [d.name for d in cur.description]
        row = cur.fetchone()
    finally:
        conn.close()
    assert row is not None
    return dict(zip(cols, row))


def _count_generators(db_url: str, run_id: int) -> int:
    conn = psycopg.connect(db_url, autocommit=True)
    try:
        row = conn.execute(
            "SELECT COUNT(*) FROM primitive_generator WHERE search_run_id=%s",
            (run_id,),
        ).fetchone()
    finally:
        conn.close()
    return int(row[0])


def _read_generator(db_url: str, run_id: int) -> dict | None:
    conn = psycopg.connect(db_url, autocommit=True)
    try:
        cur = conn.execute(
            "SELECT pis_se, pis_gaps, pis_computed_at, char_poly, "
            "hamming_weight FROM primitive_generator "
            "WHERE search_run_id=%s LIMIT 1",
            (run_id,),
        )
        cols = [d.name for d in cur.description]
        row = cur.fetchone()
    finally:
        conn.close()
    if row is None:
        return None
    return dict(zip(cols, row))


def _patch_worker(monkeypatch, *, gen_full_period: bool,
                  analysis: dict | None) -> None:
    """Stub the worker's external calls so it runs deterministically.

    ``Generator.create`` returns a single ``_FakeGen`` (max_tries=1 in
    the row keeps the loop to one iteration). ``analyze_single_generator``
    returns ``analysis``. ``introspection.get_gen_param_specs`` is also
    stubbed because the worker calls it up front.
    """
    fake_gen = _FakeGen(gen_full_period, params={"k": 31, "s": 3, "q": 13})

    def _create(family, L, **kw):
        return fake_gen

    monkeypatch.setattr(
        "regpoly.core.generator.Generator.create",
        classmethod(lambda cls, *a, **kw: fake_gen),
    )
    monkeypatch.setattr(
        primitive_task._introspection,
        "get_gen_param_specs",
        lambda family: [],
    )
    # ``analyze_single_generator`` is lazy-imported inside _run_random.
    # Patch the source module so the import resolves to our stub.
    monkeypatch.setattr(
        "regpoly.analyses.pis.analyze_single_generator",
        lambda gen: analysis,
    )


@pytest.fixture
def tmp_db_url(tmp_db_url):
    """Re-export the conftest fixture for clarity."""
    return tmp_db_url


def test_reject_path_increments_counter_no_insert(
    tmp_db_url: str, monkeypatch
) -> None:
    """A full-period candidate with SE > max_se → rejected_count++,
    no row in primitive_generator."""
    run_id = _insert_se_run(tmp_db_url, max_se=0)

    _patch_worker(
        monkeypatch,
        gen_full_period=True,
        analysis={
            "gaps": [1, 2],
            "se": 3,                # > max_se=0
            "elapsed": 0.01,
            "char_poly_int": 0xDEAD,
            "hamming_weight": 7,
            "k": 31, "L": 64,
        },
    )

    primitive_task.run_primitive_search(tmp_db_url, run_id)

    row = _read_row(tmp_db_url, run_id)
    assert row["rejected_count"] == 1, row
    assert row["found_count"] == 0, row
    assert _count_generators(tmp_db_url, run_id) == 0


def test_accept_path_populates_pis_columns(
    tmp_db_url: str, monkeypatch
) -> None:
    """SE ≤ max_se → row inserted with pis_se / pis_gaps / pis_computed_at
    populated, so the async analysis worker skips it."""
    run_id = _insert_se_run(tmp_db_url, max_se=10)

    _patch_worker(
        monkeypatch,
        gen_full_period=True,
        analysis={
            "gaps": [0, 0, 0],
            "se": 0,                # ≤ max_se=10
            "elapsed": 0.02,
            "char_poly_int": 0xBEEF,
            "hamming_weight": 5,
            "k": 31, "L": 64,
        },
    )

    primitive_task.run_primitive_search(tmp_db_url, run_id)

    row = _read_row(tmp_db_url, run_id)
    assert row["rejected_count"] == 0, row
    assert row["found_count"] == 1, row

    gen_row = _read_generator(tmp_db_url, run_id)
    assert gen_row is not None
    assert gen_row["pis_se"] == 0
    assert gen_row["pis_computed_at"] is not None
    assert gen_row["char_poly"] == hex(0xBEEF)
    assert gen_row["hamming_weight"] == 5
    # pis_gaps is JSONB; psycopg returns it as a Python list already.
    assert gen_row["pis_gaps"] == [0, 0, 0]


def test_rejected_count_survives_resume(
    tmp_db_url: str, monkeypatch
) -> None:
    """Insert a paused run with rejected_count=5, then run the worker.
    Even with a fresh reject this iteration, the persisted total must
    not reset to 1 — the worker must read the prior value first."""
    run_id = _insert_se_run(
        tmp_db_url, max_se=0, rejected_count=5, status="pending"
    )

    _patch_worker(
        monkeypatch,
        gen_full_period=True,
        analysis={
            "gaps": [1], "se": 1,   # > max_se=0
            "elapsed": 0.0,
            "char_poly_int": 1, "hamming_weight": 2,
            "k": 31, "L": 64,
        },
    )

    primitive_task.run_primitive_search(tmp_db_url, run_id)

    row = _read_row(tmp_db_url, run_id)
    # Prior 5 + this iteration's 1 reject.
    assert row["rejected_count"] == 6, row
