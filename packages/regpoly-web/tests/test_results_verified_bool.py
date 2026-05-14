# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Regression: `verified` is a PG BOOLEAN, not a smallint.

The Run-a-test card on /library/{paper}/{gen} was failing with::

    column "verified" is of type boolean but expression is of type smallint

because results.py wrote ``1 if … else 0`` — fine under SQLite, but
psycopg sends Python ``int`` as smallint and PG refuses the implicit
cast to BOOLEAN. We now pass ``bool(...)``; this test pins the typed
write against a real PG and asserts the column stores a Python bool.
"""

from __future__ import annotations

import json

import psycopg


def test_save_equidistribution_writes_real_boolean(
    seeded_db_url: str, seeded_db: str,
) -> None:
    """End-to-end through ``save_typed_result``: the equidistribution
    insert must accept ``verified=True`` and the column must round-trip
    as a real Python bool (not 0/1)."""
    from regpoly.analyses.equidistribution_results import EquidistributionResults
    from regpoly_web.database import SyncConnShim
    from regpoly_web.results import save_typed_result

    L = 8
    ecart = [0] * (L + 1)
    ecart[3] = 5
    res = EquidistributionResults(
        L=L, ecart=ecart, psi12=[True] * (L + 1),
        se=5, verified=True, mse=0, meverif=True,
        delta=[0] * (L + 1),
    )
    raw = psycopg.connect(seeded_db_url)
    try:
        # Test config blob mirrors what the route sends; the helper
        # expects a JSON string (legacy SQLite path).
        conn = SyncConnShim(raw)
        rid = save_typed_result(
            conn, 4242,
            json.dumps({"type": "equidistribution", "L": L}),
            res, kg=19937, L=L, elapsed_seconds=0.5,
        )
        raw.commit()
        assert rid is not None

        # Read back and assert the PG-typed shape: real bool, not int.
        with raw.cursor() as cur:
            cur.execute(
                "SELECT verified FROM equidistribution_result WHERE id = %s",
                (rid,),
            )
            (verified,) = cur.fetchone()
        assert verified is True
        assert isinstance(verified, bool), type(verified)
    finally:
        raw.close()
