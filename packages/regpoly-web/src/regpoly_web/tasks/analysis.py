# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Per-generator analysis worker.

Runs inside a child process (via ProcessPoolExecutor).  For a single
primitive generator, computes:

  - the characteristic polynomial (as a hex string)
  - its Hamming weight (number of 1-bits)
  - the equidistribution dimension gaps for v = 1..L using the PIS
    (StackBase) method — gap_v = floor(k/v) - t_v
  - se = sum of gaps

Results are written back to the primitive_generator row.  Errors are
recorded in `pis_error` so the scheduler knows not to retry the same
row in a tight loop.
"""

from __future__ import annotations

import traceback

from regpoly.analyses.pis import analyze_single_generator
from regpoly.core.generator import Generator

from regpoly_web.database import json_dumps, json_loads, sync_connect


def analyze_generator(db_path: str, gen_id: int) -> None:
    """Entry point run by a pool worker for a single generator."""
    with sync_connect(db_path) as conn:
        row = conn.execute(
            "SELECT family, L, k, all_params, pis_computed_at "
            "FROM primitive_generator WHERE id = ?",
            (gen_id,),
        ).fetchone()
        if row is None or row["pis_computed_at"] is not None:
            return

        family = row["family"]
        L = row["L"]
        k = row["k"]
        all_params = json_loads(row["all_params"]) or {}

        try:
            gen = Generator.create(family, L, **all_params)
            res = analyze_single_generator(gen)
            cpoly_hex = hex(res["char_poly_int"])

            conn.execute(
                """
                UPDATE primitive_generator
                SET char_poly = ?, hamming_weight = ?,
                    pis_gaps = ?, pis_se = ?,
                    pis_elapsed = ?, pis_computed_at = NOW(),
                    pis_error = NULL
                WHERE id = ?
                """,
                (
                    cpoly_hex, res["hamming_weight"],
                    json_dumps(res["gaps"]), res["se"],
                    res["elapsed"], gen_id,
                ),
            )
            conn.commit()

        except Exception as exc:
            msg = f"{exc.__class__.__name__}: {exc}\n{traceback.format_exc()}"
            conn.execute(
                """
                UPDATE primitive_generator
                SET pis_error = ?, pis_computed_at = NOW()
                WHERE id = ?
                """,
                (msg[:4000], gen_id),
            )
            conn.commit()
