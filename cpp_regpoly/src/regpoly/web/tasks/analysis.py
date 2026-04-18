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

import time
import traceback

import regpoly._regpoly_cpp as _cpp
from regpoly.combinaison import Combinaison
from regpoly.generateur import Generateur
from regpoly.web.database import json_dumps, json_loads, sync_connect


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
            gen = Generateur.create(family, L, **all_params)

            cpoly_bv = gen.char_poly()
            cpoly_int = int(cpoly_bv._val)
            cpoly_hex = hex(cpoly_int)
            # The characteristic polynomial is monic of degree k, so the
            # leading x^k coefficient (not stored in _val) always counts.
            hw = bin(cpoly_int).count("1") + 1

            comb = Combinaison(J=1, Lmax=L)
            comb.components[0].add_gen(gen)
            comb.reset()

            t0 = time.time()
            cache = _cpp.PISCache(
                [comb[0]._cpp_gen], [[]], comb.k_g, comb.L
            )
            ecart = cache.compute_all()
            # ecart is 1-indexed up to L
            gaps = [int(ecart[v]) for v in range(1, comb.L + 1)]
            se = sum(gaps)
            elapsed = time.time() - t0

            conn.execute(
                """
                UPDATE primitive_generator
                SET char_poly = ?, hamming_weight = ?,
                    pis_gaps = ?, pis_se = ?,
                    pis_elapsed = ?, pis_computed_at = datetime('now'),
                    pis_error = NULL
                WHERE id = ?
                """,
                (
                    cpoly_hex, hw,
                    json_dumps(gaps), se,
                    elapsed, gen_id,
                ),
            )
            conn.commit()

        except Exception as exc:
            msg = f"{exc.__class__.__name__}: {exc}\n{traceback.format_exc()}"
            conn.execute(
                """
                UPDATE primitive_generator
                SET pis_error = ?, pis_computed_at = datetime('now')
                WHERE id = ?
                """,
                (msg[:4000], gen_id),
            )
            conn.commit()
