# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""
regpoly.analyses.pis — primary-invariant-subspace gap analysis for a
single generator (no tempering chain).

Phase 5.2 wraps the C++ PISCache so the web app's per-generator
analysis worker (`regpoly_web.tasks.analysis`) doesn't have to touch
_cpp directly.
"""

from __future__ import annotations

import time
from typing import Any

from regpoly_cpp import _regpoly_cpp as _cpp

from regpoly.core.combination import Combination
from regpoly.core.generator import Generator


def analyze_single_generator(gen: Generator) -> dict[str, Any]:
    """Compute the per-resolution gaps + summary for one Generator.

    Builds a single-component Combination around `gen` (no tempering),
    feeds it to the C++ PISCache, and returns:

        {
            "gaps":           [int]  # 1-indexed up to comb.L
            "se":             int    # sum(gaps)
            "elapsed":        float  # seconds spent in compute_all
            "char_poly_int":  int    # the BitVect-encoded char poly value
            "hamming_weight": int    # popcount(char_poly_int) + 1 (leading x^k)
            "k":              int
            "L":              int
        }

    The `char_poly_int` and `hamming_weight` are computed alongside so
    callers (typically the web's analysis worker) can persist them in
    one shot.
    """
    cpoly_bv = gen.char_poly()
    cpoly_int = int(cpoly_bv._val)
    # Char poly is monic of degree k; the leading x^k coefficient is
    # implicit (not stored in _val), so add it to the popcount.
    hw = bin(cpoly_int).count("1") + 1

    comb = Combination(J=1, Lmax=gen.L)
    comb.components[0].add_gen(gen)
    comb.reset()

    t0 = time.time()
    cache = _cpp.PISCache(
        [comb[0]._cpp_gen], [[]], comb.k_g, comb.L,
    )
    ecart = cache.compute_all()
    elapsed = time.time() - t0
    gaps = [int(ecart[v]) for v in range(1, comb.L + 1)]
    return {
        "gaps":           gaps,
        "se":             sum(gaps),
        "elapsed":        elapsed,
        "char_poly_int":  cpoly_int,
        "hamming_weight": hw,
        "k":              comb.k_g,
        "L":              comb.L,
    }


__all__ = ["analyze_single_generator"]
