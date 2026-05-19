# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Primary-invariant-subspace gap analysis for a single generator.

Thin Python shim over `regpoly_cpp.PISCache` so the web app's
per-generator analysis worker (`regpoly_web.tasks.analysis`) doesn't
have to touch `_regpoly_cpp` directly. The output dict is the canonical
shape that downstream persistence + UI code expects.
"""

from __future__ import annotations

import time
from typing import Any

from regpoly_cpp import _regpoly_cpp as _cpp

from regpoly.core.combination import Combination
from regpoly.core.generator import Generator


def analyze_single_generator(gen: Generator) -> dict[str, Any]:
    """Compute the per-resolution dimension defect of one F₂-linear generator.

    Wraps the generator in a single-component
    :class:`regpoly.core.combination.Combination` (no tempering),
    runs the C++ PISCache, and packages the result as a flat dict so
    downstream code (the web analysis worker, the
    `regpoly` CLI's `show` mode) can persist it in one shot.

    Parameters
    ----------
    gen
        A constructed :class:`regpoly.core.generator.Generator`.
        The generator does not need to be at any particular state —
        the wrapping `Combination.reset()` is called internally.

    Returns
    -------
    dict
        A dict with the following keys:

        - ``gaps`` (`list[int]`): per-resolution defect, indexed 1..L.
        - ``se`` (`int`): total dimension defect (sum of ``gaps``).
        - ``elapsed`` (`float`): wall-clock seconds spent in
          `PISCache.compute_all`.
        - ``char_poly_int`` (`int`): integer encoding of the
          characteristic polynomial (leading ``x^k`` coefficient is
          implicit).
        - ``hamming_weight`` (`int`): population count of
          ``char_poly_int`` plus 1 (for the implicit leading bit).
        - ``k`` (`int`): total state size in bits.
        - ``L`` (`int`): output resolution in bits.

    See Also
    --------
    :cpp:class:`regpoly::core::PISCache` : the C++ implementation this wraps.

    Examples
    --------
    >>> from regpoly import Generator
    >>> from regpoly.analyses.pis import analyze_single_generator
    >>> g = Generator.create(             # doctest: +SKIP
    ...     "MTGen", L=32, w=32, r=624, m=397, p=31,
    ...     a=0x9908b0df,
    ... )
    >>> res = analyze_single_generator(g)  # doctest: +SKIP
    >>> res["se"]                          # doctest: +SKIP
    0
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
