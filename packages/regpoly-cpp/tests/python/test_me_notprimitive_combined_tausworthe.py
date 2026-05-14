# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Regression: ``test_me_notprimitive`` on combined Tausworthe gens.

LFSR258 (5-component combined Tausworthe, L'Ecuyer 1999) used to fail
with::

    RuntimeError: test_me_notprimitive: every component projects to
    zero in V (check χ recovery)

Root cause: ``polychar_comb`` returned the *raw recurrence trinomial*
for each TauswortheGen, but the Tausworthe transitions by T^s per
output (s_j ∈ {10, 5, 29, 23, 8} here). The LCM of that χ_pc with the
BM-recovered χ contained both sets of irreducibles; ``select_phi``
could pick a spurious factor (one in T's poly but not T^s's), at
which point χ_ψ = χ/φ still annihilated every component.

Fix: ``test_me_notprimitive`` now relies solely on ``recover_char_poly``
(Krylov BM on the actual state evolution). This test pins the fix.
"""

from __future__ import annotations

import pytest

from regpoly_cpp import _regpoly_cpp as cpp
from regpoly.core.combination import Combination
from regpoly.core.combination_build import build_combinaison_inputs


def _build_combination(components: list[dict], Lmax: int):
    gen_lists, temperings = build_combinaison_inputs(components, Lmax)
    comb = Combination.CreateFromFiles(gen_lists, Lmax, temperings)
    comb.reset()
    gens = [comb[j]._cpp_gen for j in range(comb.J)]
    trans = [
        [t._cpp_trans for t in cp.trans if hasattr(t, "_cpp_trans")]
        for cp in comb.components
    ]
    return comb, gens, trans


def test_lfsr258_no_longer_raises():
    """LFSR258 — 5 Tausworthe components, k_j = {63,55,52,47,41},
    s_j = {10,5,29,23,8}, kg=258, L=64. Pre-fix: every component
    projected to zero in V. Post-fix: returns a valid result."""
    components = [
        {"family": "Tausworthe", "L": 64,
         "params": {"k": 63, "nb_terms": 3, "poly": [0, 1, 63],
                    "quicktaus": True, "s": 10}},
        {"family": "Tausworthe", "L": 64,
         "params": {"k": 55, "nb_terms": 3, "poly": [0, 24, 55],
                    "quicktaus": True, "s": 5}},
        {"family": "Tausworthe", "L": 64,
         "params": {"k": 52, "nb_terms": 3, "poly": [0, 3, 52],
                    "quicktaus": True, "s": 29}},
        {"family": "Tausworthe", "L": 64,
         "params": {"k": 47, "nb_terms": 3, "poly": [0, 5, 47],
                    "quicktaus": True, "s": 23}},
        {"family": "Tausworthe", "L": 64,
         "params": {"k": 41, "nb_terms": 3, "poly": [0, 3, 41],
                    "quicktaus": True, "s": 8}},
    ]
    comb, gens, trans = _build_combination(components, 64)
    assert comb.k_g == 258, comb.k_g

    res = cpp.test_me_notprimitive(
        gens, trans, comb.k_g, comb.L, comb.L,
        [0] * (comb.L + 1), 0,
    )
    # Pre-fix this raised. Post-fix any well-formed dict is fine — the
    # specific se value is an algorithmic detail and may shift with
    # future PIS refinements.
    assert isinstance(res, dict)
    assert "ecart" in res and "se" in res
    assert isinstance(res["se"], int)
    assert len(res["ecart"]) == comb.L + 1
    # ecart values must be finite (>=0) or the INT_MAX sentinel for
    # un-verified levels. Anything negative or otherwise pathological
    # would indicate a corrupted return.
    for v in res["ecart"]:
        assert v >= 0, v


def test_two_component_tausworthe_smaller():
    """Minimal repro: 2 Tausworthe components with different s. Pre-fix
    this would have hit the same χ_ψ bug whenever the picked φ was a
    spurious raw-recurrence factor. Post-fix it must produce a result.
    """
    components = [
        {"family": "Tausworthe", "L": 32,
         "params": {"k": 31, "nb_terms": 3, "poly": [0, 6, 31],
                    "quicktaus": True, "s": 18}},
        {"family": "Tausworthe", "L": 32,
         "params": {"k": 29, "nb_terms": 3, "poly": [0, 2, 29],
                    "quicktaus": True, "s": 2}},
    ]
    comb, gens, trans = _build_combination(components, 32)
    res = cpp.test_me_notprimitive(
        gens, trans, comb.k_g, comb.L, comb.L,
        [0] * (comb.L + 1), 0,
    )
    assert isinstance(res, dict) and "se" in res
