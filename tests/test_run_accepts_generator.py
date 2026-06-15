# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""
`Test(...).run(gen)` equivalence test.

After the Phase 1-3 migration, every `AbstractTest.run(...)` accepts a
bare `Generator` (Python wrapper or `_cpp.CombinedF2LinearSource`).
Single-component results should match what
`make_combined(gen)` produces; multi-component results from
`make_combined(g1, g2)` are pinned by direct comparison against running
the same test on a J=2 CombinedF2LinearSource built from a fresh pair of
the same generators.

Field-by-field comparison because none of the four `*Results` classes
define `__eq__`.
"""

from __future__ import annotations

import sys

import pytest

import regpoly._regpoly_cpp as _cpp
from regpoly import make_combined
from regpoly.analyses.collision_free_test import CollisionFreeTest
from regpoly.analyses.equidistribution_test import EquidistributionTest
from regpoly.analyses.tuplets_test import TupletsTest
from regpoly.analyses.tvalue_test import TValueTest
from regpoly.core.generator import Generator


def _assert_eq_results(r1, r2, *fields):
    for f in fields:
        v1, v2 = getattr(r1, f), getattr(r2, f)
        assert v1 == v2, f"field {f!r}: {v1!r} != {v2!r}"


# ── Single-Generator equivalence (run(gen) vs run(make_combined(gen))) ──

def _make_mt_gen():
    # Small Tausworthe LFSR — known primitive (trinomial x^31 + x^6 + 1).
    return Generator.create(
        "TauswortheGen", L=32, k=31, nb_terms=3,
        poly=[0, 6, 31], s=18, quicktaus=True,
    )


def test_equidistribution_accepts_bare_generator():
    g = _make_mt_gen()
    t = EquidistributionTest(L=32, delta=[sys.maxsize] * 33, mse=sys.maxsize,
                             method=None)
    r_gen = t.run(g)
    r_comb = t.run(make_combined(g, Lmax=32))
    _assert_eq_results(r_gen, r_comb, "ecart", "se", "verified")


def test_collision_free_accepts_bare_generator():
    g = _make_mt_gen()
    t = CollisionFreeTest(msecf=sys.maxsize)
    r_gen = t.run(g)
    r_comb = t.run(make_combined(g, Lmax=32))
    _assert_eq_results(r_gen, r_comb, "ecart_cf", "secf", "verified")


def test_collision_free_me_results_kwarg_preserved():
    g = _make_mt_gen()
    eqt = EquidistributionTest(L=32, delta=[sys.maxsize] * 33,
                               mse=sys.maxsize, method=None)
    me_res = eqt.run(g)
    cf = CollisionFreeTest(msecf=sys.maxsize)
    r_gen = cf.run(g, me_results=me_res)
    r_comb = cf.run(make_combined(g, Lmax=32), me_results=me_res)
    _assert_eq_results(r_gen, r_comb, "ecart_cf", "secf", "verified")


def test_tuplets_accepts_bare_generator():
    g = _make_mt_gen()
    t = TupletsTest(tupletsverif=True, d=2, s=[0, 10, 5],
                    mDD=0.0, testtype=0)
    r_gen = t.run(g)
    r_comb = t.run(make_combined(g, Lmax=32))
    _assert_eq_results(r_gen, r_comb,
                       "tupd", "DELTA", "firstpart_max",
                       "secondpart_max", "verified")


def test_tvalue_accepts_bare_generator():
    g = Generator.create("SobolNet", L=4, s_max=4)
    t = TValueTest(s_max=4, max_t_sum=sys.maxsize)
    r_gen = t.run(g)
    r_comb = t.run(make_combined(g, Lmax=4))
    _assert_eq_results(r_gen, r_comb, "tvals", "se", "verified")


# ── Multi-component equivalence (make_combined(g1, g2) vs fresh rebuild) ─


def _make_two_taus():
    g1 = Generator.create("TauswortheGen", L=32, k=31, nb_terms=3,
                          poly=[0, 6, 31], s=18, quicktaus=True)
    g2 = Generator.create("TauswortheGen", L=32, k=29, nb_terms=3,
                          poly=[0, 2, 29], s=2, quicktaus=True)
    return g1, g2


def test_equidistribution_J2_make_combined_is_deterministic():
    g1, g2 = _make_two_taus()
    combined_a = make_combined(g1, g2, Lmax=32)
    combined_b = make_combined(g1, g2, Lmax=32)
    t = EquidistributionTest(L=32, delta=[sys.maxsize] * 33,
                             mse=sys.maxsize, method=None)
    r_a = t.run(combined_a)
    r_b = t.run(combined_b)
    _assert_eq_results(r_a, r_b, "ecart", "se", "verified")


def test_tvalue_accepts_raw_cpp_object_via_duck_typing():
    """A raw `_cpp.SobolNet` (no Python Generator wrapping) — current
    `_to_cpp_gen` policy accepts it via duck-typing (it has
    `.k()` / `.L()` / `.get_output()`)."""
    raw = _cpp.SobolNet(4, 4)
    t = TValueTest(s_max=4, max_t_sum=sys.maxsize)
    r = t.run(raw)
    assert r.verified


def test_run_raises_on_garbage_input():
    t = EquidistributionTest(L=32, delta=[sys.maxsize] * 33, mse=sys.maxsize)
    with pytest.raises(TypeError, match="neither a regpoly.Generator"):
        t.run("not a generator")
    with pytest.raises(TypeError, match="neither a regpoly.Generator"):
        t.run(42)
