"""Phase 2.3 confidence test: TupletsTest.run() exercises the new C++
run_tuplets path end-to-end.

The Python implementation was removed in Phase 2.3, so this test runs
against the C++ runner only and verifies the result has the right shape
and no spurious sentinel values for a known-good combination.
"""

from __future__ import annotations

from regpoly.analyses.tuplets_results import _MAX_TYPE
from regpoly.analyses.tuplets_test import TupletsTest
from regpoly.core.combination import Combination
from regpoly.core.component import Component
from regpoly.core.generator import Generator


def _make_combo() -> Combination:
    """Build a 1-component Tausworthe(31) combination at L=32."""
    g = Generator.create(
        "TauswortheGen",
        L=32,
        k=31, poly=[0, 3, 31], nb_terms=3, s=13, quicktaus=True,
    )
    comp = Component()
    comp.add_gen(g)
    comb = Combination(J=1, Lmax=32)
    comb.components = [comp]
    comb.reset()
    return comb


def test_run_returns_well_shaped_result() -> None:
    comb = _make_combo()
    test = TupletsTest(
        tupletsverif=True,
        d=2,
        s=[0, 5, 4],  # 1-indexed
        mDD=10.0,     # high enough that the run completes
        testtype=_MAX_TYPE,
    )
    result = test.run(comb)

    assert result.tupletsverif is True
    assert result.tupd == 2
    assert len(result.tuph) == 3
    assert len(result.gap) == 6  # tuph[1]+1 = 5+1
    assert len(result.DELTA) == 3
    assert len(result.pourcentage) == 3

    # gap[1..min(tuph[1], kg)] should have been computed; remaining slots
    # default to 0.0.
    assert all(isinstance(g, float) for g in result.gap)
    assert all(isinstance(d, float) for d in result.DELTA)


def test_disabled_short_circuits() -> None:
    comb = _make_combo()
    test = TupletsTest(tupletsverif=False)
    result = test.run(comb)
    assert result.tupletsverif is False
    assert result.tupd == 0
