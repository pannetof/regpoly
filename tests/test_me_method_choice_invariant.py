# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Lock in agreement between the two general ME methods
(``matricial`` and ``notprimitive``) on the paper-reproduction cases
from Bhuvaneswari & Bhattacharjee 2026.

For every combined Cellular-Automata test point listed in the paper,
both methods must return the same ME boolean verdict. (They may
disagree on the ecart_sum magnitude — the matricial and BM-based
notprimitive paths track different intermediate quantities — but the
ME/NOT-ME boolean must agree.)

``simd_notprimitive`` is intentionally excluded: the kernel rejects
multi-component generators for that method, so it is not comparable on
these combined-CA cases.

This is a regression guard. If a future kernel change makes the two
methods disagree on the boolean ME verdict, this test surfaces the
mismatch in CI rather than during an ad-hoc rerun of the (deleted)
``check_me_method_choice`` script.
"""

from __future__ import annotations

import pytest

from regpoly import make_combined
from regpoly.analyses.equidistribution_test import (
    METHOD_MATRICIAL,
    METHOD_NOTPRIMITIVE,
    EquidistributionTest,
)
from regpoly.core.generator import Generator

_METHODS = (
    ("matricial", METHOD_MATRICIAL),
    ("notprimitive", METHOD_NOTPRIMITIVE),
)


def _make_comb(k1, p1, k2, p2, s):
    g1 = Generator.create("CellularAutomataGen", L=min(k1, 64),
                          k=k1, rule150_positions=p1, s=s)
    g2 = Generator.create("CellularAutomataGen", L=min(k2, 64),
                          k=k2, rule150_positions=p2, s=s)
    return make_combined(g1, g2, Lmax=64)


def _run_method(k1, p1, k2, p2, s, method_code):
    comb = _make_comb(k1, p1, k2, p2, s)
    test = EquidistributionTest(L=64, delta=[10**9] * 65,
                                mse=10**9, method=method_code)
    res = test.run(comb)
    return res.is_me()


_PAPER_CASES = [
    pytest.param(31, [10], 32, [0, 14], 7, id="Table9_31x32_s7"),
    pytest.param(31, [10], 32, [0, 14], 8, id="Table10_31x32_s8"),
    pytest.param(31, [10], 40, [7],     8, id="Table12_31x40_s8"),
    pytest.param(31, [10], 32, [0, 14], 1, id="Table3_31x32_s1"),
    pytest.param(37, [8],  42, [18],    8, id="Table5_37x42_s8"),
    pytest.param(31, [10], 33, [0],     8, id="Table6_31x33_s8"),
]


@pytest.mark.slow
@pytest.mark.parametrize(("k1", "p1", "k2", "p2", "s"), _PAPER_CASES)
def test_matricial_and_notprimitive_agree_on_me_verdict(k1, p1, k2, p2, s):
    """The matricial and notprimitive methods agree on the ME boolean."""
    verdicts = {
        name: _run_method(k1, p1, k2, p2, s, code) for name, code in _METHODS
    }
    unique = set(verdicts.values())
    assert len(unique) == 1, (
        f"methods disagree on (k1={k1}, k2={k2}, s={s}): {verdicts!r}"
    )
