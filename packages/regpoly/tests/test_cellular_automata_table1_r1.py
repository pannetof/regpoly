# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Reproduce paper Table 1 — R1, a 32-bit single CA, NOT maximally
equidistributed.

R1's rule vector is given verbatim in Bhuvaneswari & Bhattacharjee 2026
Section 3.1:

    R1 = <90,150,90,90,90,150,150,90,90,90,90,90,150,90,90,150,
          150,90,150,150,150,90,150,150,150,150,90,150,90,150,90,150>

→ 0-indexed rule-150 positions: [1, 5, 6, 12, 15, 16, 18, 19, 20, 22,
                                  23, 24, 25, 27, 29, 31]

Paper Table 1 (k=32, L=32):

    t  | l*_t | l_t | Rank | Equidistribution
    ---+------+-----+------+------------------
    2  | 16   | 16  |  18  | not (t,l)-equi
    3  | 10   | 10  |  13  | not
    4  |  8   |  8  |  12  | not
    5  |  6   |  6  |  11  | not
    6  |  5   |  5  |  11  | not
    8  |  4   |  4  |  12  | not
    10 |  3   |  3  |  13  | not
    16 |  2   |  2  |  17  | not
    32 |  1   |  1  |  32  | EQUI

We reproduce the 9-row per-t equidistribution verdict byte-for-byte:
for each row, paper-says-equi iff our matricial DE method also says equi
at (t, l*_t).
"""

from __future__ import annotations

import pytest

from regpoly.analyses.equidistribution_test import EquidistributionTest
from regpoly.core.combination import Combination
from regpoly.core.generator import Generator


_R1_POS = [1, 5, 6, 12, 15, 16, 18, 19, 20, 22,
           23, 24, 25, 27, 29, 31]

# t -> paper's "Equidistribution?" verdict at (t, l*_t)
_PAPER_TABLE_1 = {
    2:  False,
    3:  False,
    4:  False,
    5:  False,
    6:  False,
    8:  False,
    10: False,
    16: False,
    32: True,
}


# R1 is a single CA at k=32; L = min(32, 64) = 32 (same as before, but
# noted explicitly under the new convention).  Lmax bumped to 64 to
# allow combined CAs that wrap this to honor the "L = min(L_i)" rule.
_R1_L = min(32, 64)


def _make_r1():
    gen = Generator.create("CellularAutomataGen", L=_R1_L,
                           k=32, rule150_positions=_R1_POS, s=1)
    comb = Combination.CreateFromFiles([[gen]], Lmax=64, temperings=[[]])
    next(iter(comb))
    return comb, gen


def _is_equi_at(t, ecart, k_g, L_test):
    l_star = min(L_test, k_g // t)
    if l_star == 0:
        return False
    t_star = k_g // l_star
    return t <= t_star - ecart[l_star]


def test_r1_is_primitive():
    """R1's characteristic polynomial is primitive — paper Section 3.1."""
    _, gen = _make_r1()
    assert gen._cpp_gen.is_full_period()


def test_r1_overall_not_me():
    """R1 is NOT maximally equidistributed (paper Section 3.1)."""
    comb, _ = _make_r1()
    test = EquidistributionTest(L=64, delta=[10**9] * 65,
                                mse=10**9, method=None)
    res = test.run(comb)
    assert res.is_me() is False


@pytest.mark.parametrize(
    "t,paper_equi",
    sorted(_PAPER_TABLE_1.items()),
    ids=lambda v: f"t={v}" if isinstance(v, int) else str(v),
)
def test_r1_per_t_equidistribution(t, paper_equi):
    """For each row of Table 1, the (t, l*_t)-equidistribution verdict
    must match the paper."""
    comb, _ = _make_r1()
    test = EquidistributionTest(L=64, delta=[10**9] * 65,
                                mse=10**9, method=None)
    res = test.run(comb)
    our_equi = _is_equi_at(t, res.ecart, comb.k_g, comb.L)
    assert our_equi == paper_equi, (
        f"Table 1 row t={t}: paper says equi={paper_equi}, "
        f"our test says equi={our_equi}. "
        f"ecart[l*_t={min(comb.L, comb.k_g // t)}] = "
        f"{res.ecart[min(comb.L, comb.k_g // t)]}"
    )
