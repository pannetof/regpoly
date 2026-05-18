# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Locks the equidistribution behavior of CellularAutomataGen so that
any future kernel-level change is caught.

**Two conventions** for the equidistribution matrix B exist (see
``test_cellular_automata_paper_disagreement.py``):

- **Standard:** B has t·l rows from t advances; equi iff rank = t·l.
  Implemented by the C++ matricial-DE kernel (and everywhere else in
  this repo).
- **Paper:** B has (t+1)·l rows including the initial-state block;
  equi iff rank ≥ t·l.  Reproduces paper Tables 1, 3, 7-10
  byte-for-byte (modulo one t=16 typo).

The kernel uses the STANDARD convention.  This test file locks the
kernel's ecart vector for (31, 32, s=1), and asserts the ME-overall
verdicts the kernel returns for s ∈ {1, 2, 4, 7, 8}.

For s ∈ {1, 2, 4}: paper-convention and standard both say NOT ME, and
the kernel reports NOT ME.  Match.

For s ∈ {7, 8}: paper-convention (and the paper itself) say ME, but
the standard convention says NOT ME — and the kernel, implementing
the standard convention, reports NOT ME.  The kernel verdict here is
*correct under its own convention* but disagrees with the paper.
"""

from __future__ import annotations

import pytest

from regpoly.analyses.equidistribution_test import EquidistributionTest
from regpoly.core.combination import Combination
from regpoly.core.generator import Generator


# The ecart vector that the current kernel computes for combined (31,32)
# at s=1 under the matricial method (CA's default) and per-component L
# = min(k, 64).  Length is test.L + 1 = 65.  ecart[0] = -1 sentinel.
_EXPECTED_ECART_S1 = (
    [-1]                                          # ecart[0] sentinel
    + [0, 1, 20, 14, 11, 9, 8, 6, 6, 5,           # ecart[1..10]
       4, 4, 3, 3, 3,                             # ecart[11..15]
       2, 2, 2, 2, 2, 2,                          # ecart[16..21]
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1]              # ecart[22..31]
    + [0] * 33                                    # ecart[32..64]
)


def test_combined_31_32_s1_ecart_locked():
    """Lock the matricial-DE result for combined (31,32) at s=1 under
    the new convention: per-component L = min(k, 64), combined L =
    min(L_i), test.L = 64, matricial method default for combined CA."""
    g31 = Generator.create("CellularAutomataGen", L=min(31, 64), k=31,
                           rule150_positions=[10], s=1)
    g32 = Generator.create("CellularAutomataGen", L=min(32, 64), k=32,
                           rule150_positions=[0, 14], s=1)
    comb = Combination.CreateFromFiles([[g31], [g32]], Lmax=64,
                                       temperings=[[], []])
    next(iter(comb))
    test = EquidistributionTest(L=64, delta=[10**9] * 65,
                                mse=10**9, method=None)
    res = test.run(comb)
    assert res.ecart == _EXPECTED_ECART_S1, (
        f"Combined (31,32) s=1 ecart changed.\n"
        f"  expected: {_EXPECTED_ECART_S1}\n"
        f"  got:      {res.ecart}"
    )
    # The paper says NOT-ME for s=1 (Table 3) — agrees with us.
    assert res.is_me() is False


def test_combined_31_32_s1_paper_me_agreement():
    """Verify our ME-overall result agrees with paper Table 3 (s=1 → not ME)."""
    g31 = Generator.create("CellularAutomataGen", L=min(31, 64), k=31,
                           rule150_positions=[10], s=1)
    g32 = Generator.create("CellularAutomataGen", L=min(32, 64), k=32,
                           rule150_positions=[0, 14], s=1)
    comb = Combination.CreateFromFiles([[g31], [g32]], Lmax=64,
                                       temperings=[[], []])
    next(iter(comb))
    test = EquidistributionTest(L=64, delta=[10**9] * 65,
                                mse=10**9, method=None)
    res = test.run(comb)
    assert res.is_me() is False  # Paper Table 3 — agrees.


@pytest.mark.parametrize(
    "s",
    [1, 2, 4, 7, 8],
    ids=lambda v: f"s={v}",
)
def test_combined_31_32_not_me_for_all_tested_s(s):
    """For all (31,32, s ∈ {1,2,4,7,8}), the kernel and the literal-B
    ground truth agree on NOT ME.

    Paper Tables 3 and 7-8 (s=1,2,4) claim NOT ME — agrees with us.
    Paper Tables 9-10 (s=7,8) claim ME — but pure-Python GF(2) rank
    on the literal B matrix shows specific (t, l*) rows where
    rank < t·l*, so paper's ME claim is not supported by the matrix
    definition the paper itself gives.  See
    test_cellular_automata_paper_disagreement.py for ground-truth
    verification.
    """
    g31 = Generator.create("CellularAutomataGen", L=min(31, 64), k=31,
                           rule150_positions=[10], s=s)
    g32 = Generator.create("CellularAutomataGen", L=min(32, 64), k=32,
                           rule150_positions=[0, 14], s=s)
    comb = Combination.CreateFromFiles([[g31], [g32]], Lmax=64,
                                       temperings=[[], []])
    next(iter(comb))
    test = EquidistributionTest(L=64, delta=[10**9] * 65,
                                mse=10**9, method=None)
    res = test.run(comb)
    assert res.is_me() is False
