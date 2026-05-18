# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Paper Section 3.1 R2: k=1409, CA(150′), claimed primitive but NOT ME.

R2 is the largest CA in the paper (k=1409) and is referenced as
"Bhattacharjee et al. 2023 — CA emulation of Mersenne Twister".

The primitivity check is FAST because k=1409 is precomputed in
``primitive_factors.json`` (factor table covers k up to 4096).

The equidistribution check is SLOW (k=1409 ⇒ matricial DE ~ minutes
wall-clock) and is marked ``@pytest.mark.slow`` per the same precedent
as the MTToolBox cross-checks (CLAUDE.md).
"""

from __future__ import annotations

import pytest

from regpoly.analyses.equidistribution_test import EquidistributionTest
from regpoly.core.combination import Combination
from regpoly.core.generator import Generator


_R2_K = 1409
_R2_POSITIONS = list(range(1, _R2_K))  # CA(150'): rule 150 everywhere except cell 0


_R2_L = min(_R2_K, 64)


def test_r2_primitive():
    """R2 (k=1409, CA(150')) is primitive — paper Section 3.1."""
    gen = Generator.create("CellularAutomataGen", L=_R2_L,
                           k=_R2_K, rule150_positions=_R2_POSITIONS, s=1)
    assert gen._cpp_gen.is_full_period()


@pytest.mark.slow
def test_r2_not_maximally_equidistributed():
    """R2 is NOT maximally equidistributed — paper Section 3.2 Table 1
    discussion.  k=1409 makes the matricial DE pipeline minutes-long;
    excluded from the default lane via @pytest.mark.slow."""
    gen = Generator.create("CellularAutomataGen", L=_R2_L,
                           k=_R2_K, rule150_positions=_R2_POSITIONS, s=1)
    comb = Combination.CreateFromFiles([[gen]], Lmax=64, temperings=[[]])
    next(iter(comb))
    test = EquidistributionTest(L=64, delta=[10**9] * 65,
                                mse=10**9, method=None)
    res = test.run(comb)
    assert res.is_me() is False
