# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Reproduce paper Section 3.1's "weak" linear CA-based PRNGs.

The paper identifies four existing CA-based PRNGs and demonstrates that
none of them satisfy maximal equidistribution:

- R1: k=32 single CA with explicit rule vector
      (covered by test_cellular_automata_table1_r1.py)
- R2: k=1409, CA(150′) — Bhattacharjee 2023 Mersenne-twister emulation
      (covered by test_cellular_automata_r2_slow.py)
- R3: k=35,   CA(150′) — Jaleel et al. 2023 multi-stream PRNG
- R4: k=64,   rule 150 at the 3rd and 5th cells (0-indexed [2, 4])
      — L'Ecuyer & Simard reference

This file covers R3 and R4 (single-component, fast).  Each is asserted
primitive (paper claims so) and NOT maximally equidistributed (paper
Section 3.2).
"""

from __future__ import annotations

import pytest

from regpoly.analyses.equidistribution_test import EquidistributionTest
from regpoly.core.combination import Combination
from regpoly.core.generator import Generator


_WEAK = [
    # (label, k, rule150_positions)
    ("R3", 35, list(range(1, 35))),  # CA(150')
    ("R4", 64, [2, 4]),              # paper: "rule 150 at 3rd and 5th cells"
]


@pytest.mark.parametrize(
    "label,k,positions", _WEAK,
    ids=[label for label, _, _ in _WEAK],
)
def test_weak_is_primitive(label, k, positions):
    """Paper claims each component CA is primitive."""
    gen = Generator.create("CellularAutomataGen", L=min(k, 64),
                           k=k, rule150_positions=positions, s=1)
    assert gen._cpp_gen.is_full_period(), (
        f"{label}: claimed primitive but is_full_period == False"
    )


@pytest.mark.parametrize(
    "label,k,positions", _WEAK,
    ids=[label for label, _, _ in _WEAK],
)
def test_weak_not_maximally_equidistributed(label, k, positions):
    """Paper claims R3 and R4 are NOT maximally equidistributed."""
    gen = Generator.create("CellularAutomataGen", L=min(k, 64),
                           k=k, rule150_positions=positions, s=1)
    comb = Combination.CreateFromFiles([[gen]], Lmax=64, temperings=[[]])
    next(iter(comb))
    test = EquidistributionTest(L=64, delta=[10**9] * 65,
                                mse=10**9, method=None)
    res = test.run(comb)
    assert res.is_me() is False, (
        f"{label}: paper says NOT ME, our test says ME (ecart={res.ecart})"
    )
