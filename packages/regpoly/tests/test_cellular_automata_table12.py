# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Reproduce paper Table 12 (Bhuvaneswari & Bhattacharjee 2026).

Table 12 lists 51 combined Cattell-Zhang CA-PRNGs that the paper claims
achieve maximal equidistribution (ME).  For each, the paper reports:
- A period exponent (≈ k1+k2 when close-to-maximal, less when reduced).
- "ME" verdict.

What we reproduce:

1. **Component primitivity at the target s.** Each underlying CA at
   time-spacing s is primitive iff T is primitive AND gcd(s, ρ_i) = 1.
   The combined "close-to-maximal" period claim becomes
   ``is_full_period(comp1@s) and is_full_period(comp2@s)``.  Matches
   paper byte-for-byte.
2. **Period coprimality.** gcd(2^k1 - 1, 2^k2 - 1) = 1 is the paper's
   selection criterion; we re-verify it.  Matches paper.
3. **Maximally-equidistributed verdict.** Paper says ME for every
   Table 12 entry, computed under what we now believe to be a
   **non-standard convention** — `(t+1)·l` rows in B and equi iff
   rank ≥ t·l.  See ``test_cellular_automata_paper_disagreement.py``
   for the convention investigation.

   The C++ matricial-DE kernel implements the **standard** convention
   (t·l rows, equi iff rank = t·l) and returns NOT-ME for every
   Table 12 entry.  This test asserts the kernel's standard-convention
   verdict.  The kernel result is *correct under its own convention*
   but diverges from the paper's stated ME claim.

   See ``docs/generators/CellularAutomataGen.md`` § "Two conventions
   for the equidistribution matrix B" for the full explanation.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from regpoly.analyses.equidistribution_test import EquidistributionTest
from regpoly.core.combination import Combination
from regpoly.core.generator import Generator


def _load_cz_table():
    p = (Path(__file__).resolve().parents[3]
         / "packages" / "regpoly" / "src" / "regpoly" / "data" / "cellular_automata.json")
    with open(p) as f:
        return json.load(f)["cattell_zhang_1995"]


_CZ = _load_cz_table()


# Paper Table 12: list of (k1, k2, s, paper_period_kind)
# kind in {"max", "reduced"}; "max" = close-to-maximal (2^(k1+k2)),
# "reduced" = ρ / gcd(s, ρ).  All entries claim ME by the paper.
_TABLE_12 = [
    (31, 32, 7,  "max"),
    (31, 32, 8,  "max"),
    (31, 40, 8,  "max"),
    (35, 48, 8,  "max"),
    (41, 48, 8,  "max"),
    (43, 48, 8,  "max"),
    (47, 56, 8,  "max"),
    (31, 32, 5,  "reduced"),
    (31, 32, 6,  "reduced"),
    (31, 32, 9,  "reduced"),
    (31, 32, 10, "reduced"),
    (31, 40, 9,  "reduced"),
    (31, 40, 10, "reduced"),
    (33, 40, 9,  "reduced"),
    (33, 40, 10, "reduced"),
    (33, 56, 10, "reduced"),
    (35, 48, 7,  "reduced"),
    (35, 48, 9,  "reduced"),
    (35, 48, 10, "reduced"),
    (35, 64, 10, "reduced"),
    (39, 40, 9,  "reduced"),
    (39, 40, 10, "reduced"),
    (39, 56, 9,  "reduced"),
    (39, 56, 10, "reduced"),
    (41, 48, 10, "reduced"),
    (41, 64, 10, "reduced"),
    (43, 48, 9,  "reduced"),
    (43, 48, 10, "reduced"),
    (43, 56, 10, "reduced"),
    (45, 56, 9,  "reduced"),
    (45, 56, 10, "reduced"),
    (45, 64, 9,  "reduced"),
    (45, 64, 10, "reduced"),
    (47, 56, 9,  "reduced"),
    (47, 56, 10, "reduced"),
    (47, 64, 9,  "reduced"),
    (47, 64, 10, "reduced"),
    (47, 72, 10, "reduced"),
    (51, 56, 9,  "reduced"),
    (51, 56, 10, "reduced"),
    (53, 56, 9,  "reduced"),
    (53, 56, 10, "reduced"),
    (55, 56, 9,  "reduced"),
    (55, 56, 10, "reduced"),
    (59, 64, 10, "reduced"),
    (63, 64, 10, "reduced"),
    (63, 80, 10, "reduced"),
    (67, 72, 10, "reduced"),
    (71, 72, 10, "reduced"),
]

assert len(_TABLE_12) == 49, f"Got {len(_TABLE_12)} entries"  # +R1 + s=1 baseline = 51 (rest documented separately)


def _make_components(k1, k2, s):
    """Build the two component CAs with per-component L = min(k, 64).
    Combined L (via Combination._update_stats) is then min(L_1, L_2)."""
    p1 = _CZ[str(k1)]
    p2 = _CZ[str(k2)]
    g1 = Generator.create("CellularAutomataGen", L=min(k1, 64),
                          k=k1, rule150_positions=p1, s=s)
    g2 = Generator.create("CellularAutomataGen", L=min(k2, 64),
                          k=k2, rule150_positions=p2, s=s)
    return g1, g2


@pytest.mark.parametrize(
    "k1,k2,s,kind", _TABLE_12,
    ids=[f"k{k1}-k{k2}-s{s}" for k1, k2, s, _ in _TABLE_12],
)
def test_period_classification(k1, k2, s, kind):
    """Paper's period classification: max (gcd(s, ρ_i)=1 ∀i) vs reduced."""
    g1, g2 = _make_components(k1, k2, s)
    p1 = g1._cpp_gen.is_full_period()
    p2 = g2._cpp_gen.is_full_period()
    both_primitive_at_s = p1 and p2

    if kind == "max":
        assert both_primitive_at_s, (
            f"Paper says (k1={k1}, k2={k2}, s={s}) achieves close-to-maximal "
            f"period, but at least one component's T^s is not primitive "
            f"(p1={p1}, p2={p2})."
        )
    else:
        assert not both_primitive_at_s, (
            f"Paper says (k1={k1}, k2={k2}, s={s}) has reduced period, but "
            f"both components' T^s are primitive (p1={p1}, p2={p2})."
        )


@pytest.mark.parametrize(
    "k1,k2,s,kind", _TABLE_12,
    ids=[f"k{k1}-k{k2}-s{s}" for k1, k2, s, _ in _TABLE_12],
)
def test_component_periods_coprime(k1, k2, s, kind):
    """Paper's selection criterion: gcd(ρ_1, ρ_2) = 1."""
    rho1 = (1 << k1) - 1
    rho2 = (1 << k2) - 1
    assert math.gcd(rho1, rho2) == 1, (
        f"k1={k1}, k2={k2}: component periods are not coprime."
    )


@pytest.mark.parametrize(
    "k1,k2,s,kind", _TABLE_12,
    ids=[f"k{k1}-k{k2}-s{s}" for k1, k2, s, _ in _TABLE_12],
)
def test_kernel_reports_standard_convention_verdict(k1, k2, s, kind):
    """Paper claims ME for every Table 12 entry under the paper's
    non-standard convention.  Our C++ kernel implements the STANDARD
    convention (t·l rows, equi iff rank = t·l) and reports NOT-ME.

    This test asserts the kernel's standard-convention verdict.  The
    kernel result is correct under its own convention but diverges
    from the paper.  See test_cellular_automata_paper_disagreement.py
    for both conventions.
    """
    g1, g2 = _make_components(k1, k2, s)
    comb = Combination.CreateFromFiles([[g1], [g2]], Lmax=64,
                                       temperings=[[], []])
    next(iter(comb))
    test = EquidistributionTest(L=64, delta=[10**9] * 65,
                                mse=10**9, method=None)
    res = test.run(comb)
    assert res.is_me() is False, (
        f"({k1},{k2},s={s}): standard convention should say NOT-ME, "
        f"but kernel reports ME.  Did the kernel switch conventions?"
    )
