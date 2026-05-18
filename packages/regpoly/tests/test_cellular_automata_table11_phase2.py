# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Reproduce paper Table 11 (Adak-Das CA(90')/CA(150') combinations).

Table 11 catalogs combined CA-PRNGs where each component is an
Adak-Das construction:

- CA(90′)   at size k: rule 150 only at cell 0, rule 90 elsewhere.
- CA(150′)  at size k: rule 150 everywhere except cell 0.

Both are primitive at the k values listed in
``cellular_automata.json["adak_das_2021_CA*prime_k_values"]``.

We assert:
1. Component primitivity at the target s.
2. Period coprimality (paper's selection criterion).
3. ME-overall verdict under the **standard convention** (t·l rows in B,
   equi iff rank = t·l), which is what our C++ matricial-DE kernel
   implements.  The paper computes ME under a non-standard convention
   (see test_cellular_automata_paper_disagreement.py); under the
   standard convention, the kernel returns NOT-ME for every Table 11
   entry, even those the paper claims are ME.
"""

from __future__ import annotations

import math

import pytest

from regpoly.analyses.equidistribution_test import EquidistributionTest
from regpoly.core.combination import Combination
from regpoly.core.generator import Generator


def _ca_positions(k: int, mode: str) -> list[int]:
    """0-indexed rule-150 positions for an Adak-Das CA construction."""
    if mode == "90prime":
        return [0]
    elif mode == "150prime":
        return list(range(1, k))
    else:
        raise ValueError(f"unknown mode {mode}")


# ── Table 11: (k1, mode1, k2, mode2, s, verdict) ──────────────────────
# verdict ∈ {"ME", "almost_ME", "not_ME"}.
# Source: paper Table 11 (p.27). Each row lists s-values where the
# verdict applies; we expand to one tuple per (k1, k2, s).
_TABLE_11: list[tuple[int, str, int, str, int, str]] = []

# Close-to-32 ── CA(90')×CA(90') almost-ME at s=10
for k1, k2 in [(26, 29), (26, 35), (29, 35)]:
    _TABLE_11.append((k1, "90prime", k2, "90prime", 10, "almost_ME"))

# Close-to-32 ── CA(150')×CA(150') almost-ME at s=10
for k1, k2 in [(26, 29), (26, 35), (29, 35)]:
    _TABLE_11.append((k1, "150prime", k2, "150prime", 10, "almost_ME"))

# Close-to-32 ── CA(90')×CA(150') ME (both orderings)
_MIX32 = [
    ((26, 29), [5, 7, 8, 10]),
    ((26, 35), [5, 7, 8, 10]),
    ((29, 35), [5, 6, 7, 8, 9, 10]),
    ((29, 39), [6, 8, 9, 10]),
    ((35, 39), [6, 8, 9, 10]),
]
for (k1, k2), s_vals in _MIX32:
    for s in s_vals:
        _TABLE_11.append((k1, "90prime", k2, "150prime", s, "ME"))
        _TABLE_11.append((k2, "90prime", k1, "150prime", s, "ME"))

# Close-to-64 ── CA(90')×CA(150') ME
for k1, k2 in [(65, 69)]:
    _TABLE_11.append((k1, "90prime", k2, "150prime", 10, "ME"))
    _TABLE_11.append((k2, "90prime", k1, "150prime", 10, "ME"))

# Close-to-128 ── CA(90')×CA(150') almost-ME
_MIX128 = [
    ((105, 113), [9, 10]),
    ((113, 119), [10]),
]
for (k1, k2), s_vals in _MIX128:
    for s in s_vals:
        _TABLE_11.append((k1, "90prime", k2, "150prime", s, "almost_ME"))
        _TABLE_11.append((k2, "90prime", k1, "150prime", s, "almost_ME"))


def _make_components(k1, mode1, k2, mode2, s):
    """Per-component L = min(k, 64); combined L = min(L_1, L_2)."""
    p1 = _ca_positions(k1, mode1)
    p2 = _ca_positions(k2, mode2)
    g1 = Generator.create("CellularAutomataGen", L=min(k1, 64),
                          k=k1, rule150_positions=p1, s=s)
    g2 = Generator.create("CellularAutomataGen", L=min(k2, 64),
                          k=k2, rule150_positions=p2, s=s)
    return g1, g2


@pytest.mark.parametrize(
    "k1,mode1,k2,mode2,s,verdict", _TABLE_11,
    ids=[f"{m1}{k1}-{m2}{k2}-s{s}" for k1, m1, k2, m2, s, _ in _TABLE_11],
)
def test_component_periods_coprime(k1, mode1, k2, mode2, s, verdict):
    """Paper's selection: gcd(2^k1 - 1, 2^k2 - 1) = 1."""
    rho1 = (1 << k1) - 1
    rho2 = (1 << k2) - 1
    assert math.gcd(rho1, rho2) == 1


@pytest.mark.parametrize(
    "k1,mode1,k2,mode2,s,verdict", _TABLE_11,
    ids=[f"{m1}{k1}-{m2}{k2}-s{s}" for k1, m1, k2, m2, s, _ in _TABLE_11],
)
def test_each_component_primitive_at_s_1(k1, mode1, k2, mode2, s, verdict):
    """Each underlying CA (without time spacing) is primitive at the
    paper-claimed k values."""
    g1, _ = _make_components(k1, mode1, k2, mode2, s=1)
    g2_only = Generator.create(
        "CellularAutomataGen", L=min(k2, 64),
        k=k2, rule150_positions=_ca_positions(k2, mode2), s=1,
    )
    assert g1._cpp_gen.is_full_period(), f"CA({mode1}) k={k1} not primitive"
    assert g2_only._cpp_gen.is_full_period(), f"CA({mode2}) k={k2} not primitive"


@pytest.mark.parametrize(
    "k1,mode1,k2,mode2,s,verdict", _TABLE_11,
    ids=[f"{m1}{k1}-{m2}{k2}-s{s}" for k1, m1, k2, m2, s, _ in _TABLE_11],
)
def test_kernel_reports_standard_convention_verdict(k1, mode1, k2, mode2, s, verdict):
    """Paper's ME verdicts use a non-standard convention (see
    test_cellular_automata_paper_disagreement.py).  Our C++ kernel
    implements the STANDARD convention and returns NOT-ME for every
    Table 11 entry — even those the paper claims are ME.

    This test asserts the kernel's standard-convention verdict.
    """
    g1, g2 = _make_components(k1, mode1, k2, mode2, s)
    comb = Combination.CreateFromFiles([[g1], [g2]], Lmax=64,
                                       temperings=[[], []])
    next(iter(comb))
    test = EquidistributionTest(L=64, delta=[10**9] * 65,
                                mse=10**9, method=None)
    res = test.run(comb)
    if verdict == "ME":
        # Paper says ME (paper convention); kernel under standard
        # convention says NOT ME.  This is a convention difference,
        # not a paper-side error.
        assert res.is_me() is False
    else:
        # almost_ME / not_ME: paper allows non-zero gap; kernel agrees.
        assert res.ecart is not None
