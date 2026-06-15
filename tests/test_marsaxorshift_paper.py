# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""End-to-end runtime cross-check for ``MarsaXorshiftGen`` against the
paper-faithful Δ₁ values of Panneton & L'Ecuyer (2005) Tables III & IV.

The companion test
``packages/regpoly/tests/test_panneton_lecuyer_xorshift_tables.py``
feeds an abstract (rw × rw) GF(2) transition matrix through the C++
``dimension_equid`` kernel and never calls the runtime ``next()``.
These tests close the loop: they instantiate the generator via the
public Python wrapper (``Generator.create``), run the same
``EquidistributionTest`` pipeline the web app uses, and confirm Δ₁
matches the paper.  A regression in the C++ runtime's w-bit masking
would show up here.

This file lives in ``packages/regpoly/tests/`` (the right side of the
one-way ``regpoly_web → regpoly → regpoly_cpp`` boundary) because it
imports
``regpoly.core.generator`` / ``regpoly.core.combination`` /
``regpoly.analyses.equidistribution_test`` — packages above
``regpoly_cpp`` in the workspace's one-way dependency order.
"""

from __future__ import annotations

import pytest

from regpoly import make_combined
from regpoly.analyses.equidistribution_test import EquidistributionTest
from regpoly.core.generator import Generator

# Each entry: (type, r, m, p, q, expected Δ_1) — Table III (Type II).
_RUNTIME_TABLE_III = [
    (2, 2,  1, [-11, 0, 0], [-19, 13, 0],   4),
    (2, 4,  1, [-11, 7, 0], [-20, 0, 0],   13),
    (2, 12, 5, [11, -21, 0], [-6, 0, 0],   74),
    (2, 25, 2, [-10, 0, 0], [-19, 13, 0], 158),
]

# Each entry: (type=100, r, mi_positions, mi_shifts, expected Δ_1) — Table IV.
_RUNTIME_TABLE_IV = [
    (100, 12, [2, 3, 12], [-7, 11, -21],   96),
    (100, 25, [4, 10, 25], [-21, 11, -7], 186),
]


def _runtime_delta1(family_params: dict, L: int = 32) -> int:
    """Build a ``MarsaXorshiftGen`` from ``family_params`` and run the
    standard equidistribution test pipeline, returning Σ ecart_ℓ."""
    gen = Generator.create("MarsaXorshiftGen", L, **family_params)
    comb = make_combined(gen, Lmax=L)
    test = EquidistributionTest(L=L, delta=[10 ** 9] * (L + 1),
                                mse=10 ** 9, method=None)
    return test.run(comb).se


@pytest.mark.parametrize(
    "typ, r, m, p, q, expected", _RUNTIME_TABLE_III,
    ids=[f"III-r{row[1]}" for row in _RUNTIME_TABLE_III])
def test_runtime_matches_table_iii(typ, r, m, p, q, expected):
    params = dict(type=typ, w=32, r=r, m=m, p=p, q=q)
    assert _runtime_delta1(params) == expected


@pytest.mark.parametrize(
    "typ, r, mi_positions, mi_shifts, expected", _RUNTIME_TABLE_IV,
    ids=[f"IV-r{row[1]}" for row in _RUNTIME_TABLE_IV])
def test_runtime_matches_table_iv(typ, r, mi_positions, mi_shifts, expected):
    params = dict(type=typ, w=32, r=r,
                  mi_positions=mi_positions,
                  mi_counts=[1] * len(mi_positions),
                  mi_shifts=mi_shifts)
    assert _runtime_delta1(params) == expected
