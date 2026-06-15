# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Per-bit minimal-polynomial check for paper Tables 9 and 10
(Bhuvaneswari & Bhattacharjee 2026).

Generator: (k1=31, r150=[10]) ⊕ (k2=32, r150=[0,14]) at s=7 (Table 9)
or s=8 (Table 10).

Asserts via the C++ packed Berlekamp-Massey binding that every output
bit shares the same minimal polynomial of degree 63 = 31 + 32, and
that this minimal polynomial equals ``chi_1 · chi_2`` (the product of
the two component characteristic polynomials).

The paper claims these combined PRNGs are maximally equidistributed,
which would require ``chi_1 · chi_2`` to be primitive — but the
``is_full_period`` kernel says otherwise. The kernel verdict is the
regression target; the paper/kernel disagreement is discussed in
``docs/generators/CellularAutomataGen.md``.
"""

from __future__ import annotations

import pytest

from regpoly.core.generator import Generator
from regpoly_cpp._regpoly_cpp import (
    BitVect,
    CombinedF2LinearSource,
    is_full_period,
    packed_bm,
)


def _chi_int_from_bv(bv, K, Len):
    n = 1 << Len
    for j in range(Len):
        if bv.get_bit(j):
            n |= 1 << j
    return n


def _poly_mul_gf2(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
    return r


def _poly_deg(p):
    return p.bit_length() - 1 if p else -1


def _int_to_bv(m, deg):
    bv = BitVect(deg)
    for j in range(deg):
        if (m >> j) & 1:
            bv.set_bit(j, 1)
    return bv


def _build_combined(s):
    g31 = Generator.create("CellularAutomataGen", L=min(31, 64),
                           k=31, rule150_positions=[10], s=s)
    g32 = Generator.create("CellularAutomataGen", L=min(32, 64),
                           k=32, rule150_positions=[0, 14], s=s)
    cg = CombinedF2LinearSource([g31._cpp_gen, g32._cpp_gen], 64)
    init = BitVect(cg.k())
    init.set_bit(0, 1)
    init.set_bit(31, 1)
    return g31, g32, cg, init


@pytest.mark.parametrize("s", [7, 8], ids=["Table9_s7", "Table10_s8"])
def test_bit_minpoly_is_product_of_component_charpolys(s):
    """All 64 output bits share min poly = chi_1 · chi_2 of degree 63."""
    g31, g32, cg, init = _build_combined(s)
    K = cg.k()
    L = cg.L()
    assert K == 63

    assert g31._cpp_gen.is_full_period()
    assert g32._cpp_gen.is_full_period()

    p1 = _chi_int_from_bv(g31._cpp_gen.char_poly(), 31, 31)
    p2 = _chi_int_from_bv(g32._cpp_gen.char_poly(), 32, 32)
    p1p2 = _poly_mul_gf2(p1, p2)
    assert _poly_deg(p1p2) == 63

    bit_polys = []
    for bit_idx in range(L):
        Len_b, mp_bv = packed_bm(cg, init, K, bit_idx)
        bit_polys.append((Len_b, _chi_int_from_bv(mp_bv, K, Len_b)))

    lens = {ln for ln, _ in bit_polys}
    polys = {m for _, m in bit_polys}
    assert lens == {63}, f"bits have differing linear complexities: {lens}"
    assert polys == {p1p2}, (
        f"bits do not share min poly = chi_1*chi_2 (s={s})"
    )


@pytest.mark.parametrize("s", [7, 8], ids=["Table9_s7", "Table10_s8"])
def test_chi1_chi2_is_not_primitive_at_s_7_8(s):
    """chi_1 · chi_2 is reducible at s=7 and s=8 — kernel says NOT
    primitive. (Paper claims ME, which would require primitivity; the
    paper-vs-kernel disagreement is the documented investigation in
    ``docs/generators/CellularAutomataGen.md``.)"""
    g31, g32, cg, _init = _build_combined(s)
    p1 = _chi_int_from_bv(g31._cpp_gen.char_poly(), 31, 31)
    p2 = _chi_int_from_bv(g32._cpp_gen.char_poly(), 32, 32)
    p1p2 = _poly_mul_gf2(p1, p2)
    deg = _poly_deg(p1p2)
    assert deg == 63

    bv = _int_to_bv(p1p2, deg)
    assert is_full_period(bv, deg) is False, (
        f"chi_1*chi_2 unexpectedly primitive at s={s}; kernel verdict changed"
    )
