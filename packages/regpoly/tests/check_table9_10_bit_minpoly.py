"""Per-bit minimal-polynomial check for Tables 9 and 10 of Bhuvaneswari
& Bhattacharjee 2026.

Generator: (k1=31, r150=[10]) ⊕ (k2=32, r150=[0,14]) at s=7 (T9) or s=8 (T10).

Method:
  - Use the C++ packed Berlekamp-Massey (`packed_bm`) directly via the
    pybind11 binding.
  - The BitVect it returns follows the same convention as
    Generator::char_poly(): bit j = coefficient of z^j in the recurrence
    characteristic polynomial chi(z), with an implicit leading z^K term.
    This matches char_poly_to_ntl (primitivity.cpp).
  - For BM-returned linear complexity Len, the actual chi(z) has degree
    Len, with implicit leading z^Len.  The bits at positions [0, K-Len)
    of the K-bit output BitVect are zero (the polynomial does not reach
    degree K).
  - Init: every component must be non-zero (else its CA stays at 0).
    CombinedGenerator::init splits the K-bit init by prefix offsets, so
    we set the first bit of each component's slice.
"""

from __future__ import annotations

from regpoly_cpp._regpoly_cpp import (
    BitVect, CombinedGenerator, packed_bm, is_full_period
)
from regpoly.core.generator import Generator


def chi_int_from_bv(bv, K, Len):
    """Extract chi(z) of degree Len as an integer (bit i = coeff of z^i).

    The BitVect from packed_bm has K bits and follows char_poly_to_ntl's
    convention: bit j = coeff of z^j of chi (degree-K poly with implicit
    leading z^K).  For Len < K, the actual chi has degree Len and the
    first K-Len bits of bv (the high-degree end) are zero.

    Returns an integer with bit i = coefficient of z^i (and bit Len = 1
    set explicitly for the leading term).
    """
    n = 1 << Len
    for j in range(Len):
        # For Len < K, the chi(z) of degree Len lives in bv bits [0, Len).
        # The high-end bits [Len, K) of bv are the "would-be" coefficients
        # of z^Len ... z^(K-1) — all zero when Len is the true degree.
        if bv.get_bit(j):
            n |= 1 << j
    return n


def primitive_full_bv(bv, K):
    """For Len=K BM output, the BitVect is exactly the format
    is_full_period(bv, K) expects.  Return (bv, K) unchanged for
    convenience."""
    return bv, K


def poly_mul_gf2(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
    return r


def poly_div_gf2(a, b):
    """Return (q, r) with a = b*q + r in GF(2)[x]."""
    db = b.bit_length() - 1
    q = 0
    r = a
    while r.bit_length() - 1 >= db and r != 0:
        shift = (r.bit_length() - 1) - db
        q ^= 1 << shift
        r ^= b << shift
    return q, r


def poly_gcd_gf2(a, b):
    while b:
        _, r = poly_div_gf2(a, b)
        a, b = b, r
    return a


def poly_deg(p):
    return p.bit_length() - 1 if p else -1


def int_to_bv(m, deg):
    """Build a `deg`-bit BitVect such that is_full_period(bv, deg) tests
    primitivity of poly m (with bit i = coeff of z^i, leading bit deg)."""
    bv = BitVect(deg)
    for j in range(deg):
        if (m >> j) & 1:
            bv.set_bit(j, 1)
    return bv


def analyse(s):
    print(f"\n{'=' * 72}")
    print(f"Table {'9' if s == 7 else '10'} — combined (k1=31,r150=[10]) ⊕ (k2=32,r150=[0,14]) at s={s}")
    print(f"{'=' * 72}")

    g31 = Generator.create("CellularAutomataGen", L=min(31, 64),
                           k=31, rule150_positions=[10], s=s)
    g32 = Generator.create("CellularAutomataGen", L=min(32, 64),
                           k=32, rule150_positions=[0, 14], s=s)

    cg31, cg32 = g31._cpp_gen, g32._cpp_gen
    cp31_bv = cg31.char_poly()
    cp32_bv = cg32.char_poly()
    p1 = chi_int_from_bv(cp31_bv, 31, 31)
    p2 = chi_int_from_bv(cp32_bv, 32, 32)
    prim31 = cg31.is_full_period()
    prim32 = cg32.is_full_period()
    print(f"  Component 1 (k=31): primitive = {prim31}, p1(z) = 0x{p1:x}")
    print(f"  Component 2 (k=32): primitive = {prim32}, p2(z) = 0x{p2:x}")

    cg = CombinedGenerator([cg31, cg32], 64)
    K = cg.k()
    L = cg.L()
    print(f"  Combined: k = {K}, L = {L}")

    # Init that makes BOTH components non-zero.  CombinedGenerator::init
    # uses prefix_k offsets: bits [0, k1) → comp1, bits [k1, k1+k2) → comp2.
    init = BitVect(K)
    init.set_bit(0, 1)        # comp1 cell 0 = 1
    init.set_bit(31, 1)       # comp2 cell 0 = 1

    p1p2 = poly_mul_gf2(p1, p2)
    print(f"  p1·p2 (expected combined chi) = 0x{p1p2:x}, degree = {poly_deg(p1p2)}")

    # Per-bit BM
    bit_polys = []
    bit_complexities = []
    for bit_idx in range(L):
        Len_b, mp_bv = packed_bm(cg, init, K, bit_idx)
        m = chi_int_from_bv(mp_bv, K, Len_b)
        bit_polys.append(m)
        bit_complexities.append(Len_b)

    unique = set(bit_polys)
    all_same = len(unique) == 1
    print(f"\n  All {L} output bits share the same minimal polynomial: {all_same}")
    if all_same:
        common = bit_polys[0]
        common_complexity = bit_complexities[0]
        common_deg = poly_deg(common)
        print(f"  Common bit-min-poly         = 0x{common:x}")
        print(f"  Common linear complexity    = {common_complexity}")
        print(f"  Polynomial degree           = {common_deg}")
        print(f"  bit-min-poly == p1·p2 ?     = {common == p1p2}")

        # Test primitivity of the common min poly directly, on its actual degree.
        bv_for_prim = int_to_bv(common, common_deg)
        is_prim = is_full_period(bv_for_prim, common_deg)
        print(f"  bit-min-poly primitive ?    = {is_prim}")

        # If common != p1·p2, factor it.
        if common != p1p2:
            # gcd with p1·p2: the min poly should divide p1·p2 (since
            # both components contribute to the output recurrence).
            g = poly_gcd_gf2(common, p1p2)
            print(f"  gcd(common, p1·p2)          = 0x{g:x}, degree = {poly_deg(g)}")
            q1, r1 = poly_div_gf2(p1p2, common)
            print(f"  p1·p2 / common: q=0x{q1:x}, deg(q)={poly_deg(q1)}, r=0x{r1:x}")

        # Factor against p1 and p2 individually.
        gp1 = poly_gcd_gf2(common, p1)
        gp2 = poly_gcd_gf2(common, p2)
        print(f"  gcd(common, p1) = 0x{gp1:x}, deg={poly_deg(gp1)}")
        print(f"  gcd(common, p2) = 0x{gp2:x}, deg={poly_deg(gp2)}")
    else:
        print(f"  Distinct min polys ({len(unique)}):")
        for b, (p, c) in enumerate(zip(bit_polys, bit_complexities)):
            print(f"    bit {b:3d}: Len={c}, mp=0x{p:x}")

    # ── Interpretation ──
    print("\n  Interpretation:")
    if all_same:
        common = bit_polys[0]
        common_deg = poly_deg(common)
        bv_for_prim = int_to_bv(common, common_deg)
        is_prim = is_full_period(bv_for_prim, common_deg)
        print(f"    • All {L} output bits share the SAME minimal polynomial of degree {common_deg}.")
        if is_prim:
            print(f"    • That polynomial IS primitive.")
            print(f"    • Period of each output bit = 2^{common_deg} - 1.")
        else:
            print(f"    • That polynomial is NOT primitive (reducible product of two primitives).")
            print(f"    • Period of each output bit ≤ 2^{common_deg} - 1.")


if __name__ == "__main__":
    analyse(s=7)
    analyse(s=8)
