# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Reproduce Tables III and IV from Panneton & L'Ecuyer (2005), "On the
Xorshift Random Number Generators".

The paper's recurrences (Type II for Table III, Type III for Table IV)
operate on w=32 unsigned-int words: every shift composes with a w-bit
mask before the next operation, so a left-then-right composition cannot
leak high bits back into the low w.  ``MarsaXorshiftGen.next()`` in C++
masks intermediate values to w bits after every shift to honour this
(see ``packages/regpoly-cpp/src/generators/marsaxorshift.cpp``).

This test is the paper-faithful oracle: it constructs the (rw x rw)
GF(2) transition matrix ``A`` directly from the paper's matrices (with
strict 32-bit masking on every intermediate result) and feeds the
unrolled state-trajectory into the C++ ``GaussMatrix.dimension_equid``
kernel.  Cross-checking against the runtime generator (instantiated via
the catalog YAML) is left to the smoke test in the regpoly-web suite.
"""

from __future__ import annotations

import pytest

import regpoly._regpoly_cpp as cpp

W = 32
MASK = (1 << W) - 1


# ── Single-word xorshift matrices, applied to a 32-bit integer ──────────

def _apply_shift(v: int, s: int) -> int:
    """Apply (I + H^a) where ``s = +a`` for right-shift, ``s = -a`` for left."""
    if s > 0:
        return (v ^ (v >> s)) & MASK
    return (v ^ ((v << (-s)) & MASK)) & MASK


def _apply_chain(v: int, shifts: list[int]) -> int:
    """Apply a product of xorshift matrices, written **right-to-left**.

    ``shifts = [s1, s2, s3]`` realises the matrix ``M3 M2 M1``: the
    *rightmost* matrix M1 is applied first (so its shift comes first in
    the list), then M2, then M3.  Each shift is followed by a 32-bit
    mask so cross-side compositions match the paper's unsigned-int
    semantics exactly.
    """
    for s in shifts:
        v = _apply_shift(v, s)
    return v


# ── Transition-matrix builders ─────────────────────────────────────────

def _build_A_type2(r: int, m: int, H_shifts: list[int], G_shifts: list[int]):
    """Build a (rw x rw) bit-matrix for ``v_n = G v_{n-m} + H v_{n-r}``.

    State layout (most-recent at the top of the block, indexed by the
    columns of A): ``x_n = [v_{n-r+1} | ... | v_{n-1} | v_n]``.  Each
    block is w bits, so block k spans columns ``k*w .. k*w+w-1``.
    """
    k = r * W
    # rows[i] is a Python int of width k; bit (k-1-j) is column j.
    A = [0] * k

    # Identity blocks: for block k = 0..r-2, A[k][k+1] = I.  In bit terms
    # row (k*W + b) of A has the column (k+1)*W + b set.
    for kb in range(r - 1):
        for b in range(W):
            row = kb * W + b
            col = (kb + 1) * W + b
            A[row] |= 1 << (k - 1 - col)

    # Bottom block row: A[r-1][0] = H, A[r-1][r-m] = G.
    # Bit row (r-1)*W + b is the b-th output bit of the (I + ...) applied
    # to its source block.  We construct that row by applying the chain
    # to each w-bit unit basis vector e_c (only bit c set) and reading
    # the b-th output bit.
    src_blocks = {
        0: H_shifts,             # H acts on v_{n-r+1} = block 0
        r - m: G_shifts,         # G acts on v_{n-m+1} = block r-m
    }
    for src_block, shifts in src_blocks.items():
        if shifts is None:
            continue
        for c in range(W):
            e_c = 1 << (W - 1 - c)            # column c is bit (W-1-c)
            out = _apply_chain(e_c, shifts)
            # For each output bit b, set A[(r-1)*W + b][src_block*W + c].
            for b in range(W):
                if (out >> (W - 1 - b)) & 1:
                    row = (r - 1) * W + b
                    col = src_block * W + c
                    A[row] |= 1 << (k - 1 - col)
    return A, k


def _build_A_type3(r: int, m1: int, m2: int, shift1: int, shift2: int, shift3: int):
    """Build A for ``v_n = (I+S1) v_{n-m1} + (I+S2) v_{n-m2} + (I+S3) v_{n-r}``."""
    k = r * W
    A = [0] * k

    for kb in range(r - 1):
        for b in range(W):
            row = kb * W + b
            col = (kb + 1) * W + b
            A[row] |= 1 << (k - 1 - col)

    src_blocks = [(0, [shift3]), (r - m2, [shift2]), (r - m1, [shift1])]
    for src_block, shifts in src_blocks:
        for c in range(W):
            e_c = 1 << (W - 1 - c)
            out = _apply_chain(e_c, shifts)
            for b in range(W):
                if (out >> (W - 1 - b)) & 1:
                    row = (r - 1) * W + b
                    col = src_block * W + c
                    A[row] |= 1 << (k - 1 - col)
    return A, k


# ── Apply A to a state vector and run the equidistribution kernel ──────


def _delta1_from_A(A_rows: list[int], k: int) -> tuple[int, list[int]]:
    """Build the matricial-DE matrix and run ``dimension_equid`` for every
    ℓ in [1, W].  Returns ``(Delta_1, ecart_list)``.

    Bit layout convention: each ``A_rows[i]`` is a k-bit Python int with
    bit ``k-1-j`` corresponding to column j (MSB-first).  We build the
    column-packed representation ``A_col[j]`` of the same matrix so that
    the recurrence on the output sequence ``C_{t+1} = C_t @ A`` reduces
    to an XOR over set bits.
    """
    L = W
    indice_max = k
    total_cols = indice_max * L
    out_nwords = (total_cols + 63) // 64

    # B_b (b = 0..W-1) selects column k-W+b from the state, i.e. the
    # b-th bit of the last w-bit block.  As a k-bit MSB-first int that
    # is ``1 << (W-1-b)`` (bit at position k-1-(k-W+b) = W-1-b).
    C = [1 << ((W - 1) - b) for b in range(W)]

    M_words: list[list[int]] = [[0] * out_nwords for _ in range(k)]

    for t in range(indice_max):
        # Emit: column (t*L + b) of M holds C_t[b], read bit-by-bit.
        for b in range(W):
            col = t * L + b
            wid = col // 64
            mbit = 63 - (col % 64)
            mbit_mask = 1 << mbit
            row = C[b]
            while row:
                low = row & -row
                bitnum = low.bit_length() - 1
                s = (k - 1) - bitnum
                M_words[s][wid] |= mbit_mask
                row ^= low

        # C_{t+1}[b] = C_t[b] * A (row-times-matrix), implemented as an
        # XOR over the set columns of C_t[b] of the corresponding rows of A.
        new_C = [0] * W
        for b in range(W):
            row = C[b]
            acc = 0
            while row:
                low = row & -row
                bitnum = low.bit_length() - 1
                j = (k - 1) - bitnum
                acc ^= A_rows[j]
                row ^= low
            new_C[b] = acc
        C = new_C

    mat = cpp.GaussMatrix(k, total_cols)
    for s in range(k):
        mat.set_row_from_words(s, M_words[s])

    ecart: list[int] = [0]
    se = 0
    for ell in range(1, L + 1):
        t_star = k // ell
        mat_copy = mat.copy()
        t_ell = mat_copy.dimension_equid(k, ell, L)
        gap = t_star - t_ell
        ecart.append(gap)
        se += gap
    return se, ecart


# ── Table III ───────────────────────────────────────────────────────────

# Each entry: (num, r, m, H_shifts, G_shifts, expected_Delta_1)
# Shifts use the convention: list applied right-to-left through the
# matrix product (rightmost matrix first).  Positive = right shift,
# negative = left shift.
#
# Translation of paper notation to the lists below:
#   H = (I + L^a)(I + R^b)  ↔  H_shifts = [b, -a]
#   H = (I + L^a)           ↔  H_shifts = [-a]
#   G = (I + R^c)           ↔  G_shifts = [c]
#   etc.
TABLE_III = [
    # Paper notation                                            -> shifts (rightmost factor first)
    # Num  r   m   H                        G                       H_shifts     G_shifts     Δ_1
    ( 1,  2, 1, [-11],     [-19, 13],         4),  # H=(I+L^11)            G=(I+R^13)(I+L^19)
    ( 2,  2, 1, [-11],     [13, -19],         7),  # H=(I+L^11)            G=(I+L^19)(I+R^13)
    ( 3,  2, 1, [-9, 8],   [-22],             7),  # H=(I+R^8)(I+L^9)      G=(I+L^22)
    ( 4,  2, 1, [9, -11],  [-17],             7),  # H=(I+L^11)(I+R^9)     G=(I+L^17)
    ( 5,  3, 1, [-23],     [-13, 4],         11),  # H=(I+L^23)            G=(I+R^4)(I+L^13)
    ( 6,  3, 1, [-23],     [7, -11],         12),  # H=(I+L^23)            G=(I+L^11)(I+R^7)
    ( 7,  3, 1, [-23],     [-11, 7],         12),  # H=(I+L^23)            G=(I+R^7)(I+L^11)
    ( 8,  3, 1, [-18],     [-13, 5],         13),  # H=(I+L^18)            G=(I+R^5)(I+L^13)
    ( 9,  4, 1, [-11, 7],  [-20],            13),  # H=(I+R^7)(I+L^11)     G=(I+L^20)
    (10,  4, 1, [-11, 7],  [-19],            17),  # H=(I+R^7)(I+L^11)     G=(I+L^19)
    (11,  4, 1, [-17],     [-12, 5],         17),  # H=(I+L^17)            G=(I+R^5)(I+L^12)
    (12,  4, 1, [-19],     [-7, 15],         19),  # H=(I+L^19)            G=(I+R^15)(I+L^7)
    (13,  4, 1, [7, -11],  [-19],            19),  # H=(I+L^11)(I+R^7)     G=(I+L^19)
    (14,  5, 1, [-11, 7],  [-20],            18),  # H=(I+R^7)(I+L^11)     G=(I+L^20)
    (15,  5, 1, [-11, 6],  [-20],            19),  # H=(I+R^6)(I+L^11)     G=(I+L^20)
    (16,  5, 3, [6, -9],   [-20],            25),  # H=(I+L^9)(I+R^6)      G=(I+L^20)
    (17,  5, 3, [-9, 6],   [-20],            25),  # H=(I+R^6)(I+L^9)      G=(I+L^20)
    (18,  8, 3, [-19, 13], [-8],             45),  # H=(I+R^13)(I+L^19)    G=(I+L^8)
    (19,  8, 3, [-17, 14], [-8],             48),  # H=(I+R^14)(I+L^17)    G=(I+L^8)
    (20,  8, 1, [-15, 7],  [-10],            52),  # H=(I+R^7)(I+L^15)     G=(I+L^10)
    (21,  8, 1, [-8, 11],  [-21],            54),  # H=(I+R^11)(I+L^8)     G=(I+L^21)
    (22, 12, 5, [11, -21], [-6],             74),  # H=(I+L^21)(I+R^11)    G=(I+L^6)
    (23, 12, 5, [-7, 6],   [-22],            79),  # H=(I+R^6)(I+L^7)      G=(I+L^22)
    (24, 12, 1, [-7, 8],   [-18],            84),  # H=(I+R^8)(I+L^7)      G=(I+L^18)
    (25, 12, 5, [6, -7],   [-22],            90),  # H=(I+L^7)(I+R^6)      G=(I+L^22)
    (26, 25, 9, [-11, 8],  [-18],           123),  # H=(I+R^8)(I+L^11)     G=(I+L^18)
    (27, 25, 9, [8, -11],  [-18],           137),  # H=(I+L^11)(I+R^8)     G=(I+L^18)
    (28, 25, 7, [-20],     [-5, 13],        155),  # H=(I+L^20)            G=(I+R^13)(I+L^5)
    (29, 25, 2, [-10],     [-19, 13],       158),  # H=(I+L^10)            G=(I+R^13)(I+L^19)
]


# ── Table IV ────────────────────────────────────────────────────────────

# Each entry: (num, r, m2, m1, H1^a1, H2^a2, H3^a3, expected_Delta_1)
# Recurrence: v_n = (I + H1^a1) v_{n-m1}
#                 + (I + H2^a2) v_{n-m2}
#                 + (I + H3^a3) v_{n-r}
# Shift convention as above.
TABLE_IV = [
    # Paper Table IV columns: r  m_2  m_1  H_2^{a_2}  H_1^{a_1}  H_3^{a_3}  Δ_1
    # _build_A_type3 takes:    r, m1, m2, shift1, shift2, shift3
    # where shift_k is the H_k^{a_k} term.
    # Num  r   m2  m1  shift1   shift2   shift3   Δ_1
    (31, 12, 2,  3,  11,  -7, -21,  96),  # H_2=L^7   H_1=R^{11}  H_3=L^{21}
    (32, 12, 5, 11, -18,  -5,  11, 100),  # H_2=L^5   H_1=L^{18}  H_3=R^{11}
    (33, 12, 7,  9,  11, -18,  -5, 102),  # H_2=L^{18} H_1=R^{11} H_3=L^5
    (34, 12, 5, 10, -18,  -5,  11, 103),  # H_2=L^5   H_1=L^{18}  H_3=R^{11}
    (35, 25, 4, 10,  11, -21,  -7, 186),  # H_2=L^{21} H_1=R^{11} H_3=L^7
    (36, 25, 5, 24,  11,  -5, -18, 188),  # H_2=L^5   H_1=R^{11}  H_3=L^{18}
    (37, 25, 7, 24,  11,  -5, -18, 190),  # H_2=L^5   H_1=R^{11}  H_3=L^{18}
    (38, 25, 5, 16,  11, -19,  -5, 219),  # H_2=L^{19} H_1=R^{11} H_3=L^5
]


# ── Test entry points ───────────────────────────────────────────────────

@pytest.mark.parametrize("row", TABLE_III, ids=lambda row: f"III-{row[0]}")
def test_table_iii_delta1(row):
    num, r, m, H_shifts, G_shifts, expected = row
    A, k = _build_A_type2(r, m, H_shifts, G_shifts)
    se, _ = _delta1_from_A(A, k)
    assert se == expected, f"Table III row {num}: got Delta_1={se}, expected {expected}"


@pytest.mark.parametrize("row", TABLE_IV, ids=lambda row: f"IV-{row[0]}")
def test_table_iv_delta1(row):
    num, r, m2, m1, s1, s2, s3, expected = row
    A, k = _build_A_type3(r, m1, m2, s1, s2, s3)
    se, _ = _delta1_from_A(A, k)
    assert se == expected, f"Table IV row {num}: got Delta_1={se}, expected {expected}"


# The runtime cross-check (which imports `regpoly.*`) lives in
# packages/regpoly/tests/test_marsaxorshift_paper.py so this test file
# remains within the regpoly-cpp layer's dependency closure.


if __name__ == "__main__":
    print("=== Table III ===")
    for row in TABLE_III:
        num, r, m, H_shifts, G_shifts, expected = row
        A, k = _build_A_type2(r, m, H_shifts, G_shifts)
        se, _ = _delta1_from_A(A, k)
        flag = "OK" if se == expected else "**MISMATCH**"
        print(f" {num:2d}  r={r:2d} m={m:2d}  Delta_1 = {se:4d}  (expected {expected:4d})  {flag}")
    print()
    print("=== Table IV ===")
    for row in TABLE_IV:
        num, r, m2, m1, s1, s2, s3, expected = row
        A, k = _build_A_type3(r, m1, m2, s1, s2, s3)
        se, _ = _delta1_from_A(A, k)
        flag = "OK" if se == expected else "**MISMATCH**"
        print(f" {num:2d}  r={r:2d} m1={m1:2d} m2={m2:2d}  Delta_1 = {se:4d}  "
              f"(expected {expected:4d})  {flag}")
