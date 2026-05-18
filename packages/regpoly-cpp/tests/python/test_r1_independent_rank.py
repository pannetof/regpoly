# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Independent rank check for paper Table 1, row t=2.

The paper's Section 2.3 procedure is the standard L'Ecuyer convention:
build a t·l × k matrix B whose j-th column is the first l output bits
at each of t steps starting from x_0 = e_j, then read (t,l)-equi off
rank(B) = t·l (Proposition 2).

This file reimplements R1 from scratch — no `regpoly`, no C++ kernel —
and computes rank(B) for (t=2, l=16) in two ways:

  (A) Algebraically from the k×k transition matrix T:
      B = [first l rows of T;  first l rows of T^2].

  (B) By simulation: for each j = 0..k-1, seed the CA with e_j, advance
      it twice, and stack (state_1[:l], state_2[:l]) as column j.

Both methods must agree.  The paper's published rank for this row is
18, the standard-convention rank is 17 — a +1 discrepancy that has no
support in the paper's text.
"""

from __future__ import annotations

# ── R1 spec from paper §3.1 ──────────────────────────────────────────
K = 32
RULE150_POSITIONS = [1, 5, 6, 12, 15, 16, 18, 19, 20, 22, 23, 24, 25, 27, 29, 31]
R150 = set(RULE150_POSITIONS)

# ── Transition matrix built from the per-cell rules ──────────────────
# Null boundary: cells at i = -1 and i = K are pinned to 0.
# Rule 90  at cell i: S_i' = S_{i-1} ⊕ S_{i+1}.
# Rule 150 at cell i: S_i' = S_{i-1} ⊕ S_i ⊕ S_{i+1}.


def build_T() -> list[list[int]]:
    """k×k transition matrix over GF(2)."""
    T = [[0] * K for _ in range(K)]
    for i in range(K):
        if i - 1 >= 0:
            T[i][i - 1] = 1          # left neighbor contributes
        if i + 1 < K:
            T[i][i + 1] = 1          # right neighbor contributes
        if i in R150:
            T[i][i] = 1              # self only when cell i is rule 150
    return T


# ── GF(2) helpers (rows-as-lists representation) ─────────────────────

def matmul(A: list[list[int]], B: list[list[int]]) -> list[list[int]]:
    m, k = len(A), len(A[0])
    n = len(B[0])
    out = [[0] * n for _ in range(m)]
    for i in range(m):
        Ai = A[i]
        for j in range(n):
            s = 0
            for l in range(k):
                s ^= Ai[l] & B[l][j]
            out[i][j] = s
    return out


def gf2_rank(rows: list[list[int]], ncols: int) -> int:
    """Rank of (possibly non-square) binary matrix over GF(2)."""
    rows = [r[:] for r in rows]
    r = 0
    for c in range(ncols):
        pivot = None
        for i in range(r, len(rows)):
            if rows[i][c]:
                pivot = i
                break
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                rows[i] = [a ^ b for a, b in zip(rows[i], rows[r])]
        r += 1
    return r


# ── Simulation: advance one CA step on a length-k bit list ───────────

def ca_step(state: list[int]) -> list[int]:
    new = [0] * K
    for i in range(K):
        left  = state[i - 1] if i - 1 >= 0 else 0
        right = state[i + 1] if i + 1 < K else 0
        if i in R150:
            new[i] = left ^ state[i] ^ right
        else:
            new[i] = left ^ right
    return new


# ── Tests ────────────────────────────────────────────────────────────

def test_T_is_invertible_over_gf2():
    """Primitive char poly ⇒ T has full rank."""
    T = build_T()
    assert gf2_rank(T, K) == K


def test_R1_rule_count():
    """Sanity: R1 has 16 rule-150 cells out of 32."""
    assert len(RULE150_POSITIONS) == 16
    assert RULE150_POSITIONS == sorted(RULE150_POSITIONS)


def test_table1_t2_l16_rank_algebraic():
    """B = [T[:l]; T^2[:l]] at (t=2, l=16).  Standard convention.

    Matrix is 32×32; we compute rank over GF(2) by Gaussian elimination.
    """
    T = build_T()
    T2 = matmul(T, T)
    L = 16
    B = T[:L] + T2[:L]
    assert len(B) == 2 * L == 32
    assert len(B[0]) == K == 32

    rank = gf2_rank(B, K)
    assert rank == 17, (
        f"standard-convention rank at (t=2, l=16) is {rank}; "
        f"expected 17.  (Paper Table 1 publishes 18 — off by +1.)"
    )


def test_table1_t2_l16_rank_simulation():
    """Build B by simulating the CA from x_0 = e_j for j = 0..K-1.

    Column j is the first l bits of state after step 1, stacked with
    the first l bits of state after step 2.  This is the paper's
    Section 2.3 procedure read literally: t=2 steps, l=16 bits per step.
    """
    L = 16
    cols: list[list[int]] = []
    for j in range(K):
        state = [0] * K
        state[j] = 1                 # x_0 = e_j
        state = ca_step(state)       # x_1 = T · e_j   ← step 1
        bits = state[:L]
        state = ca_step(state)       # x_2 = T^2 · e_j ← step 2
        bits = bits + state[:L]
        cols.append(bits)

    # Transpose cols (k columns of 2l bits each) into rows for the rank routine.
    nrows = 2 * L
    B = [[cols[j][i] for j in range(K)] for i in range(nrows)]
    rank = gf2_rank(B, K)
    assert rank == 17, (
        f"simulation rank at (t=2, l=16) is {rank}; "
        f"expected 17 (matches algebraic computation)."
    )


def test_algebraic_and_simulation_agree():
    """Both routes must produce the same matrix — sanity for the test."""
    T = build_T()
    T2 = matmul(T, T)
    L = 16
    B_alg = T[:L] + T2[:L]

    cols = []
    for j in range(K):
        s = [0] * K
        s[j] = 1
        s = ca_step(s)
        bits = s[:L]
        s = ca_step(s)
        bits = bits + s[:L]
        cols.append(bits)
    nrows = 2 * L
    B_sim = [[cols[j][i] for j in range(K)] for i in range(nrows)]

    assert B_alg == B_sim


def test_paper_value_is_off_by_one():
    """The paper publishes rank=18 for this row.  Document the +1 gap."""
    T = build_T()
    T2 = matmul(T, T)
    L = 16
    rank = gf2_rank(T[:L] + T2[:L], K)
    paper_published_rank = 18
    assert rank == paper_published_rank - 1
