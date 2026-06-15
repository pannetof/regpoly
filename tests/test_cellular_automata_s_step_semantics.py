# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Verify CellularAutomataGen's `s` time-spacing parameter via direct
simulation:

1. State equivalence: CA(s) after N steps == CA(1) after s*N steps.
2. Output decimation: bit 0 of CA(s) outputs == every s-th bit-0 sample
   of CA(1).
3. Primitivity under T^s: gcd(s, 2^k-1) = 1 iff CA(s) is primitive
   whenever CA(1) is.
"""

from __future__ import annotations

import math

import pytest

from regpoly.core.generator import Generator
from regpoly_cpp._regpoly_cpp import BitVect


def _make_ca(k, positions, s):
    return Generator.create(
        "CellularAutomataGen",
        L=min(k, 64), k=k, rule150_positions=positions, s=s,
    )


def _init_bit0(g):
    init = BitVect(g.k)
    init.set_bit(0, 1)
    g._cpp_gen.init(init)


def _state_to_int(g):
    st = g._cpp_gen.state()
    k = g.k
    n = 0
    for i in range(k):
        if st.get_bit(i):
            n |= 1 << (k - 1 - i)
    return n


_STATE_EQUIV_CASES = [
    (31, [10], 1),
    (31, [10], 2),
    (31, [10], 3),
    (31, [10], 5),
    (31, [10], 7),
    (31, [10], 8),
    (31, [10], 10),
    (32, [0, 14], 1),
    (32, [0, 14], 2),
    (32, [0, 14], 4),
    (32, [0, 14], 7),
    (32, [0, 14], 8),
    (32, [0, 14], 10),
    (65, [0], 1),
    (65, [0], 5),
    (65, [0], 10),
    (35, list(range(1, 35)), 1),
    (35, list(range(1, 35)), 7),
    (35, list(range(1, 35)), 10),
]


@pytest.mark.parametrize(("k", "positions", "s"), _STATE_EQUIV_CASES)
def test_state_equivalence_under_s(k, positions, s):
    """N steps under CA(s) == s*N steps under CA(1)."""
    N = 20
    g_s = _make_ca(k, positions, s=s)
    g_1 = _make_ca(k, positions, s=1)
    _init_bit0(g_s)
    _init_bit0(g_1)
    for _ in range(N):
        g_s._cpp_gen.next()
    for _ in range(s * N):
        g_1._cpp_gen.next()
    assert _state_to_int(g_s) == _state_to_int(g_1)


_DECIMATION_CASES = [
    (31, [10], 7),
    (32, [0, 14], 8),
    (32, [0, 14], 4),
    (31, [10], 2),
]


@pytest.mark.parametrize(("k", "positions", "s"), _DECIMATION_CASES)
def test_output_decimation_under_s(k, positions, s):
    """bit-0 sequence under CA(s) == every s-th bit-0 sample of CA(1)."""
    N = 100
    g_s = _make_ca(k, positions, s=s)
    g_1 = _make_ca(k, positions, s=1)
    _init_bit0(g_s)
    _init_bit0(g_1)
    seq_s = []
    seq_1 = []
    for _ in range(N):
        g_s._cpp_gen.next()
        seq_s.append(g_s._cpp_gen.get_output().get_bit(0))
    for _ in range(s * N):
        g_1._cpp_gen.next()
        seq_1.append(g_1._cpp_gen.get_output().get_bit(0))
    decimated = seq_1[s - 1::s]
    assert seq_s == decimated


@pytest.mark.parametrize("k,positions", [(31, [10]), (32, [0, 14])])
@pytest.mark.parametrize("s", list(range(1, 11)))
def test_primitivity_under_s_step(k, positions, s):
    """CA(s) is primitive iff CA(1) is primitive AND gcd(s, 2^k-1) = 1."""
    g1 = _make_ca(k, positions, s=1)
    prim_1 = g1._cpp_gen.is_full_period()

    g_s = _make_ca(k, positions, s=s)
    prim_s = g_s._cpp_gen.is_full_period()

    rho = (1 << k) - 1
    expected = (math.gcd(s, rho) == 1) and prim_1
    assert prim_s == expected
