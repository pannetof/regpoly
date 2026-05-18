# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Smoke tests for CellularAutomataGen.

Verifies construction, state evolution, copy semantics, and that the
characteristic polynomial recovered by Berlekamp-Massey matches the
expected degree.
"""

from __future__ import annotations

import pytest

from regpoly.core.generator import Generator
from regpoly_cpp._regpoly_cpp import BitVect


def _ca_L(k: int) -> int:
    """Default L for a single CellularAutomataGen: min(k, 64)."""
    return min(k, 64)


def _make(k, rule150_positions, s=1, L=None):
    if L is None:
        L = _ca_L(k)
    return Generator.create(
        "CellularAutomataGen", L=L,
        k=k, rule150_positions=rule150_positions, s=s,
    )


def test_construct_and_advance():
    """Build a k=31 CA with single rule-150 at cell 10; run next()."""
    gen = _make(k=31, rule150_positions=[10])
    bv = BitVect(31)
    bv.set_bit(0, 1)
    gen._cpp_gen.init(bv)

    # Initial state: cell 0 = 1, all others 0
    assert gen._cpp_gen.state().get_bit(0) == 1
    for i in range(1, 31):
        assert gen._cpp_gen.state().get_bit(i) == 0

    # After one step (rule 90 everywhere except cell 10):
    # Cell 1 = S_0 ^ S_2 = 1 ^ 0 = 1; all other cells = 0
    gen._cpp_gen.next()
    assert gen._cpp_gen.state().get_bit(0) == 0
    assert gen._cpp_gen.state().get_bit(1) == 1
    for i in range(2, 31):
        assert gen._cpp_gen.state().get_bit(i) == 0


def test_run_many_steps_state_evolves():
    """100 steps from a non-zero state must change the state at least once."""
    gen = _make(k=32, rule150_positions=[0, 14])
    bv = BitVect(32)
    bv.set_bit(15, 1)  # middle-ish bit set
    gen._cpp_gen.init(bv)

    initial_words = list(bytes(gen._cpp_gen.state().nwords() * 8))
    changed = False
    for _ in range(100):
        gen._cpp_gen.next()
        # Check at least one bit is now different from "first cell only"
        st = gen._cpp_gen.state()
        ones = sum(st.get_bit(i) for i in range(32))
        if ones != 1:
            changed = True
            break
    assert changed, "State did not evolve over 100 steps"


def test_copy_preserves_state():
    """copy() snapshots the state independently."""
    gen = _make(k=31, rule150_positions=[10])
    bv = BitVect(31)
    bv.set_bit(5, 1)
    gen._cpp_gen.init(bv)
    for _ in range(10):
        gen._cpp_gen.next()

    snapshot = gen._cpp_gen.copy()
    s_orig = [gen._cpp_gen.state().get_bit(i) for i in range(31)]
    s_copy = [snapshot.state().get_bit(i) for i in range(31)]
    assert s_orig == s_copy

    # Advance the original; snapshot should be unchanged.
    gen._cpp_gen.next()
    s_after = [gen._cpp_gen.state().get_bit(i) for i in range(31)]
    s_copy_after = [snapshot.state().get_bit(i) for i in range(31)]
    assert s_orig != s_after, "next() didn't change anything"
    assert s_copy == s_copy_after, "Snapshot was mutated by next() on the original"


def test_char_poly_degree_equals_k():
    """For primitive CAs the recovered char_poly should be degree-k."""
    gen = _make(k=31, rule150_positions=[10])
    cp = gen._cpp_gen.char_poly()
    # char_poly stores k bits (the non-leading coefficients of x^k); the
    # leading x^k term is implicit. The polynomial degree is k.
    assert cp.nbits() == 31


def test_time_spacing_s_advances_s_steps():
    """A CA with s=3 should be equivalent to three s=1 steps."""
    g1 = _make(k=31, rule150_positions=[10], s=1)
    g3 = _make(k=31, rule150_positions=[10], s=3)

    bv = BitVect(31)
    bv.set_bit(0, 1)
    g1._cpp_gen.init(bv)
    g3._cpp_gen.init(bv)

    # Advance g1 three times, g3 once
    for _ in range(3):
        g1._cpp_gen.next()
    g3._cpp_gen.next()

    s1 = [g1._cpp_gen.state().get_bit(i) for i in range(31)]
    s3 = [g3._cpp_gen.state().get_bit(i) for i in range(31)]
    assert s1 == s3, "s=3 step differs from three s=1 steps"


def test_invalid_position_raises():
    """Position outside [0, k) must raise."""
    with pytest.raises(Exception):
        _make(k=31, rule150_positions=[31])  # out of range
    with pytest.raises(Exception):
        _make(k=31, rule150_positions=[-1])
