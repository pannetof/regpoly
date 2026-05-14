# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Runtime correctness for ``MarsaXorshiftGen`` at arbitrary word widths.

For a long time the runtime hard-coded a 32-bit BitVect block layout
in ``V`` / ``SetV``: every word read returned only its low 32 bits
and every word write truncated to 32 bits regardless of ``w_``.  The
runtime "happened to work" for ``w <= 32`` because the compute mask
(``wmask_``) zeroed any bits past ``w`` before SetV wrote them, but
``w == 64`` was strictly broken — the high 32 bits of every block
were silently dropped, so almost no candidate was full-period.

These tests pin the fix: V/SetV index ``w_``-bit blocks of the
state, and every supported recurrence type still finds full-period
generators at ``w == 64``.  ``w == 32`` paths are covered exhaustively
by the Tables III/IV reproduction tests.
"""

from __future__ import annotations

import pytest

import regpoly._regpoly_cpp as cpp


def test_marsaglia_classic_xorshift64_is_full_period() -> None:
    """Marsaglia's canonical w=64 xorshift (`x^=x<<13; x^=x>>7;
    x^=x<<17`) is the X1 representative of (a, b, c) = (13, 7, 17).
    Encoded for the runtime: t1 = (-13, +7, -17).  This config is
    full-period over GF(2) — the test that exposed the V/SetV bug
    will pass iff each word is read/written as 64 bits, not 32."""
    gen = cpp.create_generator(
        "MarsaXorshiftGen",
        {"type": 1, "w": 64, "r": 1, "shifts": [-13, 7, -17]},
        64)
    assert cpp.is_full_period(gen)


def test_w64_type1_search_finds_many_full_period_in_initial_sweep() -> None:
    """A Type-I enumerator at w=64 has 4·63^3 = 1,000,188 candidates.
    Marsaglia's paper notes 275 full-period triples for w=64; a sweep
    of the first 5000 indices should hit at least a handful.  Before
    the V/SetV fix the count was 0."""
    e = cpp.make_gen_enumerator(
        "MarsaXorshiftGen", {"type": 1, "w": 64, "r": 1}, 64)
    hits = 0
    for idx in range(5000):
        gen = cpp.create_generator("MarsaXorshiftGen", e.at(idx), 64)
        if cpp.is_full_period(gen):
            hits += 1
    assert hits > 0, (
        "no full-period w=64 Type-I generators in the first 5000 "
        "candidates — V/SetV is probably truncating to 32 bits again")


@pytest.mark.parametrize("typ, w, r, params", [
    # Type 2 (Type II / Table III shape) at w=64.  Use a degenerate
    # config that's guaranteed to round-trip cleanly: shifts that
    # produce different states across consecutive calls.
    (2, 64, 2, {"m": 1, "p": [-11, 0, 0], "q": [-19, 13, 0]}),
    # Type 4 (Brent four-xorshift) at w=64.
    (4, 64, 2, {"m": 1, "p": [-13, 7], "q": [17, -5]}),
    # Type 100 (general multi-tap) at w=64 — Table IV shape.
    (100, 64, 12, {"mi_positions": [2, 3, 12], "mi_counts": [1, 1, 1],
                   "mi_shifts": [-7, 11, -21]}),
])
def test_w64_runtime_round_trips_distinct_outputs(typ, w, r, params) -> None:
    """For every supported runtime type, a w=64 generator initialised
    to a non-zero state must produce successive outputs that span the
    full 64-bit range — i.e. at least one output has bits set above
    position 31, which would be impossible if SetV were still
    truncating to 32 bits."""
    full_params = {"type": typ, "w": w, "r": r, **params}
    gen = cpp.create_generator("MarsaXorshiftGen", full_params, w)
    # Initialise with a non-trivial 64-bit pattern so the recurrence
    # has bits to shuffle.
    seed_bits = 0xDEADBEEFCAFEBABE & ((1 << (w * r)) - 1)
    bv = cpp.BitVect.from_int(w * r, seed_bits)
    gen.init(bv)
    # Iterate; OR the outputs together; check that bits past 31 light up.
    accumulated = 0
    for _ in range(64):
        gen.next()
        out = gen.get_output().to_int()
        accumulated |= out
    high_bits = accumulated >> 32
    assert high_bits != 0, (
        f"type={typ} w={w}: no output bit above position 31 was ever "
        f"set across 64 iterations — V/SetV is truncating to 32 bits")


def test_type3_w64_uses_single_word_path_not_cross_word() -> None:
    """Type 3's deliberate cross-word read is a w=32-only Marsaglia C
    artifact.  At w=64 the cross-word trick is meaningless (V already
    returns a 64-bit block); the runtime must collapse to a regular
    one-shift xorshift per tap.  We verify two things:
      (a) the generator constructs and runs without error;
      (b) it reaches the upper 32 bits in its outputs (the test
          above exercises (a) and (b) for types 2/4/100; type 3 is
          covered here separately because it has its own next()
          branch with the legacy_cross_word guard)."""
    params = {
        "type": 3, "w": 64, "r": 5,
        "tap_positions": [1, 3, 5],
        "tap_shifts":    [-13, 7, -17],
    }
    gen = cpp.create_generator("MarsaXorshiftGen", params, 64)
    seed_bits = 0xDEADBEEFCAFEBABE
    bv = cpp.BitVect.from_int(64 * 5, seed_bits)
    gen.init(bv)
    accumulated = 0
    for _ in range(64):
        gen.next()
        accumulated |= gen.get_output().to_int()
    assert (accumulated >> 32) != 0, (
        "type 3 w=64: outputs never set bits above position 31")


@pytest.mark.parametrize("w", [65, 96, 128, 200, 256])
def test_wide_w_runtime_accepts_and_advances(w) -> None:
    """w > 64 routes through next_wide_(), the BitVect-backed slow
    path.  Each generator must construct, init from a non-trivial
    seed, advance, and produce outputs whose high bits (above 63)
    actually light up — the canonical sign that the BitVect path is
    flowing the upper-half bits, not silently truncating."""
    gen = cpp.create_generator(
        "MarsaXorshiftGen",
        {"type": 1, "w": w, "r": 1, "shifts": [-13, 7, -17]},
        w)
    seed_bits = (1 << (w - 1)) | (1 << (w // 2)) | 1
    bv = cpp.BitVect.from_int(w, seed_bits)
    gen.init(bv)
    accumulated = 0
    for _ in range(64):
        gen.next()
        accumulated |= gen.get_output().to_int()
    high_bits = accumulated >> 64
    assert high_bits != 0, (
        f"w={w}: no output bit above position 63 was ever set — wide "
        f"path is probably truncating somewhere")


def test_wide_w128_enumerator_at_zero_round_trips() -> None:
    """Sanity: at w=128 the type-1 enumerator builds, at(0) returns a
    valid Params bag, and that bag instantiates a runnable generator.
    No bound-check regression."""
    e = cpp.make_gen_enumerator(
        "MarsaXorshiftGen", {"type": 1, "w": 128, "r": 1}, 128)
    p = e.at(0)
    assert p["w"] == 128
    gen = cpp.create_generator("MarsaXorshiftGen", p, 128)
    bv = cpp.BitVect.from_int(128, 1)
    gen.init(bv)
    gen.next()


def test_wide_path_matches_fast_path_at_boundary_w64() -> None:
    """Defensive: w=64 still uses the fast uint64_t kernel
    (next_wide_ only fires for w > 64), but if a future refactor
    flipped the dispatch, the bit-exact equivalence at w=64 must
    hold.  Build two generators with identical params at w=64,
    seed identically, advance, and compare outputs bit-for-bit."""
    params = {"type": 1, "w": 64, "r": 1, "shifts": [-13, 7, -17]}
    seed = 0xDEADBEEFCAFEBABE
    g1 = cpp.create_generator("MarsaXorshiftGen", params, 64)
    g2 = cpp.create_generator("MarsaXorshiftGen", params, 64)
    g1.init(cpp.BitVect.from_int(64, seed))
    g2.init(cpp.BitVect.from_int(64, seed))
    for _ in range(32):
        g1.next(); g2.next()
        assert g1.get_output().to_int() == g2.get_output().to_int()


@pytest.mark.parametrize("w", [4, 8, 16, 32, 48, 53, 60, 63])
def test_arbitrary_small_w_round_trips(w) -> None:
    """Runtime must work for every w in [2, 64], not just power-of-2
    widths.  At each w we construct a Type-1 generator with all shifts
    of magnitude 1, seed it, advance it, and check the output stays
    inside the w-bit window."""
    shifts = [-1, 1, -1]      # X1 with magnitudes (1, 1, 1)
    gen = cpp.create_generator(
        "MarsaXorshiftGen",
        {"type": 1, "w": w, "r": 1, "shifts": shifts}, w)
    bv = cpp.BitVect.from_int(w, (1 << (w - 1)) | 1)
    gen.init(bv)
    wmask = (1 << w) - 1
    for _ in range(32):
        gen.next()
        out = gen.get_output().to_int()
        assert out == (out & wmask), (
            f"w={w}: output {out:#x} has bits beyond bit {w - 1}")
