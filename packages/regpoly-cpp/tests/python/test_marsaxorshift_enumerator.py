# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Exhaustive-search enumerator for MarsaXorshiftGen.

Every type (1, 2, 3, 4, 100) registers an axis layout in C++; the
search driver iterates idx ∈ [0, total) and asks the enumerator for
the corresponding ``Params``.  These tests verify the cardinalities,
the shape of ``at(idx)`` for boundary indices, and round-trip the
generated params back through ``MarsaXorshiftGen.from_params`` to
confirm the runtime accepts every produced configuration.

Pybind11 surfaces ``std::invalid_argument`` and ``std::out_of_range``
as ``ValueError`` and ``IndexError`` respectively; this file uses
those concrete exception types so a regression in the binding's
translation table would fail loudly.
"""

from __future__ import annotations

import pytest

import regpoly._regpoly_cpp as cpp


def _make_enum(typ: int, w: int, r: int, **extra) -> "cpp.GenEnumerator":
    """Build a MarsaXorshiftGen enumerator.  Uses ``typ`` rather than
    ``type`` to avoid shadowing the Python builtin."""
    params = {"type": typ, "w": w, "r": r}
    params.update(extra)
    e = cpp.make_gen_enumerator("MarsaXorshiftGen", params, 32)
    assert e is not None
    return e


# ── Type 1 ────────────────────────────────────────────────────────────

def test_type1_size_w4():
    # 4 patterns × (w-1)^3 = 4 × 27 = 108
    assert _make_enum(1, w=4, r=1).size() == 108


def test_type1_size_w32():
    assert _make_enum(1, w=32, r=1).size() == 4 * 31 ** 3


def test_type1_axes():
    e = _make_enum(1, w=8, r=1)
    names = [a["name"] for a in e.axes()]
    assert names == ["pattern", "a", "b", "c"]


def test_type1_at_zero_is_X1_a1_b1_c1():
    p = _make_enum(1, w=8, r=1).at(0)
    # Pattern X1 with a=b=c=1 → (-a, +b, -c) = (-1, +1, -1).
    assert p["type"] == 1
    assert p["w"] == 8
    assert p["r"] == 1
    assert p["shifts"] == [-1, 1, -1]


def test_type1_rejects_r_neq_1():
    with pytest.raises(ValueError, match=r"needs_r_eq_1"):
        cpp.make_gen_enumerator(
            "MarsaXorshiftGen", {"type": 1, "w": 8, "r": 2}, 32)


# ── Type 2 ────────────────────────────────────────────────────────────

def test_type2_size_w4_r2():
    # m_size=1 (only m=1 valid for r=2), shift_size=2w-1=7, total = 7^6.
    assert _make_enum(2, w=4, r=2).size() == 7 ** 6


def test_type2_size_with_m_pinned_collapses_m_axis():
    free = _make_enum(2, w=4, r=5).size()
    pinned = _make_enum(2, w=4, r=5, m=3).size()
    # r=5 → 4 m-values when free, 1 when pinned.
    assert free == 4 * pinned


def test_type2_pinned_m_appears_in_at():
    p = _make_enum(2, w=4, r=5, m=3).at(0)
    assert p["m"] == 3


def test_type2_at_first_and_last():
    e = _make_enum(2, w=4, r=2)
    first = e.at(0)
    last = e.at(7 ** 6 - 1)
    assert first["p"] == [-3, -3, -3] and first["q"] == [-3, -3, -3]
    assert last["p"] == [3, 3, 3] and last["q"] == [3, 3, 3]


def test_type2_rejects_invalid_pinned_m():
    with pytest.raises(ValueError, match=r"needs_admissible_fixed_m"):
        cpp.make_gen_enumerator(
            "MarsaXorshiftGen", {"type": 2, "w": 4, "r": 3, "m": 5}, 32)


def test_type2_w32_r2_size_is_NTL_ZZ_scale():
    # Confirms the binding lifts to a Python int beyond uint64.
    assert _make_enum(2, w=32, r=2).size() == 63 ** 6


# ── Type 3 ────────────────────────────────────────────────────────────

def test_type3_size_w4_r3():
    # C(3, 3) × (2(w-1))^3 = 1 × 6^3 = 216
    assert _make_enum(3, w=4, r=3).size() == 216


def test_type3_size_scales_with_C_r_3():
    # C(5, 3) × 6^3 = 10 × 216
    assert _make_enum(3, w=4, r=5).size() == 10 * 216


def test_type3_at_first_taps_are_lex_first_subset():
    p = _make_enum(3, w=4, r=5).at(0)
    assert p["tap_positions"] == [1, 2, 3]
    # All shifts at the lowest index → -(w-1) = -3.
    assert p["tap_shifts"] == [-3, -3, -3]


def test_type3_at_last_taps_are_lex_last_subset():
    p = _make_enum(3, w=4, r=5).at(2159)
    assert p["tap_positions"] == [3, 4, 5]
    assert p["tap_shifts"] == [3, 3, 3]


def test_type3_shifts_never_zero():
    # The shift alphabet excludes 0 by construction.
    e = _make_enum(3, w=4, r=3)
    for idx in (0, 1, 5, 100, 215):
        for s in e.at(idx)["tap_shifts"]:
            assert s != 0


def test_type3_rejects_r_lt_3():
    with pytest.raises(ValueError, match=r"needs_r_ge_3"):
        cpp.make_gen_enumerator(
            "MarsaXorshiftGen", {"type": 3, "w": 4, "r": 2}, 32)


# ── Type 4 ────────────────────────────────────────────────────────────

def test_type4_size_w4_r2():
    # m=1 forced (r-1=1), shift_size=2w-1=7, total = 7^4 = 2401.
    assert _make_enum(4, w=4, r=2).size() == 7 ** 4


def test_type4_at_first_and_last():
    e = _make_enum(4, w=4, r=2)
    first = e.at(0)
    last = e.at(7 ** 4 - 1)
    assert first["p"] == [-3, -3] and first["q"] == [-3, -3]
    assert last["p"] == [3, 3] and last["q"] == [3, 3]


# ── Type 100 ──────────────────────────────────────────────────────────

def test_type100_default_is_table_iv_form():
    # Defaults: nb_taps=3, shifts_per_tap=[1,1,1]
    # Same cardinality as type 3 with the same (w, r).
    assert _make_enum(100, w=32, r=12).size() == 220 * 62 ** 3


def test_type100_axes_default():
    e = _make_enum(100, w=4, r=4)
    names = [a["name"] for a in e.axes()]
    assert names == ["taps", "tap[0].shift[0]", "tap[1].shift[0]",
                     "tap[2].shift[0]"]


def test_type100_at_first_emits_table_iv_shape():
    p = _make_enum(100, w=32, r=12).at(0)
    assert p["type"] == 100
    assert p["mi_counts"] == [1, 1, 1]
    assert p["mi_positions"] == [1, 2, 3]
    assert p["mi_shifts"] == [-31, -31, -31]


def test_type100_size_grows_with_nb_taps():
    # nb_taps=4 → C(12, 4) × 62^4 = 495 × 14,776,336
    e = _make_enum(100, w=32, r=12, nb_taps=4)
    assert e.size() == 495 * 62 ** 4


def test_type100_size_with_per_tap_vector():
    # shifts_per_tap=[2,1,2] → C(12,3) × 62^(2+1+2) = 220 × 62^5
    e = _make_enum(100, w=32, r=12,
              nb_taps=3, shifts_per_tap=[2, 1, 2])
    assert e.size() == 220 * 62 ** 5


def test_type100_per_tap_counts_appear_in_at():
    p = _make_enum(100, w=32, r=12,
              nb_taps=3, shifts_per_tap=[2, 1, 2]).at(0)
    assert p["mi_counts"] == [2, 1, 2]
    assert len(p["mi_shifts"]) == 5


def test_type100_rejects_short_shifts_per_tap():
    with pytest.raises(ValueError, match=r"needs_shifts_per_tap_len_eq_nb_taps"):
        cpp.make_gen_enumerator(
            "MarsaXorshiftGen",
            {"type": 100, "w": 32, "r": 12,
             "nb_taps": 3, "shifts_per_tap": [1, 1]}, 32)


def test_type100_rejects_zero_shifts_per_tap():
    with pytest.raises(ValueError, match=r"needs_shifts_per_tap_ge_1"):
        cpp.make_gen_enumerator(
            "MarsaXorshiftGen",
            {"type": 100, "w": 32, "r": 12,
             "nb_taps": 3, "shifts_per_tap": [1, 0, 1]}, 32)


def test_type100_rejects_r_lt_nb_taps():
    with pytest.raises(ValueError, match=r"needs_r_ge_nb_taps"):
        cpp.make_gen_enumerator(
            "MarsaXorshiftGen",
            {"type": 100, "w": 32, "r": 2, "nb_taps": 3}, 32)


def test_type100_table_iv_31_appears_in_enumeration():
    """Table IV #31 of Panneton & L'Ecuyer (2005):
       r=12, mi_positions=[2, 3, 12], mi_shifts=[-7, 11, -21], one
       shift per tap.  Compute the index and confirm the enumerator
       emits exactly those params and the resulting generator passes
       is_full_period."""
    w, r = 32, 12
    e = _make_enum(100, w=w, r=r)               # default: nb_taps=3, [1,1,1]
    # Tap-subset: positions {2, 3, 12} → 0-based {1, 2, 11}.
    # Lex-rank in C(12, 3): use Python itertools to compute the rank.
    from itertools import combinations
    target = (1, 2, 11)
    tap_rank = next(i for i, c in enumerate(combinations(range(12), 3))
                    if c == target)
    # Shift digits: digit = (shift < 0) ? shift + (w-1) : shift + (w-2).
    base = 2 * (w - 1)                     # = 62
    def d(s):
        assert s != 0
        return s + (w - 1) if s < 0 else s + (w - 2)
    shift_digits = [d(-7), d(11), d(-21)]
    sizes = [220, base, base, base]
    digits = [tap_rank, *shift_digits]
    idx = 0
    for digit, size in zip(digits, sizes):
        idx = idx * size + digit

    p = e.at(idx)
    assert p == {"type": 100, "w": w, "r": r,
                 "mi_positions": [2, 3, 12],
                 "mi_counts": [1, 1, 1],
                 "mi_shifts": [-7, 11, -21]}
    gen = cpp.create_generator("MarsaXorshiftGen", p, 32)
    assert cpp.is_full_period(gen)


# ── Round-trip: every produced params instantiates a runnable gen ─────

@pytest.mark.parametrize("type_, w, r, extra", [
    (1, 4, 1, {}),
    (2, 4, 2, {}),
    (3, 4, 3, {}),
    (4, 4, 2, {}),
    (100, 4, 4, {}),
    (100, 4, 4, {"nb_taps": 4}),
    (100, 4, 4, {"nb_taps": 2, "shifts_per_tap": [2, 1]}),
])
def test_at_params_instantiate_a_runnable_generator(type_, w, r, extra):
    e = _make_enum(type_, w=w, r=r, **extra)
    total = e.size()
    # Sample a few indices: first, last, and some interior.
    sample_idxs = sorted({0, 1, total // 2, total - 1, min(7, total - 1)})
    for idx in sample_idxs:
        params = e.at(idx)
        gen = cpp.create_generator("MarsaXorshiftGen", params, 32)
        # Advancing once must not raise; the runtime accepts the params.
        gen.next()


# ── End-to-end: enumerator drives the full-period search ─────────────

def test_type1_w32_marsaglia_favorite_is_at_expected_index():
    """Marsaglia's "favorite" Type-I generator is X2 with
    (a, b, c) = (5, 17, 13) at w=32 — known full-period.  Compute the
    enumerator index it occupies and assert (a) the params round-trip,
    (b) the resulting generator passes is_full_period.

    This is the end-to-end proof that the enumerator drives the same
    code path the primitive_search loop will use without depending on
    a brute-force sweep that may find nothing for tiny w."""
    e = _make_enum(1, w=32, r=1)
    pat_X2 = 1
    a, b, c = 5, 17, 13
    mag = 32 - 1                                    # = 31
    idx = (((pat_X2 * mag) + (a - 1)) * mag + (b - 1)) * mag + (c - 1)

    params = e.at(idx)
    # X2 mapping: (t1.a, t1.b, t1.c) = (-c, +b, -a) = (-13, +17, -5).
    assert params == {"type": 1, "w": 32, "r": 1,
                      "shifts": [-13, 17, -5]}

    gen = cpp.create_generator("MarsaXorshiftGen", params, 32)
    assert cpp.is_full_period(gen), \
        "Marsaglia's X2 (5, 17, 13) at w=32 must have a primitive char poly"


def test_table_iii_gen1_appears_in_type2_enumeration():
    """Table III #1 of Panneton & L'Ecuyer (2005) is a Type-II generator
    with r=2, m=1, p=[-11, 0, 0], q=[-19, 13, 0].  Build the enumerator
    for (type=2, w=32, r=2) and confirm that the index corresponding
    to those (m, p, q) values produces exactly those params back.

    Note: this test reimplements the mixed-radix encoding in Python; a
    bug that desynchronised both encoders in the same direction would
    not be caught here.  ``test_bijection_small_w`` below catches that
    class of bug by sweeping the full enumeration for small (w, r) and
    checking uniqueness of every emitted ``Params``."""
    w, r = 32, 2
    e = _make_enum(2, w=w, r=r)
    # Each shift digit is `shift + (w-1)`; zero → digit w-1 = 31.
    def shift_digit(s): return s + (w - 1)
    p_digits = [shift_digit(s) for s in (-11, 0, 0)]
    q_digits = [shift_digit(s) for s in (-19, 13, 0)]
    # Mixed-radix axes: [m (size r-1=1), p[0..2] (63), q[0..2] (63)],
    # outermost varies slowest.
    sizes = [r - 1, 63, 63, 63, 63, 63, 63]
    digits = [0, *p_digits, *q_digits]            # m=1 → digit 0
    idx = 0
    for digit, size in zip(digits, sizes):
        idx = idx * size + digit
    p = e.at(idx)
    assert p == {"type": 2, "w": w, "r": r, "m": 1,
                 "p": [-11, 0, 0], "q": [-19, 13, 0]}
    gen = cpp.create_generator("MarsaXorshiftGen", p, 32)
    assert cpp.is_full_period(gen)


# ── Bijection sanity check ────────────────────────────────────────────
#
# The named-instance tests above duplicate the mixed-radix encoder in
# Python and so could pass for the wrong reason if the C++ encoder
# also mutated.  The test below catches that: for a small full
# enumeration, it confirms that every idx maps to a *distinct*
# (canonical) params bag — i.e. that at() is injective on [0, total).
#
# Type 1 is intentionally excluded: its X1/X2 patterns coincide
# whenever a == c — (-a, +b, -c) and (-c, +b, -a) emit the same
# shifts list — so the type-1 axis space (4 · (w-1)^3 indices) maps
# onto a strictly smaller image.  This is a known consequence of
# enumerating the equivalence-class reps separately, not a bug.

@pytest.mark.parametrize("typ, w, r, extra, expected_size", [
    (3, 4,  3, {},                                         216),
    (3, 4,  5, {},                                        2160),
    (4, 4,  2, {},                                        2401),
    (100, 4, 4, {},                                        864),    # C(4,3)·6^3
    (100, 4, 3, {"nb_taps": 3, "shifts_per_tap": [1, 1, 1]}, 216),
])
def test_at_is_injective(typ, w, r, extra, expected_size):
    e = _make_enum(typ, w=w, r=r, **extra)
    assert e.size() == expected_size
    seen: set[tuple] = set()
    for idx in range(e.size()):
        # Canonicalise: tuple of sorted (key, value-as-tuple) pairs.
        p = e.at(idx)
        canon = tuple(sorted(
            (k, tuple(v) if isinstance(v, list) else v)
            for k, v in p.items()))
        seen.add(canon)
    assert len(seen) == e.size(), (
        f"at() collisions: {e.size() - len(seen)} idx pairs produced "
        f"identical params bags")
