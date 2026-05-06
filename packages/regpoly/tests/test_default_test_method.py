# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Tests for Generator.default_test_method binding and the YAML-omission
path through EquidistributionTest.run().

The C++ rules are exercised exhaustively in
packages/regpoly-cpp/tests/test_default_test_method.cpp; this file
verifies that:
  - the pybind11 binding round-trips str | None,
  - EquidistributionTest accepts method=None,
  - run() resolves the default by asking the underlying C++ generator,
  - an explicit method= still wins over the generator default.
"""

import pytest

from regpoly import BitVect, Generator
from regpoly._regpoly_cpp import CombinedGenerator as _CppCombined
from regpoly.analyses.equidistribution_test import (
    METHOD_HARASE,
    METHOD_MATRICIAL,
    METHOD_NOTPRIMITIVE,
    METHOD_SIMD_NOTPRIMITIVE,
    EquidistributionTest,
)


# ── Binding round-trip ──────────────────────────────────────────────────


def test_binding_returns_str_for_known_test_type():
    gen = Generator.create("Tausworthe", 32, poly=[0, 3, 31], s=13, quicktaus=True)
    assert gen._cpp_gen.default_test_method("equidistribution") == "matricial"


def test_binding_returns_none_for_unknown_test_type():
    gen = Generator.create("Tausworthe", 32, poly=[0, 3, 31], s=13, quicktaus=True)
    assert gen._cpp_gen.default_test_method("collision_free") is None
    assert gen._cpp_gen.default_test_method("nonsense") is None


def test_binding_full_period_large_k_returns_harase():
    # Tausworthe with k=127 (Mersenne prime), primitive trinomial x^127+x+1.
    gen = Generator.create("Tausworthe", 32, poly=[0, 1, 127], s=1, quicktaus=False)
    assert gen.is_full_period
    assert gen._cpp_gen.k() > 100
    assert gen._cpp_gen.default_test_method("equidistribution") == "harase"


def test_binding_combined_returns_notprimitive():
    a = Generator.create("Tausworthe", 32, poly=[0, 3, 31], s=13, quicktaus=True)
    b = Generator.create("Tausworthe", 32, poly=[0, 1, 127], s=1, quicktaus=False)
    combined = _CppCombined([a._cpp_gen, b._cpp_gen], [[], []], 32)
    assert combined.default_test_method("equidistribution") == "notprimitive"
    assert combined.default_test_method("tuplets") is None


# ── EquidistributionTest accepts method=None ────────────────────────────


def test_equidistribution_test_accepts_method_none():
    t = EquidistributionTest(L=32, delta=[0] * 33, mse=0, method=None)
    assert t.method is None


def test_equidistribution_test_rejects_unknown_method():
    with pytest.raises(ValueError, match="Unknown method"):
        EquidistributionTest(L=32, delta=[0] * 33, mse=0, method=999)


# ── _from_params handles missing method ─────────────────────────────────


def test_from_params_missing_method_yields_none():
    t = EquidistributionTest._from_params({"max_gap_sum": 0}, Lmax=32)
    assert t.method is None


def test_from_params_explicit_method_resolves_to_constant():
    t = EquidistributionTest._from_params(
        {"max_gap_sum": 0, "method": "harase"}, Lmax=32,
    )
    assert t.method == METHOD_HARASE


def test_from_params_simd_notprimitive_resolves():
    t = EquidistributionTest._from_params(
        {"max_gap_sum": 0, "method": "simd_notprimitive"}, Lmax=32,
    )
    assert t.method == METHOD_SIMD_NOTPRIMITIVE


def test_from_params_notprimitive_resolves():
    t = EquidistributionTest._from_params(
        {"max_gap_sum": 0, "method": "notprimitive"}, Lmax=32,
    )
    assert t.method == METHOD_NOTPRIMITIVE


# ── run() resolves None via _resolve_method() ───────────────────────────


def _build_combination(*gens):
    """Minimal Combination wrapping the given Python Generator objects."""
    from regpoly.core.combination import Combination
    comb = Combination(J=len(gens), Lmax=32)
    for j, g in enumerate(gens):
        comb.components[j].add_gen(g)
    assert comb.reset()
    return comb


def test_resolve_method_for_full_period_small_returns_matricial():
    gen = Generator.create("Tausworthe", 32, poly=[0, 3, 31], s=13, quicktaus=True)
    bv = BitVect.zeros(31); bv.put_bit(0, 1)
    gen.initialize_state(bv)
    comb = _build_combination(gen)
    t = EquidistributionTest(L=32, delta=[0] * 33, mse=0, method=None)
    assert t._resolve_method(comb) == METHOD_MATRICIAL


def test_resolve_method_for_full_period_large_returns_harase():
    gen = Generator.create("Tausworthe", 32, poly=[0, 1, 127], s=1, quicktaus=False)
    bv = BitVect.zeros(127); bv.put_bit(0, 1)
    gen.initialize_state(bv)
    comb = _build_combination(gen)
    t = EquidistributionTest(L=32, delta=[0] * 33, mse=0, method=None)
    assert t._resolve_method(comb) == METHOD_HARASE


def test_resolve_method_for_combined_returns_notprimitive():
    a = Generator.create("Tausworthe", 32, poly=[0, 3, 31], s=13, quicktaus=True)
    b = Generator.create("Tausworthe", 32, poly=[0, 1, 127], s=1, quicktaus=False)
    bva = BitVect.zeros(31); bva.put_bit(0, 1)
    bvb = BitVect.zeros(127); bvb.put_bit(0, 1)
    a.initialize_state(bva)
    b.initialize_state(bvb)
    comb = _build_combination(a, b)
    t = EquidistributionTest(L=32, delta=[0] * 33, mse=0, method=None)
    assert t._resolve_method(comb) == METHOD_NOTPRIMITIVE


# ── Hex/decimal string param round-trip via dict_to_params ──────────────


def test_create_generator_parses_hex_string_param():
    """The pybind11 dict_to_params converts hex/decimal strings to int
    so YAML-loaded values like 'a': '0x9908B0DF' reach the generator
    factory as integers (not strings, which factories would read as 0).
    Regression for the mt19937 is_full_period=False bug."""
    import regpoly._regpoly_cpp as _cpp
    g = _cpp.create_generator(
        "MTGen",
        {"w": 32, "r": 624, "m": 397, "p": 31, "a": "0x9908B0DF"},
        32,
    )
    assert g.is_full_period() is True


def test_create_generator_parses_uppercase_hex():
    import regpoly._regpoly_cpp as _cpp
    g_lower = _cpp.create_generator(
        "MTGen",
        {"w": 32, "r": 624, "m": 397, "p": 31, "a": "0x9908b0df"},
        32,
    )
    g_upper = _cpp.create_generator(
        "MTGen",
        {"w": 32, "r": 624, "m": 397, "p": 31, "a": "0X9908B0DF"},
        32,
    )
    # Both must yield the same characteristic polynomial
    cp_l = g_lower.char_poly()
    cp_u = g_upper.char_poly()
    assert cp_l.nbits() == cp_u.nbits()
    for i in range(cp_l.nbits()):
        assert cp_l.get_bit(i) == cp_u.get_bit(i)


def test_create_generator_parses_decimal_string_param():
    """Decimal numeric strings should also convert to int."""
    import regpoly._regpoly_cpp as _cpp
    g_str = _cpp.create_generator(
        "MTGen",
        {"w": 32, "r": 624, "m": 397, "p": 31, "a": "2567483615"},
        32,
    )
    g_int = _cpp.create_generator(
        "MTGen",
        {"w": 32, "r": 624, "m": 397, "p": 31, "a": 2567483615},
        32,
    )
    assert g_str.is_full_period() == g_int.is_full_period()


def test_create_generator_invalid_hex_falls_back_to_string():
    """A string that begins with '0x' but isn't valid hex should not
    crash dict_to_params; it falls through to set_string. The generator
    factory then reads 'a' as 0 (default for missing int)."""
    import regpoly._regpoly_cpp as _cpp
    # 'a' = '0xZZZ' is not parseable as hex; factory will see no int.
    g = _cpp.create_generator(
        "MTGen",
        {"w": 32, "r": 624, "m": 397, "p": 31, "a": "0xZZZ"},
        32,
    )
    # Just verify construction succeeded; the resulting char_poly is
    # whatever the factory produces with a=0 (degenerate, not crashy).
    assert g.k() == 19937
