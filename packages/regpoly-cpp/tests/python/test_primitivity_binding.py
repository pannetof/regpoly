"""Phase 2.4: pybind11 binding for primitivity entries.

The C++ algorithm (Mersenne fast path + Cunningham-style factor table)
is exercised in detail by the GoogleTest suite (test_primitivity.cpp).
Here we only check the binding shape: that the three entries are
exposed and behave correctly on a handful of well-known cases.
"""

from __future__ import annotations


def test_is_mersenne_prime_exponent() -> None:
    from regpoly_cpp._regpoly_cpp import is_mersenne_prime_exponent

    assert is_mersenne_prime_exponent(19937) is True
    assert is_mersenne_prime_exponent(31) is True
    assert is_mersenne_prime_exponent(20) is False
    assert is_mersenne_prime_exponent(1) is False


def test_get_primitive_factors_for_k() -> None:
    from regpoly_cpp._regpoly_cpp import get_primitive_factors_for_k

    # Aggregated factors of 2^12 - 1 = 3 * 3 * 5 * 7 * 13.
    assert get_primitive_factors_for_k(12) == ["3", "5", "7", "13"]

    # Mersenne prime exponent: single factor 2^k - 1.
    assert get_primitive_factors_for_k(31) == ["2147483647"]

    # k well beyond the embedded table -> None.
    assert get_primitive_factors_for_k(2_000_000) is None


def test_is_full_period_via_generator() -> None:
    from regpoly_cpp._regpoly_cpp import (
        create_generator,
        is_full_period,
    )

    # TGFSR with the published TT800 parameters is non-primitive
    # (k=96, full-period status depends on the chosen `a`); this just
    # ensures the binding accepts a Generator and returns a bool.
    g = create_generator(
        "TGFSRGen",
        {"w": 32, "r": 3, "m": 1, "a": 0x9908b0df},
        32,
    )
    result = is_full_period(g)
    assert isinstance(result, bool)
