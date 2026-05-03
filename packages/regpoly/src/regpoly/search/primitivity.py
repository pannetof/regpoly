"""
primitivity.py — Thin shim over the C++ primitivity entry points.

The full algorithm (Mersenne fast path + Cunningham-style factor
table for non-Mersenne k) lives in C++ since Phase 2.4. The
Cunningham factors are embedded in the regpoly_cpp extension via
src/algebra/primitive_factors_data.cpp, which is generated from
packages/regpoly/src/regpoly/data/primitive_factors.json.
"""

from __future__ import annotations

import regpoly._regpoly_cpp as _cpp


def is_mersenne_prime_exponent(k: int) -> bool:
    """Return True if 2^k - 1 is a Mersenne prime."""
    return _cpp.is_mersenne_prime_exponent(k)


def is_full_period(cpp_gen, k: int | None = None) -> bool:
    """
    Test whether a generator's characteristic polynomial is primitive.

    Parameters
    ----------
    cpp_gen : C++ generator object
    k : int — accepted for backwards compatibility; not used (the
        generator carries its own degree)

    Returns True if the characteristic polynomial is primitive (period
    2^k - 1). Raises RuntimeError when the factorisation of 2^k - 1
    is unavailable (extend primitive_factors.json and regenerate
    primitive_factors_data.cpp).
    """
    del k  # signature kept for backwards-compatibility
    return _cpp.is_full_period(cpp_gen)
