"""
primitivity.py — Primitivity testing for characteristic polynomials.

Tests whether a characteristic polynomial over GF(2) is primitive
(i.e., generates a maximum-length sequence of period 2^k - 1).

For Mersenne prime exponents (where 2^k - 1 is itself prime),
irreducibility implies primitivity.  For other k, a full
factorization of 2^k - 1 is needed.  This module provides
precomputed factors from the Cunningham tables.
"""

from __future__ import annotations

import json
import os

import regpoly._regpoly_cpp as _cpp


# ═══════════════════════════════════════════════════════════════════════════
# Precomputed prime factors of Φ_k(2) for primitivity testing
# ═══════════════════════════════════════════════════════════════════════════

_FACTORS_DATA: dict | None = None
_FACTORS_FILE = os.path.join(os.path.dirname(__file__), "data",
                              "primitive_factors.json")


def _load_factors() -> dict:
    global _FACTORS_DATA
    if _FACTORS_DATA is None:
        if os.path.exists(_FACTORS_FILE):
            with open(_FACTORS_FILE) as f:
                _FACTORS_DATA = json.load(f)
        else:
            _FACTORS_DATA = {}
    return _FACTORS_DATA


def _get_factors_for_k(k: int) -> list[int] | None:
    """
    Return all prime factors of 2^k - 1, or None if the factorization
    is incomplete.

    Uses precomputed factors of Φ_d(2) for every divisor d of k.
    If any divisor has an incomplete factorization, returns None.
    """
    data = _load_factors()
    if not data:
        return None

    # For Mersenne prime exponents, 2^k - 1 is itself prime.
    # No factorization from the table needed.
    if is_mersenne_prime_exponent(k):
        return [2**k - 1]

    # Collect factors of Phi_d(2) for every divisor d of k.
    # Phi_1(2) = 1 (no prime factors), so d=1 is always safe to skip.
    factors = set()
    for d in _divisors(k):
        if d == 1:
            continue
        entry = data.get(str(d))
        if entry is None:
            return None
        if not entry["complete"]:
            return None
        for p in entry["factors"]:
            factors.add(p)
    return sorted(factors)


_MERSENNE_PRIMES = {
    2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607,
    1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213,
    19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091,
    756839, 859433, 1257787,
}


def is_mersenne_prime_exponent(k: int) -> bool:
    """Return True if 2^k - 1 is a Mersenne prime."""
    return k in _MERSENNE_PRIMES


def _divisors(n: int) -> list[int]:
    """Return all divisors of n in ascending order."""
    divs = []
    for i in range(1, int(n**0.5) + 1):
        if n % i == 0:
            divs.append(i)
            if i != n // i:
                divs.append(n // i)
    return sorted(divs)


def is_full_period(cpp_gen, k: int) -> bool:
    """
    Test whether a generator's characteristic polynomial is primitive.

    Parameters
    ----------
    cpp_gen : C++ generator object with a char_poly() method
    k : int — degree of the characteristic polynomial

    Returns True if the characteristic polynomial is primitive
    (period 2^k - 1).  Raises ValueError if the factorization
    of 2^k - 1 is unavailable.
    """
    cp = cpp_gen.char_poly()

    # Fast path: for Mersenne prime exponents, irreducible = primitive
    if is_mersenne_prime_exponent(k):
        return _cpp.is_irreducible(cp, k)

    factors = _get_factors_for_k(k)
    if factors is None:
        raise ValueError(
            f"Cannot test primitivity for k={k}: the complete "
            f"factorization of 2^{k} - 1 is not available. "
            f"Regenerate primitive_factors.json with a larger range "
            f"or provide the missing factors."
        )
    return _cpp.is_primitive_with_factors(
        cp, k, [str(p) for p in factors]
    )
