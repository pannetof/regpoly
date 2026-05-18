# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Reproduce Tables 2 and 5 of Blackman & Vigna (2022).

Each row pins a (engine, w, r, A, B, C-or-None) tuple and a published
characteristic-polynomial weight.  We instantiate the corresponding
linear engine via ``Generator.create``, ask for its characteristic
polynomial (Berlekamp-Massey under the hood — exact for the published
parameter sets, since they are all full-period), and assert the weight
matches.

Weight here is the number of non-zero coefficients of the
characteristic polynomial including the leading ``x^k`` term (k is the
state size in bits, = r*w).  The ``char_poly`` BitVect carries only
the non-leading coefficients (bits 0..k-1 = coefficients of
``x^0..x^{k-1}``), so the comparison is
``popcount(char_poly) + 1 == published_weight``.
"""

from __future__ import annotations

import pytest

from regpoly.core.generator import Generator


def _weight(gen: Generator) -> int:
    cp = gen.char_poly()
    # cp._val is the integer encoding of the k non-leading coefficients.
    return bin(cp._val).count("1") + 1


# ── Table 2 (64-bit engines) ─────────────────────────────────────────────

@pytest.mark.parametrize(
    "family, r, A, B, C, expected_weight, label",
    [
        ("XoroshiroGen",  2, 24, 16, 37,  53, "xoroshiro128"),
        ("XoroshiroGen",  2, 49, 21, 28,  63, "xoroshiro128++"),
        ("XoshiroGen",    4, 17, 45, None, 115, "xoshiro256"),
        ("XoshiroGen",    8, 11, 21, None, 251, "xoshiro512"),
        ("XoroshiroGen", 16, 25, 27, 36, 439, "xoroshiro1024"),
    ],
    ids=lambda v: str(v) if not isinstance(v, str) else v,
)
def test_table_2_weight(family, r, A, B, C, expected_weight, label):
    """Table 2: 64-bit linear engines."""
    kwargs = {"w": 64, "r": r, "A": A, "B": B}
    if C is not None:
        kwargs["C"] = C
    gen = Generator.create(family, L=64, **kwargs)
    assert _weight(gen) == expected_weight, (
        f"{label}: expected weight {expected_weight}, got {_weight(gen)}"
    )


# ── Table 5 (32-bit engines) ─────────────────────────────────────────────

@pytest.mark.parametrize(
    "family, r, A, B, C, expected_weight, label",
    [
        ("XoroshiroGen", 2, 26, 9, 13, 31, "xoroshiro64"),
        ("XoshiroGen",   4,  9, 11, None, 55, "xoshiro128"),
    ],
    ids=lambda v: str(v) if not isinstance(v, str) else v,
)
def test_table_5_weight(family, r, A, B, C, expected_weight, label):
    """Table 5: 32-bit linear engines."""
    kwargs = {"w": 32, "r": r, "A": A, "B": B}
    if C is not None:
        kwargs["C"] = C
    gen = Generator.create(family, L=32, **kwargs)
    assert _weight(gen) == expected_weight, (
        f"{label}: expected weight {expected_weight}, got {_weight(gen)}"
    )
