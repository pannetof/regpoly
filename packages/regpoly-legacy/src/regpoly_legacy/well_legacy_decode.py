# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Python port of well_legacy_decode.cpp.

Helpers for decoding the legacy 3-mask form of WELL paper M6 back to the
paper's (q, t, s, a) signature.

The original WELLGen accepted three full 32-bit masks per M6 instance
(paramsulong[0..2]) instead of the paper's parametric (s, t) integers +
single 32-bit constant `a`. For every published WELL generator
(Panneton et al. 2006, Table II), the masks decompose losslessly:

    pu[0] = a                   (XOR constant)
    pu[1] = ~(1 << (31 - s))    (the d_s mask: zero bit at MSB-position s)
    pu[2] =  (1 << (31 - t))    (the test mask: one bit at MSB-position t)

These helpers extract s and t from the masks. They raise RuntimeError if
the mask is not in single-bit / single-zero-bit form.
"""

from __future__ import annotations

_M32 = 0xFFFFFFFF


def _find_unique_msb_position(mask32: int) -> int:
    """Return the MSB-indexed position of the *single* set bit in ``mask32``
    (so bit at position 0 is the leftmost). Returns -1 if mask32 does not
    have exactly one bit set.

    Mirrors C++ ``find_unique_msb_position`` verbatim.
    """
    mask32 &= _M32
    if mask32 == 0:
        return -1
    if bin(mask32).count("1") != 1:
        return -1
    # Lowest-set-bit index (LSB-indexed); convert to MSB-indexed.
    lsb = 0
    while ((mask32 >> lsb) & 1) == 0:
        lsb += 1
    return 31 - lsb


def decode_d_s_mask(mask: int, context: str) -> int:
    """Return s ∈ [0..31] such that ``mask == ~(1 << (31 - s)) & 0xFFFFFFFF``.

    Raises RuntimeError if no such s exists.
    """
    m = mask & _M32
    inv = (~m) & _M32
    s = _find_unique_msb_position(inv)
    if s < 0:
        nz = bin(inv).count("1")
        raise RuntimeError(
            f"{context}: cannot decode WELL paper M6 d_s mask 0x{m:x} "
            f"— expected exactly one zero bit (i.e. ~(1 << (31 - s))), "
            f"but got {nz} zero bits. The legacy 3-mask M6 form falls "
            f"outside the paper's parametric family and must be re-encoded "
            f"by hand."
        )
    return s


def decode_test_mask(mask: int, context: str) -> int:
    """Return t ∈ [0..31] such that ``mask == (1 << (31 - t))``.

    Raises RuntimeError if no such t exists.
    """
    m = mask & _M32
    t = _find_unique_msb_position(m)
    if t < 0:
        ns = bin(m).count("1")
        raise RuntimeError(
            f"{context}: cannot decode WELL paper M6 test mask 0x{m:x} "
            f"— expected exactly one set bit (i.e. (1 << (31 - t))), "
            f"but got {ns} set bits. The legacy 3-mask M6 form falls "
            f"outside the paper's parametric family and must be re-encoded "
            f"by hand."
        )
    return t
