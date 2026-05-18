# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Port of ``packages/regpoly-cpp/tests/test_well_legacy_decode.cpp``.

Each bit position 0..31 is tested both for ``decode_d_s_mask`` and
``decode_test_mask``. Plus the error-path checks that the original C++
test verified.
"""

from __future__ import annotations

import pytest

from regpoly_legacy.well_legacy_decode import (
    _find_unique_msb_position,
    decode_d_s_mask,
    decode_test_mask,
)


@pytest.mark.parametrize("t", range(32))
def test_decode_test_mask_round_trip(t):
    """For every t in [0..31], (1 << (31 - t)) round-trips to t."""
    mask = 1 << (31 - t)
    assert decode_test_mask(mask, "ctx") == t


@pytest.mark.parametrize("s", range(32))
def test_decode_d_s_mask_round_trip(s):
    """For every s in [0..31], ~(1 << (31 - s)) (mod 2^32) round-trips to s."""
    mask = (~(1 << (31 - s))) & 0xFFFFFFFF
    assert decode_d_s_mask(mask, "ctx") == s


def test_decode_test_mask_rejects_zero():
    with pytest.raises(RuntimeError, match="expected exactly one set bit"):
        decode_test_mask(0, "ctx")


def test_decode_test_mask_rejects_multi_bit():
    # 0b11 has two set bits.
    with pytest.raises(RuntimeError, match="2 set bits"):
        decode_test_mask(0b11, "ctx")


def test_decode_d_s_mask_rejects_all_ones():
    # 0xFFFFFFFF inverts to 0 → zero set bits, decoder rejects.
    with pytest.raises(RuntimeError, match="expected exactly one zero bit"):
        decode_d_s_mask(0xFFFFFFFF, "ctx")


def test_decode_d_s_mask_rejects_multi_zero():
    # 0b...11111100 (two trailing zeros) inverts to 0b00000011 → 2 zero bits.
    mask = 0xFFFFFFFC
    with pytest.raises(RuntimeError, match="2 zero bits"):
        decode_d_s_mask(mask, "ctx")


def test_find_unique_msb_position_internal():
    """Private helper: empty / single / multi / max-bit cases."""
    assert _find_unique_msb_position(0) == -1
    assert _find_unique_msb_position(1) == 31      # bit at LSB → MSB-position 31
    assert _find_unique_msb_position(1 << 31) == 0  # bit at MSB → MSB-position 0
    assert _find_unique_msb_position(0b11) == -1    # multi-bit


def test_error_message_mentions_context_string():
    """Context label is preserved verbatim in the error message."""
    with pytest.raises(RuntimeError, match=r"my-custom-ctx"):
        decode_test_mask(0, "my-custom-ctx")
