# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 0 smoke test for the regpoly_cpp pybind11 binding.

Asserts the extension loads and exposes the minimal surface we depend on. Phase 1+
will add ABI/round-trip tests for `CombinedGenerator`, the new `SearchDriver` family,
and result structs; this file only checks that the .so is present and importable.
"""

from __future__ import annotations


def test_extension_imports() -> None:
    from regpoly_cpp import _regpoly_cpp

    assert hasattr(_regpoly_cpp, "create_generator")
    assert hasattr(_regpoly_cpp, "BitVect")


def test_bitvect_basic() -> None:
    from regpoly_cpp._regpoly_cpp import BitVect

    bv = BitVect(8)
    assert bv.nbits() == 8
