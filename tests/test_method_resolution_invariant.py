# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""
Phase 1 — method-resolution invariant test.

The plan migrates `EquidistributionTest._resolve_method` from a
J-branched implementation (J=1 reads `C[0]._cpp_gen.default_test_method`
directly; J>=2 builds a `_cpp.CombinedF2LinearSource` and asks it) to a
single-path implementation that always asks the C++ generator handed
in. The two paths must agree for the J=1 case — `CombinedF2LinearSource`'s
override is expected to short-circuit to the inner component's
answer when J=1.

This test pins that invariant so a future regression in
`CombinedF2LinearSource::compute_default_test_method` cannot silently
change the default method.
"""

from __future__ import annotations

import pytest

import regpoly._regpoly_cpp as _cpp
from regpoly.core.generator import Generator

# (family, structural params, L) tuples covering every family that:
# - the C++ factory can construct,
# - has a meaningful `default_test_method("equidistribution")` answer.
# Chosen to keep individual constructions cheap (small k where possible).
_FIXTURES = [
    ("MTGen",            {"w": 32, "r": 3, "m": 1, "p": 0, "a": 0xC0FFEE01}, 32),
    ("WELLGen",          {"w": 32, "r": 5, "m1": 2, "m2": 3, "m3": 4,
                          "M0_type": 2, "M0_p": 0,
                          "M1_type": 0, "M1_p": 0,
                          "M2_type": 0, "M2_p": 0,
                          "M3_type": 0, "M3_p": 0}, 32),
    ("TGFSRGen",         {"w": 32, "r": 3, "m": 1, "a": 0x9908B0DF}, 32),
    ("TauswortheGen",    {"k": 31, "nb_terms": 3,
                          "poly": [0, 6, 31], "s": 18,
                          "quicktaus": True}, 32),
]


@pytest.mark.parametrize("family,params,L", _FIXTURES)
def test_default_test_method_invariant_under_J1_wrap(family, params, L):
    """For every supported family, the resolution returned by the bare
    primitive and by a J=1 CombinedF2LinearSource wrap must match."""
    try:
        gen = Generator.create(family, L=L, **params)
    except (RuntimeError, ValueError) as exc:
        pytest.skip(f"Family {family} construction skipped: {exc}")

    bare_method = gen._cpp_gen.default_test_method("equidistribution")
    combined = _cpp.CombinedF2LinearSource([gen._cpp_gen], [[]], L)
    wrapped_method = combined.default_test_method("equidistribution")

    assert bare_method == wrapped_method, (
        f"{family}: bare={bare_method!r} vs J=1 CombinedF2LinearSource wrap"
        f"={wrapped_method!r}. CombinedF2LinearSource's J=1 short-circuit"
        " must defer to the inner component's answer."
    )
