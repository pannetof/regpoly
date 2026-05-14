# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Integration tests for ``regpoly.core.parametric.build_gen_enumerator``
against ``MarsaXorshiftGen``.

The C++ enumerator's behaviour is covered by
``packages/regpoly-cpp/tests/python/test_marsaxorshift_enumerator.py``.
Here we cover the Python wrapper specifically — the ``NotEnumerable``
exception-massaging in
``regpoly.core.parametric.build_gen_enumerator`` strips the pybind11
class-name prefix from the C++ message and re-raises with a
``needs_*`` reason tag.  If the pybind translation table or the
prefix-stripping logic regresses, the primitive-search web route
would silently route ``RuntimeError`` straight to a 500 response
instead of the intended 422 ``needs_*`` error code.

These tests run against the runtime built from the marsaxorshift
ctor invariants and the per-type enumerator subclasses landed in the
same patch series, so they double as a regression suite for the
exception-name normalization.
"""

from __future__ import annotations

import pytest

from regpoly.core.parametric import Enumerator, NotEnumerable, build_gen_enumerator


def test_returns_enumerator_for_complete_type1_inputs():
    e = build_gen_enumerator("MarsaXorshiftGen", L=32,
                             resolved={"type": 1, "w": 32, "r": 1})
    assert isinstance(e, Enumerator)
    assert e.total == 4 * 31 ** 3
    assert [a["name"] for a in e.axes] == ["pattern", "a", "b", "c"]


def test_returns_enumerator_for_type100_with_default_pins():
    e = build_gen_enumerator("MarsaXorshiftGen", L=32,
                             resolved={"type": 100, "w": 32, "r": 12})
    assert isinstance(e, Enumerator)
    # Default nb_taps=3, shifts_per_tap=[1,1,1] → Table IV cardinality.
    assert e.total == 220 * 62 ** 3


@pytest.mark.parametrize("resolved, expected_reason", [
    # type 1 must have r=1 (or omit r)
    ({"type": 1, "w": 32, "r": 2},          "needs_r_eq_1_for_type1"),
    # type 3 needs r >= 3
    ({"type": 3, "w": 32, "r": 2},          "needs_r_ge_3_for_type3"),
    # type 100 with nb_taps > r
    ({"type": 100, "w": 32, "r": 2, "nb_taps": 3},
                                            "needs_r_ge_nb_taps"),
    # type 100 with mismatched shifts_per_tap length
    ({"type": 100, "w": 32, "r": 12, "nb_taps": 3,
      "shifts_per_tap": [1, 1]},            "needs_shifts_per_tap_len_eq_nb_taps"),
    # type 2 with admissible-fixed-m violation
    ({"type": 2, "w": 32, "r": 3, "m": 5},  "needs_admissible_fixed_m"),
    # missing required structural input
    ({"type": 2, "w": 32},                   "needs_r"),
])
def test_NotEnumerable_with_needs_reason(resolved, expected_reason):
    with pytest.raises(NotEnumerable) as excinfo:
        build_gen_enumerator("MarsaXorshiftGen", L=32, resolved=resolved)
    assert excinfo.value.reason == expected_reason


def test_unsupported_family_returns_None():
    # WELLGen has no enumerator registered.
    e = build_gen_enumerator("WELLGen", L=32, resolved={})
    assert e is None
