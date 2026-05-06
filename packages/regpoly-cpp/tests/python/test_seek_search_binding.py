# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 2.4b: pybind11 binding for the SeekDriver.

C++ semantics are exercised in detail by test_seek_search.cpp. Here
we only check the binding shape — that callbacks fire, the iter
result struct exposes its fields, and the driver iterates the C++
Combination correctly from Python.
"""

from __future__ import annotations


def _make_tgfsr(a: int):
    from regpoly_cpp._regpoly_cpp import create_generator
    return create_generator(
        "TGFSRGen",
        {"w": 32, "r": 3, "m": 1, "a": a},
        32,
    )


def test_run_seek_search_smoke() -> None:
    from regpoly_cpp._regpoly_cpp import (
        Combination,
        SeekTestKind,
        SeekTestSpec,
        run_seek_search,
    )

    c = Combination(2, 32)
    c.component(0).add_gen(_make_tgfsr(1))
    c.component(0).add_gen(_make_tgfsr(2))
    c.component(1).add_gen(_make_tgfsr(3))
    c.component(1).add_gen(_make_tgfsr(4))
    assert c.reset() is True

    # ME-mode equidist test. With mse=INT_MAX every combo is selected.
    spec = SeekTestSpec()
    spec.kind = SeekTestKind.EquidistributionMatricial
    spec.eq_L_max_test = 32
    spec.eq_delta = [2**31 - 1] * 33
    spec.eq_mse = 2**31 - 1

    iter_count = [0]

    def on_iter(comb, iter_result):
        iter_count[0] += 1
        assert iter_result.me_ran
        assert iter_result.selected
        assert isinstance(iter_result.me_se, int)
        assert isinstance(iter_result.me_ecart, list)

    result = run_seek_search(
        c, [spec], 1, 100,
        on_prep=None, on_iter=on_iter, on_progress=None)

    assert result.nbgen == 4  # 2 * 2 cartesian
    assert result.nb_select == 4
    assert iter_count[0] == 4


def test_run_seek_search_nbtries() -> None:
    from regpoly_cpp._regpoly_cpp import (
        Combination,
        SeekTestKind,
        SeekTestSpec,
        run_seek_search,
    )

    c = Combination(1, 32)
    c.component(0).add_gen(_make_tgfsr(7))
    assert c.reset() is True

    spec = SeekTestSpec()
    spec.kind = SeekTestKind.EquidistributionMatricial
    spec.eq_L_max_test = 32
    spec.eq_delta = [2**31 - 1] * 33
    spec.eq_mse = 2**31 - 1

    prep_calls = []

    def on_prep(comb, is_retry):
        prep_calls.append(is_retry)

    result = run_seek_search(
        c, [spec], 3, 100,  # nbtries=3
        on_prep=on_prep, on_iter=None, on_progress=None)

    # 1 combo * 3 retries = 3 iterations.
    assert result.nbgen == 3
    assert prep_calls == [False, True, True]
