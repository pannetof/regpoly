# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 2.4c: pybind11 binding shape test for run_tempering_search.

C++ semantics are exercised in detail by test_tempering_search.cpp.
Here we just verify the binding shape — that callbacks fire from C++,
TemperingTryOutcome round-trips correctly, and the iteration counts
line up.
"""

from __future__ import annotations


def _make_tgfsr(a: int):
    from regpoly_cpp._regpoly_cpp import create_generator
    return create_generator(
        "TGFSRGen",
        {"w": 32, "r": 3, "m": 1, "a": a},
        32,
    )


def test_run_tempering_search_smoke() -> None:
    from regpoly_cpp._regpoly_cpp import (
        Combination,
        TemperingSearchConfig,
        TemperingTryOutcome,
        run_tempering_search,
    )

    c = Combination(2, 32)
    c.component(0).add_gen(_make_tgfsr(1))
    c.component(0).add_gen(_make_tgfsr(2))
    c.component(1).add_gen(_make_tgfsr(3))
    c.component(1).add_gen(_make_tgfsr(4))
    assert c.reset() is True

    cfg = TemperingSearchConfig()
    cfg.nb_tries = 3
    cfg.progress_interval = 0

    starts = []
    tries = []
    dones = []

    def on_combo_start(comb, idx):
        starts.append(idx)

    def on_try(comb, combo_idx, try_idx, is_first):
        tries.append((combo_idx, try_idx, is_first))
        out = TemperingTryOutcome()
        out.got_result = True
        out.score = (try_idx + 1) * 10  # decreasing-better not triggered
        return out

    def on_combo_done(comb, idx, best_score, best_try):
        dones.append((idx, best_score, best_try))

    result = run_tempering_search(
        c, cfg,
        on_combo_start=on_combo_start,
        on_try=on_try,
        on_combo_done=on_combo_done,
        on_progress=None,
    )

    # 2 * 2 = 4 combos, each runs 3 tries.
    assert result.nbgen == 4
    assert result.nb_with_result == 4
    assert starts == [1, 2, 3, 4]
    assert len(tries) == 12
    # First-try flag set on tries 0, 3, 6, 9.
    for i, (_, try_idx, is_first) in enumerate(tries):
        assert is_first == (try_idx == 0)
    # Best score per combo is min over [10, 20, 30] = 10 at try_idx 0.
    assert dones == [(1, 10, 0), (2, 10, 0), (3, 10, 0), (4, 10, 0)]


def test_run_tempering_search_early_stop_on_zero() -> None:
    from regpoly_cpp._regpoly_cpp import (
        Combination,
        TemperingSearchConfig,
        TemperingTryOutcome,
        run_tempering_search,
    )

    c = Combination(1, 32)
    c.component(0).add_gen(_make_tgfsr(7))
    assert c.reset() is True

    cfg = TemperingSearchConfig()
    cfg.nb_tries = 100
    cfg.progress_interval = 0

    try_count = [0]

    def on_try(comb, combo_idx, try_idx, is_first):
        try_count[0] += 1
        out = TemperingTryOutcome()
        out.got_result = True
        out.score = 0 if try_idx == 4 else 7
        return out

    dones = []

    def on_combo_done(comb, idx, best_score, best_try):
        dones.append((best_score, best_try))

    result = run_tempering_search(
        c, cfg,
        on_combo_start=None,
        on_try=on_try,
        on_combo_done=on_combo_done,
        on_progress=None,
    )

    assert result.nbgen == 1
    assert try_count[0] == 5
    assert dones == [(0, 4)]


def test_run_tempering_search_no_combo_done_when_no_result() -> None:
    from regpoly_cpp._regpoly_cpp import (
        Combination,
        TemperingSearchConfig,
        TemperingTryOutcome,
        run_tempering_search,
    )

    c = Combination(1, 32)
    c.component(0).add_gen(_make_tgfsr(5))
    assert c.reset() is True

    cfg = TemperingSearchConfig()
    cfg.nb_tries = 4
    cfg.progress_interval = 0

    def on_try(comb, combo_idx, try_idx, is_first):
        out = TemperingTryOutcome()
        out.got_result = False
        return out

    dones = []

    def on_combo_done(comb, idx, best_score, best_try):
        dones.append(idx)

    result = run_tempering_search(
        c, cfg, on_combo_start=None, on_try=on_try,
        on_combo_done=on_combo_done, on_progress=None)

    assert result.nbgen == 1
    assert result.nb_with_result == 0
    assert dones == []
