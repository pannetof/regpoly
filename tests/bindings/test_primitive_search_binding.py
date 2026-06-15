# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Phase 2.4: pybind11 binding for the PrimitiveSearchDriver.

The C++ loop semantics are exercised in detail by
test_primitive_search.cpp. Here we only check the binding shape:
that the config struct accepts dicts, the callback signatures match,
and an end-to-end run produces TestedGenerator objects with the
expected attributes.
"""

from __future__ import annotations


def test_run_primitive_search_smoke() -> None:
    from regpoly_cpp._regpoly_cpp import (
        PrimitiveSearchConfig,
        run_primitive_search,
    )

    cfg = PrimitiveSearchConfig()
    cfg.family = "TGFSRGen"
    cfg.L = 32
    cfg.structural_params = {"w": 32, "r": 3}
    cfg.fixed_params = {"m": 1}
    cfg.max_tries = 100
    cfg.max_seconds = 0.0
    cfg.progress_interval = 25
    cfg.random_seed = 2024

    hits = []
    progress_events = []

    def on_hit(tg) -> None:
        hits.append(tg)

    def on_progress(sp) -> None:
        progress_events.append((sp.tries, sp.elapsed_seconds))

    total = run_primitive_search(cfg, on_hit, on_progress)
    assert total == 100

    # Every progress event should report tries that's a multiple of
    # progress_interval, except the final event which reports the
    # total.
    assert progress_events, "expected at least one progress event"
    assert progress_events[-1][0] == 100

    for tg in hits:
        assert tg.family == "TGFSRGen"
        assert tg.k > 0
        assert tg.L == 32
        assert "a" in tg.params
        assert isinstance(tg.tries_at_hit, int)
        assert 1 <= tg.tries_at_hit <= 100
