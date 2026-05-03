"""Phase 2.4d: pybind11 binding for the TemperingOptimizerDriver.

The C++ recursion is exercised by test_tempering_optimizer.cpp. Here
we only check the binding shape: that the config + result classes
round-trip, and that an end-to-end ME-mode run on a small fixture
returns sensible values + sync's the locator's current value back.
"""

from __future__ import annotations


def test_run_tempering_optimizer_once_smoke() -> None:
    from regpoly_cpp._regpoly_cpp import (
        TemperingOptimizerConfig,
        TemperOptCache,
        create_generator,
        create_transformation,
        run_tempering_optimizer_once,
    )

    gen = create_generator(
        "TauswortheGen",
        {
            "k": 31,
            "nb_terms": 3,
            "s": 12,
            "quicktaus": True,
            "poly": [0, 13, 31],
        },
        32,
    )
    trans = create_transformation(
        "tempMK",
        {
            "w": 32,
            "eta": 7,
            "mu": 15,
            "b": 0xE46E1700,
            "c": 0x9B868000,
        },
    )
    cache = TemperOptCache([gen], [[trans]], 31, 32)

    initial_b = 0xE46E1700
    locators = [(trans, "b", 32, initial_b)]
    # Naive safe masks: at resolution v, mask = (1 << (32 - v + 1)) - 1.
    safe_masks = [[0]]  # v=0 unused
    for v in range(1, 33):
        nbits = 32 - v + 1
        mask = (1 << nbits) - 1 if nbits > 0 else 0
        safe_masks.append([mask])

    cfg = TemperingOptimizerConfig()
    cfg.max_essais = 50
    cfg.delta = [0] * 33
    cfg.mse = 0
    cfg.n_restarts = 1
    cfg.random_seed = 42

    result, final_values = run_tempering_optimizer_once(
        cfg, cache, locators, safe_masks)

    assert result.se >= 0
    assert result.essais <= cfg.max_essais
    assert len(result.gaps) == 33
    assert len(final_values) == 1
    assert isinstance(final_values[0], int)
