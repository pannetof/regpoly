// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "temper_optimizer.h"
#include "transformation.h"

#include <cstdint>
#include <string>
#include <vector>

/**
 * @file tempering_optimizer.h
 * @brief Bitmask-search driver for the per-component tempering optimiser.
 * @ingroup core
 *
 * Phase 2.4d: tempering bitmask optimisation driver in C++. Replaces
 * the recursive `optimize(v)` + `run_minimize()` body of
 * `packages/regpoly/src/regpoly/search/tempering_optimizer.py` with a
 * pure-C++ implementation that operates on:
 *
 *   - a `TemperOptCache` (already-existing C++ stateful kernel that
 *     knows how to compute one resolution's gap incrementally), and
 *   - a list of `TemperParamLocator` entries — one per optimisable
 *     bitmask parameter (currently `b` for tempMK and `b` for
 *     `laggedTempering`).
 *
 * Safe-mask construction stays in Python (it is a small, structural
 * computation that depends on `mu` / width and changes infrequently).
 * The driver receives the precomputed `safe_masks` and just consumes
 * them.
 *
 * Three operating modes are selected by the caller through the
 * `(delta, mse, n_restarts)` triple in `TemperingOptimizerConfig`:
 *
 *  - **ME mode** — `delta` all-zero and `mse == 0`: find a bitmask
 *    with zero equidistribution defect.
 *  - **Targeted mode** — caller-supplied `delta` / `mse`: find a
 *    bitmask satisfying both caps.
 *  - **Minimisation mode** — `n_restarts > 1`: iteratively tighten
 *    the `delta` budget after every successful run.
 *
 * @see :py:class:`regpoly.search.tempering_optimizer.TemperingOptimizer`
 */

namespace regpoly::core {

/**
 * @brief Locator describing one optimisable bitmask parameter.
 *
 * Bound to a non-owning `Transformation*` plus the field's name,
 * width, and current value. The driver mutates `current_value`
 * in place and writes the best-found value back into `trans` via
 * `Transformation::update` at the end of each accepted iteration.
 *
 * @ingroup core
 */
struct TemperParamLocator {
    Transformation* trans;       ///< Target transformation (not owned).
    std::string param_name;      ///< Parameter name to update (e.g. `"b"`).
    int width;                   ///< Bit width consumed by the perturbation step.
    /**
     * @brief Shadow of the transformation's current bitmask value.
     *
     * Updated in place by the driver and readable by the caller after
     * the driver returns.
     */
    int64_t current_value;
};

/**
 * @brief Configuration for the tempering bitmask optimiser.
 *
 * @ingroup core
 */
struct TemperingOptimizerConfig {
    int max_essais;              ///< Perturbation budget for one `run_once()` call.
    /**
     * @brief Per-resolution gap budget; size `L + 1` (`delta[0]` unused).
     *
     * Set every entry to `INT_MAX` for "minimise without caps".
     */
    std::vector<int> delta;
    int mse;                     ///< Upper bound on the cumulative equidistribution defect.
    /**
     * @brief Number of `run_once` restarts.
     *
     * `1` runs one `run_once` only; `> 1` enables iterative
     * delta-tightening via `run_minimize`.
     */
    int n_restarts;
    uint64_t random_seed;        ///< Seed for perturbation; `0` = OS entropy.
};

/**
 * @brief Result reported by `run_tempering_optimizer_*`.
 *
 * @ingroup core
 */
struct TemperingOptResult {
    int se;                      ///< Sum of gaps for `v = 1..L`.
    std::vector<int> gaps;       ///< Per-resolution gaps (size `L + 1`; `gaps[0]` unused).
    double elapsed_seconds;      ///< Wall-clock time.
    /**
     * @brief Iterations consumed.
     *
     * For `run_once`: number of perturbation iterations. For
     * `run_minimize`: number of `run_once` calls.
     */
    int essais;
};

/**
 * @brief Run a single recursive bitmask optimisation pass.
 *
 * `safe_masks` is shaped `[L+1][P]` where `safe_masks[v][p]` is the
 * bitmask of bits that are "safe" to perturb in parameter `p` at
 * resolution `v`. `safe_masks[0]` is unused.
 *
 * On entry, `params[p].current_value` must reflect the live value of
 * `trans->update(...)` so the perturbation step can XOR against it.
 * On return, the locators carry the best-found values and the
 * underlying `Transformation` objects have those values applied.
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   TemperingOptimizerConfig cfg;
 *   cfg.max_essais = 400;
 *   cfg.delta.assign(L + 1, 0);   // ME mode (zero-defect target)
 *   cfg.mse = 0;
 *   cfg.n_restarts = 1;
 *   auto result = run_tempering_optimizer_once(cfg, cache, params, masks);
 * @endcode
 *
 * @param cfg         Optimiser configuration.
 * @param cache       Stateful tempering cache (one per `Combination`).
 * @param params      Locators for every optimisable bitmask parameter.
 * @param safe_masks  Per-resolution per-parameter safe-bit masks.
 * @return            Best-found `(se, gaps)` plus runtime metadata.
 *
 * @see :py:class:`regpoly.search.tempering_optimizer.TemperingOptimizer`
 */
TemperingOptResult run_tempering_optimizer_once(
    const TemperingOptimizerConfig& cfg,
    TemperOptCache& cache,
    std::vector<TemperParamLocator>& params,
    const std::vector<std::vector<uint64_t>>& safe_masks);

/**
 * @brief Iterative delta-tightening loop (`n_restarts > 1`).
 *
 * Re-randomises every locator's `current_value`, runs `run_once` with
 * the previous best as the new `(delta, mse)` cap, accepts on
 * improvement, counts consecutive non-improvements as failures, and
 * stops when `failures == n_restarts`.
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   TemperingOptimizerConfig cfg;
 *   cfg.max_essais = 400;
 *   cfg.delta.assign(L + 1, INT_MAX);   // unrestricted; we will tighten
 *   cfg.mse = INT_MAX;
 *   cfg.n_restarts = 8;
 *   auto result = run_tempering_optimizer_minimize(cfg, cache, params, masks);
 * @endcode
 *
 * @param cfg         Optimiser configuration (`n_restarts > 1`).
 * @param cache       Stateful tempering cache (one per `Combination`).
 * @param params      Locators for every optimisable bitmask parameter.
 * @param safe_masks  Per-resolution per-parameter safe-bit masks.
 * @return            Best-found `(se, gaps)` plus runtime metadata.
 *
 * @see :py:class:`regpoly.search.tempering_optimizer.TemperingOptimizer`
 */
TemperingOptResult run_tempering_optimizer_minimize(
    const TemperingOptimizerConfig& cfg,
    TemperOptCache& cache,
    std::vector<TemperParamLocator>& params,
    const std::vector<std::vector<uint64_t>>& safe_masks);

}  // namespace regpoly::core
