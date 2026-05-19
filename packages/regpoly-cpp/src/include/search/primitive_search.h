// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "params.h"
#include "search_types.h"

#include <cstdint>
#include <functional>
#include <string>

/**
 * @file primitive_search.h
 * @brief Full-period (primitive characteristic polynomial) parameter-search driver.
 * @ingroup core
 *
 * Phase 2.4: full-period (primitive char-poly) search driver.
 *
 * Replaces the inner loop of
 * `packages/regpoly/src/regpoly/search/search_primitive.py`:
 *
 * @verbatim
 *   while not (max_tries reached or max_seconds reached):
 *       full_params = randomize_missing(structural + fixed, family_specs)
 *       gen = create_generator(family, full_params, L)        # may throw
 *       if gen.is_full_period():
 *           emit TestedGenerator{family, gen.k(), L, full_params}
 *       periodically emit SearchProgress{tries, elapsed_seconds}
 * @endverbatim
 *
 * Python remains responsible for: YAML parsing, `.partial` recovery
 * I/O, dedup against the existing-results set, the merged output
 * file, header/footer printing. The driver owns the random sampling
 * (via `regpoly::random::sample_param_into`), the factory call, and
 * the `is_full_period` check.
 *
 * @see :py:class:`regpoly.search.search_primitive.PrimitiveSearch`
 */

namespace regpoly::core {

/**
 * @brief Driver configuration for `run_primitive_search`.
 *
 * Captures the family + L of the search, any user-pinned structural
 * or fixed parameters, the stop conditions (try count or wall clock),
 * the progress cadence, and (for WELL families) the optional cost
 * ceiling that gates the per-iteration matrix sampler.
 *
 * @ingroup core
 */
struct PrimitiveSearchConfig {
    std::string family;        ///< Canonical generator family name.
    int L;                     ///< Output word width in bits.
    Params structural_params;  ///< User-pinned structural inputs (define `k`).
    Params fixed_params;       ///< User-pinned non-structural inputs.
    int64_t max_tries;         ///< Max iterations; `0` = no try-count limit.
    double max_seconds;        ///< Wall-time budget; `0.0` = no time limit.
    int progress_interval;     ///< Emit `on_progress` every N tries (default 100).
    uint64_t random_seed;      ///< Seed for sampling; `0` = OS entropy.
    /**
     * @brief WELL-only: upper bound on the sum of per-Mi costs across
     *        the 8 algorithm slots T0..T7.
     *
     * `0` disables the cap. When `> 0` and `family` is a WELL alias,
     * the search loop overrides `matrices` per-iteration with a
     * freshly-sampled config whose total cost <= `max_cost`. Throws
     * `std::invalid_argument` at search start if the user has also
     * pinned `matrices` in `fixed_params` (pin + search are mutually
     * exclusive).
     */
    int max_cost = 0;
};

/// Per-hit callback signature consumed by `run_primitive_search`.
using OnHitFn = std::function<void(const TestedGenerator&)>;
/// Periodic progress callback signature consumed by `run_primitive_search`.
using OnProgressFn = std::function<void(const SearchProgress&)>;

/**
 * @brief Run the full-period parameter search loop.
 *
 * Iterates until `max_tries` or `max_seconds` is exhausted. Each
 * iteration samples missing parameters, builds the candidate
 * generator via the factory, and emits a `TestedGenerator` to
 * `on_hit` when the candidate's characteristic polynomial is
 * primitive. Progress is reported through `on_progress` every
 * `cfg.progress_interval` tries.
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   PrimitiveSearchConfig cfg;
 *   cfg.family = "MTGen";
 *   cfg.L = 32;
 *   cfg.max_tries = 10000;
 *   cfg.progress_interval = 100;
 *   run_primitive_search(
 *       cfg,
 *       [](const TestedGenerator& g) {
 *           // Persist or display g.params.
 *       },
 *       [](const SearchProgress& p) {
 *           // Periodic progress.
 *       });
 * @endcode
 *
 * @param cfg          Search configuration.
 * @param on_hit       Callback fired for every accepted candidate.
 * @param on_progress  Periodic progress callback.
 * @return             The total number of tries actually executed.
 *
 * @see :py:class:`regpoly.search.search_primitive.PrimitiveSearch`
 */
int64_t run_primitive_search(
    const PrimitiveSearchConfig& cfg,
    const OnHitFn& on_hit,
    const OnProgressFn& on_progress);

}  // namespace regpoly::core
