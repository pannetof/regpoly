// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "combination.h"
#include "search_types.h"

#include <cstdint>
#include <functional>

/**
 * @file tempering_search.h
 * @brief Full-stack combined-generator tempering search driver.
 * @ingroup core
 *
 * Phase 2.4c: tempering-search driver in C++.
 *
 * Replaces the per-combo / per-try loop body of
 * `regpoly.search.tempering_search.TemperingSearch.run()` with a
 * pure-C++ loop that drives a `Combination` and fires Python
 * callbacks for the per-try work. The per-try work (re-randomising
 * tempering parameters via `Transformation::randomize_params`,
 * optionally invoking the optimiser, running the test) stays in
 * Python because `randomize_params` currently lives on Python
 * `Transformation` objects. The C++ side owns:
 *
 *   - `Combination` iteration (`reset` / `next`),
 *   - the `nb_tries` inner loop with best-score tracking,
 *   - early exit on `score == 0`,
 *   - combo lifecycle callbacks (`on_combo_start` /
 *     `on_combo_done`),
 *   - progress reporting at fixed intervals.
 *
 * The *generator* pool comes from the catalog or inline YAML, while
 * the *tempering* pool is the parameter space the optimiser
 * searches. `on_combo_done` fires only when at least one try in the
 * combo succeeded (`got_result == true`). This matches Python's
 * behaviour of skipping combos with no usable results.
 *
 * @see :py:class:`regpoly.search.tempering_search.TemperingSearch`
 */

namespace regpoly::core {

/**
 * @brief Configuration for `run_tempering_search`.
 *
 * @ingroup core
 */
struct TemperingSearchConfig {
    int nb_tries = 1;                ///< Number of tries per combo.
    int progress_interval = 0;       ///< Tries between `on_progress` callbacks; `0` = never fire.
};

/**
 * @brief One try's outcome reported by `on_try`.
 *
 * `score` is what the driver minimises (lower = better). `got_result
 * == false` means the try produced nothing usable; the driver ignores
 * its score and continues to the next try.
 *
 * @ingroup core
 */
struct TemperingTryOutcome {
    bool got_result = false;     ///< True iff the try produced a usable score.
    int score = 0;               ///< Score (lower is better).
};

/**
 * @brief Loop-level totals returned by `run_tempering_search`.
 *
 * @ingroup core
 */
struct TemperingSearchResult {
    int64_t nbgen = 0;               ///< Total combos iterated.
    int64_t nb_with_result = 0;      ///< Combos with >= 1 successful try.
    double elapsed_seconds = 0.0;    ///< Wall-clock time.
};

/// Combo-start callback. Fires once at the start of each combo.
using TempSearchOnComboStartFn =
    std::function<void(Combination&, int /*combo_idx*/)>;
/// Per-try callback. Fires once per try inside each combo; must return the try's outcome.
using TempSearchOnTryFn =
    std::function<TemperingTryOutcome(
        Combination&, int /*combo_idx*/, int /*try_idx*/, bool /*is_first*/)>;
/// Combo-done callback. Fires once at the end of each combo that had >= 1 successful try.
using TempSearchOnComboDoneFn =
    std::function<void(Combination&, int /*combo_idx*/,
                       int /*best_score*/, int /*best_try*/)>;
/// Periodic progress callback.
using TempSearchOnProgressFn =
    std::function<void(const SearchProgress&)>;

/**
 * @brief Run the tempering-search loop over a configured `Combination`.
 *
 * `comb` must be `reset()`-able by the caller (e.g. via
 * `comb.reset()` returning true). The loop:
 *
 * @verbatim
 *   for each combination:
 *     fire on_combo_start(comb, combo_idx)         # if non-null
 *     best_score = INT_MAX; best_try = -1
 *     for t in 0..nb_tries-1:
 *       outcome = on_try(comb, combo_idx, t, t==0)
 *       fire on_progress every progress_interval tries  # if interval > 0
 *       if outcome.got_result and outcome.score < best_score:
 *         best_score = outcome.score; best_try = t
 *         if best_score == 0: break
 *     if best_try >= 0:
 *       fire on_combo_done(comb, combo_idx, best_score, best_try)
 *   advance via comb.next()
 * @endverbatim
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   TemperingSearchConfig cfg;
 *   cfg.nb_tries = 16;
 *   cfg.progress_interval = 100;
 *   run_tempering_search(comb, cfg,
 *       [](Combination&, int) { },
 *       [](Combination&, int, int, bool) -> TemperingTryOutcome {
 *           // ... run the per-try work (randomise + optimise + test) ...
 *           return {true, 0};  // got_result = true, score = 0
 *       },
 *       [](Combination&, int, int, int) { },
 *       [](const SearchProgress&) { });
 * @endcode
 *
 * @param comb            Configured iterator (must be `reset()`-able).
 * @param cfg             Search configuration.
 * @param on_combo_start  Optional per-combo entry callback.
 * @param on_try          Per-try callback returning the outcome.
 * @param on_combo_done   Optional per-combo done callback (only on success).
 * @param on_progress     Periodic progress callback.
 * @return                Loop-level totals.
 *
 * @see :py:class:`regpoly.search.tempering_search.TemperingSearch`
 */
TemperingSearchResult run_tempering_search(
    Combination& comb,
    const TemperingSearchConfig& cfg,
    const TempSearchOnComboStartFn& on_combo_start,
    const TempSearchOnTryFn& on_try,
    const TempSearchOnComboDoneFn& on_combo_done,
    const TempSearchOnProgressFn& on_progress);

}  // namespace regpoly::core
