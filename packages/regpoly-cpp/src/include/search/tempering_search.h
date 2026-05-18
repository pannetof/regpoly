// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "combination.h"
#include "search_types.h"

#include <cstdint>
#include <functional>

// Phase 2.4c: tempering-search driver in C++.
//
// Replaces the per-combo / per-try loop body of
// regpoly.search.tempering_search.TemperingSearch.run() with a pure-C++
// loop that drives a Combination and fires Python callbacks for the
// per-try work. The per-try work (re-randomizing tempering parameters
// via Transformation::randomize_params, optionally invoking the
// optimizer, running the test) stays in Python because randomize_params
// currently lives on Python Transformation objects. The C++ side owns:
//
//   * Combination iteration (reset / next),
//   * the nb_tries inner loop with best-score tracking,
//   * early-exit on score == 0,
//   * combo lifecycle callbacks (on_combo_start / on_combo_done),
//   * progress reporting at fixed intervals.
//
// on_combo_done fires only when at least one try in the combo
// succeeded (got_result == true). This matches Python's behaviour of
// skipping combos with no usable results.

namespace regpoly::core {

struct TemperingSearchConfig {
    int nb_tries = 1;
    int progress_interval = 0;       // 0 = never fire on_progress
};

// One try's outcome. score is what the driver minimizes (lower = better).
// got_result == false means the try produced nothing usable; the driver
// ignores its score and continues to the next try.
struct TemperingTryOutcome {
    bool got_result = false;
    int score = 0;
};

struct TemperingSearchResult {
    int64_t nbgen = 0;               // total combos iterated
    int64_t nb_with_result = 0;      // combos with >=1 successful try
    double elapsed_seconds = 0.0;
};

using TempSearchOnComboStartFn =
    std::function<void(Combination&, int /*combo_idx*/)>;
using TempSearchOnTryFn =
    std::function<TemperingTryOutcome(
        Combination&, int /*combo_idx*/, int /*try_idx*/, bool /*is_first*/)>;
using TempSearchOnComboDoneFn =
    std::function<void(Combination&, int /*combo_idx*/,
                       int /*best_score*/, int /*best_try*/)>;
using TempSearchOnProgressFn =
    std::function<void(const SearchProgress&)>;

// Run the tempering-search loop. `comb` must be reset()-able by the
// caller (e.g. via comb.reset() returning true). The loop:
//
//   for each combination:
//     fire on_combo_start(comb, combo_idx)         # if non-null
//     best_score = INT_MAX; best_try = -1
//     for t in 0..nb_tries-1:
//       outcome = on_try(comb, combo_idx, t, t==0)
//       fire on_progress every progress_interval tries  # if interval > 0
//       if outcome.got_result and outcome.score < best_score:
//         best_score = outcome.score; best_try = t
//         if best_score == 0: break
//     if best_try >= 0:
//       fire on_combo_done(comb, combo_idx, best_score, best_try)
//   advance via comb.next()
//
// Returns the loop-level totals.
TemperingSearchResult run_tempering_search(
    Combination& comb,
    const TemperingSearchConfig& cfg,
    const TempSearchOnComboStartFn& on_combo_start,
    const TempSearchOnTryFn& on_try,
    const TempSearchOnComboDoneFn& on_combo_done,
    const TempSearchOnProgressFn& on_progress);

}  // namespace regpoly::core
