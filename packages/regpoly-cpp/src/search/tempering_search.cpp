// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "tempering_search.h"

#include <chrono>
#include <climits>

TemperingSearchResult run_tempering_search(
    Combination& comb,
    const TemperingSearchConfig& cfg,
    const TempSearchOnComboStartFn& on_combo_start,
    const TempSearchOnTryFn& on_try,
    const TempSearchOnComboDoneFn& on_combo_done,
    const TempSearchOnProgressFn& on_progress)
{
    const int nb_tries = (cfg.nb_tries < 1) ? 1 : cfg.nb_tries;
    const int progress_interval = cfg.progress_interval;

    int64_t nbgen = 0;
    int64_t nb_with_result = 0;
    int64_t total_tries = 0;
    int combo_idx = 0;
    const auto t_start = std::chrono::steady_clock::now();

    while (true) {
        ++combo_idx;
        ++nbgen;

        if (on_combo_start) on_combo_start(comb, combo_idx);

        int best_score = INT_MAX;
        int best_try = -1;

        for (int t = 0; t < nb_tries; ++t) {
            const bool is_first = (t == 0);
            TemperingTryOutcome outcome =
                on_try(comb, combo_idx, t, is_first);
            ++total_tries;

            if (outcome.got_result) {
                if (best_try < 0 || outcome.score < best_score) {
                    best_score = outcome.score;
                    best_try = t;
                }
            }

            if (on_progress && progress_interval > 0
                && (total_tries % progress_interval == 0)) {
                const auto now = std::chrono::steady_clock::now();
                SearchProgress p;
                p.tries = total_tries;
                p.elapsed_seconds =
                    std::chrono::duration<double>(now - t_start).count();
                on_progress(p);
            }

            if (best_try >= 0 && best_score == 0) break;
        }

        if (best_try >= 0) {
            ++nb_with_result;
            if (on_combo_done) {
                on_combo_done(comb, combo_idx, best_score, best_try);
            }
        }

        if (!comb.next()) break;
    }

    const auto t_end = std::chrono::steady_clock::now();
    TemperingSearchResult result;
    result.nbgen = nbgen;
    result.nb_with_result = nb_with_result;
    result.elapsed_seconds =
        std::chrono::duration<double>(t_end - t_start).count();
    return result;
}
