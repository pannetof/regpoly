// Phase 2.4c (TDD): TemperingSearchDriver — per-combo / per-try search loop.
//
// This driver iterates a Combination, fires a Python-supplied on_try
// callback nb_tries times per combo (callback re-randomizes tempering
// params, optionally invokes the optimizer, runs a test, returns the
// score), tracks the best score per combo, and emits on_combo_start /
// on_combo_done / on_progress callbacks. Test invocation lives in the
// callback because randomize_params currently lives on Python
// Transformations; everything else is C++-driven loop scaffolding.

#include <gtest/gtest.h>

#include <climits>
#include <vector>

#include "combination.h"
#include "factory.h"
#include "params.h"
#include "tempering_search.h"

namespace {

std::unique_ptr<Generator> make_tgfsr(uint32_t a) {
    Params p;
    p.set_int("w", 32);
    p.set_int("r", 3);
    p.set_int("m", 1);
    p.set_int("a", static_cast<int64_t>(a));
    return create_generator("TGFSRGen", p, /*L=*/32);
}

}  // namespace

TEST(TemperingSearchDriver, IteratesEveryComboNbTriesTimes) {
    Combination c(2, 32);
    c.component(0).add_gen(*make_tgfsr(1));
    c.component(0).add_gen(*make_tgfsr(2));
    c.component(1).add_gen(*make_tgfsr(3));
    c.component(1).add_gen(*make_tgfsr(4));
    ASSERT_TRUE(c.reset());

    TemperingSearchConfig cfg;
    cfg.nb_tries = 5;
    cfg.progress_interval = 0;

    int combo_starts = 0;
    int try_calls = 0;
    int combo_dones = 0;

    auto on_combo_start = [&](Combination&, int) { ++combo_starts; };
    auto on_try = [&](Combination&, int, int, bool) {
        ++try_calls;
        TemperingTryOutcome out;
        out.got_result = true;
        out.score = 1;        // constant score; never breaks early
        return out;
    };
    auto on_combo_done = [&](Combination&, int, int, int) { ++combo_dones; };

    auto result = run_tempering_search(
        c, cfg, on_combo_start, on_try, on_combo_done, nullptr);

    // Cartesian: 2 * 2 = 4 combos; each runs 5 tries.
    EXPECT_EQ(result.nbgen, 4);
    EXPECT_EQ(result.nb_with_result, 4);
    EXPECT_EQ(combo_starts, 4);
    EXPECT_EQ(try_calls, 4 * 5);
    EXPECT_EQ(combo_dones, 4);
}

TEST(TemperingSearchDriver, EarlyStopOnZeroScore) {
    Combination c(1, 32);
    c.component(0).add_gen(*make_tgfsr(7));
    ASSERT_TRUE(c.reset());

    TemperingSearchConfig cfg;
    cfg.nb_tries = 100;
    cfg.progress_interval = 0;

    int try_calls = 0;
    int reported_best_score = -1;
    int reported_best_try = -1;

    auto on_try = [&](Combination&, int, int t, bool) {
        ++try_calls;
        TemperingTryOutcome out;
        out.got_result = true;
        // Return score 0 on the third try -> driver should break early.
        out.score = (t == 2) ? 0 : 5;
        return out;
    };
    auto on_combo_done =
        [&](Combination&, int, int best_score, int best_try) {
            reported_best_score = best_score;
            reported_best_try = best_try;
        };

    auto result = run_tempering_search(
        c, cfg, nullptr, on_try, on_combo_done, nullptr);

    EXPECT_EQ(result.nbgen, 1);
    EXPECT_EQ(try_calls, 3);          // stopped after the zero
    EXPECT_EQ(reported_best_score, 0);
    EXPECT_EQ(reported_best_try, 2);
}

TEST(TemperingSearchDriver, TracksMinimumScorePerCombo) {
    Combination c(1, 32);
    c.component(0).add_gen(*make_tgfsr(1));
    c.component(0).add_gen(*make_tgfsr(2));
    ASSERT_TRUE(c.reset());

    TemperingSearchConfig cfg;
    cfg.nb_tries = 4;
    cfg.progress_interval = 0;

    // Hand-rolled scores per (combo, try). Combo 1: [9, 3, 7, 4]; min=3.
    // Combo 2: [8, 5, 2, 6]; min=2.
    int combo_seen = 0;
    std::vector<int> per_try_scores_combo1 = {9, 3, 7, 4};
    std::vector<int> per_try_scores_combo2 = {8, 5, 2, 6};

    std::vector<int> best_scores;
    std::vector<int> best_tries;

    auto on_combo_start = [&](Combination&, int) { ++combo_seen; };
    auto on_try = [&](Combination&, int, int t, bool) {
        TemperingTryOutcome out;
        out.got_result = true;
        out.score = (combo_seen == 1) ? per_try_scores_combo1[t]
                                       : per_try_scores_combo2[t];
        return out;
    };
    auto on_combo_done =
        [&](Combination&, int, int best_score, int best_try) {
            best_scores.push_back(best_score);
            best_tries.push_back(best_try);
        };

    auto result = run_tempering_search(
        c, cfg, on_combo_start, on_try, on_combo_done, nullptr);

    EXPECT_EQ(result.nbgen, 2);
    ASSERT_EQ(best_scores.size(), 2u);
    EXPECT_EQ(best_scores[0], 3);
    EXPECT_EQ(best_tries[0], 1);
    EXPECT_EQ(best_scores[1], 2);
    EXPECT_EQ(best_tries[1], 2);
}

TEST(TemperingSearchDriver, FailedTriesNotCounted) {
    Combination c(1, 32);
    c.component(0).add_gen(*make_tgfsr(1));
    ASSERT_TRUE(c.reset());

    TemperingSearchConfig cfg;
    cfg.nb_tries = 3;
    cfg.progress_interval = 0;

    bool combo_done_fired = false;
    int reported_best_try = -99;

    auto on_try = [&](Combination&, int, int, bool) {
        TemperingTryOutcome out;
        out.got_result = false;     // every try fails
        out.score = 0;
        return out;
    };
    auto on_combo_done =
        [&](Combination&, int, int, int best_try) {
            combo_done_fired = true;
            reported_best_try = best_try;
        };

    auto result = run_tempering_search(
        c, cfg, nullptr, on_try, on_combo_done, nullptr);

    EXPECT_EQ(result.nbgen, 1);
    EXPECT_EQ(result.nb_with_result, 0);
    // on_combo_done is NOT fired when no try produced a result.
    EXPECT_FALSE(combo_done_fired);
    EXPECT_EQ(reported_best_try, -99);
}

TEST(TemperingSearchDriver, IsFirstFlagOnFirstTryOfEachCombo) {
    Combination c(1, 32);
    c.component(0).add_gen(*make_tgfsr(1));
    c.component(0).add_gen(*make_tgfsr(2));
    ASSERT_TRUE(c.reset());

    TemperingSearchConfig cfg;
    cfg.nb_tries = 3;
    cfg.progress_interval = 0;

    std::vector<bool> is_first_flags;
    auto on_try = [&](Combination&, int, int, bool is_first) {
        is_first_flags.push_back(is_first);
        TemperingTryOutcome out;
        out.got_result = true;
        out.score = 1;
        return out;
    };
    auto on_combo_done = [&](Combination&, int, int, int) {};

    run_tempering_search(c, cfg, nullptr, on_try, on_combo_done, nullptr);

    // 2 combos * 3 tries = 6 calls. is_first true at indices 0, 3.
    ASSERT_EQ(is_first_flags.size(), 6u);
    EXPECT_TRUE(is_first_flags[0]);
    EXPECT_FALSE(is_first_flags[1]);
    EXPECT_FALSE(is_first_flags[2]);
    EXPECT_TRUE(is_first_flags[3]);
    EXPECT_FALSE(is_first_flags[4]);
    EXPECT_FALSE(is_first_flags[5]);
}

TEST(TemperingSearchDriver, ProgressFiresOnInterval) {
    Combination c(1, 32);
    for (uint32_t a : {1u, 2u, 3u, 4u, 5u}) {
        c.component(0).add_gen(*make_tgfsr(a));
    }
    ASSERT_TRUE(c.reset());

    TemperingSearchConfig cfg;
    cfg.nb_tries = 4;
    cfg.progress_interval = 6;       // every 6 tries

    auto on_try = [](Combination&, int, int, bool) {
        TemperingTryOutcome out;
        out.got_result = true;
        out.score = 1;
        return out;
    };
    auto on_combo_done = [](Combination&, int, int, int) {};

    std::vector<int64_t> progress_tries;
    auto on_progress = [&](const SearchProgress& p) {
        progress_tries.push_back(p.tries);
    };

    auto result = run_tempering_search(
        c, cfg, nullptr, on_try, on_combo_done, on_progress);

    EXPECT_EQ(result.nbgen, 5);
    // Total tries: 5 * 4 = 20. Progress fires at multiples of 6: 6, 12, 18.
    ASSERT_GE(progress_tries.size(), 3u);
    EXPECT_EQ(progress_tries[0], 6);
    EXPECT_EQ(progress_tries[1], 12);
    EXPECT_EQ(progress_tries[2], 18);
}
