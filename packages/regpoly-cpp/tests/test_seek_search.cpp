// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

// Phase 2.4b (TDD): SeekDriver — equidistribution search loop in C++.

#include <gtest/gtest.h>

#include <climits>
#include <vector>

#include "combination.h"
#include "factory.h"
#include "params.h"
#include "seek_search.h"

using namespace regpoly::core;


namespace {

std::unique_ptr<Generator> make_tgfsr(uint32_t a) {
    Params p;
    p.set_int("w", 32);
    p.set_int("r", 3);
    p.set_int("m", 1);
    p.set_int("a", static_cast<int64_t>(a));
    return create_generator("TGFSRGen", p, /*L=*/32);
}

SeekTestSpec make_matricial(int Lmax_test, int mse) {
    SeekTestSpec s;
    s.kind = SeekTestKind::EquidistributionMatricial;
    s.eq_L_max_test = Lmax_test;
    s.eq_delta.assign(Lmax_test + 1, INT_MAX);
    s.eq_mse = mse;
    return s;
}

}  // namespace

TEST(SeekDriver, EnumeratesEveryCombo) {
    Combination c(2, 32);
    c.component(0).add_gen(*make_tgfsr(1));
    c.component(0).add_gen(*make_tgfsr(2));
    c.component(1).add_gen(*make_tgfsr(3));
    c.component(1).add_gen(*make_tgfsr(4));
    ASSERT_TRUE(c.reset());

    std::vector<SeekTestSpec> tests{make_matricial(32, INT_MAX)};
    int64_t calls = 0;
    auto on_iter = [&](Combination&, const SeekIterResult& r) {
        EXPECT_TRUE(r.me_ran);
        EXPECT_EQ(static_cast<int>(r.me_ecart.size()), 33);
        ++calls;
    };
    auto on_progress = [](const SeekProgress&) {};

    auto result = run_seek_search(
        c, tests, /*nbtries=*/1, /*progress_interval=*/100,
        nullptr, on_iter, on_progress);

    // Cartesian: 2 * 2 = 4 combos, all selected (mse=INT_MAX).
    EXPECT_EQ(result.nbgen, 4);
    EXPECT_EQ(result.nb_select, 4);
    EXPECT_EQ(calls, 4);
}

TEST(SeekDriver, NbtriesRetriesEveryComboNTimes) {
    Combination c(1, 32);
    c.component(0).add_gen(*make_tgfsr(1));
    c.component(0).add_gen(*make_tgfsr(2));
    c.component(0).add_gen(*make_tgfsr(3));
    ASSERT_TRUE(c.reset());

    std::vector<SeekTestSpec> tests{make_matricial(32, INT_MAX)};
    auto on_iter = [](Combination&, const SeekIterResult&) {};
    auto on_progress = [](const SeekProgress&) {};

    auto result = run_seek_search(
        c, tests, /*nbtries=*/3, /*progress_interval=*/1,
        nullptr, on_iter, on_progress);

    // 3 combos * 3 retries = 9 iterations.
    EXPECT_EQ(result.nbgen, 9);
}

TEST(SeekDriver, MseRejectionShortCircuits) {
    Combination c(1, 32);
    c.component(0).add_gen(*make_tgfsr(1));
    ASSERT_TRUE(c.reset());

    // mse=0 means we accept ONLY full-equidistribution combos. With one
    // iteration the test runs and either selects or doesn't.
    std::vector<SeekTestSpec> tests{make_matricial(32, /*mse=*/0)};
    int selected = 0;
    auto on_iter = [&](Combination&, const SeekIterResult&) { ++selected; };
    auto on_progress = [](const SeekProgress&) {};

    auto result = run_seek_search(
        c, tests, 1, 1, nullptr, on_iter, on_progress);

    EXPECT_EQ(result.nbgen, 1);
    EXPECT_LE(result.nb_select, 1);
    EXPECT_EQ(selected, result.nb_select);
}

TEST(SeekDriver, ProgressFiresOnInterval) {
    Combination c(1, 32);
    for (int a : {1, 2, 3, 4, 5, 6, 7, 8, 9, 10})
        c.component(0).add_gen(*make_tgfsr(a));
    ASSERT_TRUE(c.reset());

    std::vector<SeekTestSpec> tests{make_matricial(32, INT_MAX)};
    auto on_iter = [](Combination&, const SeekIterResult&) {};

    std::vector<int64_t> progress_nbgens;
    auto on_progress = [&](const SeekProgress& p) {
        progress_nbgens.push_back(p.nbgen);
    };

    auto result = run_seek_search(
        c, tests, 1, /*progress_interval=*/3,
        nullptr, on_iter, on_progress);

    EXPECT_EQ(result.nbgen, 10);
    // Progress fires at 3, 6, 9 (multiples of 3 inside the loop).
    ASSERT_GE(progress_nbgens.size(), 3u);
    EXPECT_EQ(progress_nbgens[0], 3);
    EXPECT_EQ(progress_nbgens[1], 6);
    EXPECT_EQ(progress_nbgens[2], 9);
}

TEST(SeekDriver, OnPrepFiresPerIterationWithRetryFlag) {
    Combination c(1, 32);
    c.component(0).add_gen(*make_tgfsr(1));
    c.component(0).add_gen(*make_tgfsr(2));
    ASSERT_TRUE(c.reset());

    std::vector<SeekTestSpec> tests{make_matricial(32, INT_MAX)};
    std::vector<bool> retry_flags;
    auto on_prep = [&](Combination&, bool is_retry) {
        retry_flags.push_back(is_retry);
    };
    auto on_iter = [](Combination&, const SeekIterResult&) {};
    auto on_progress = [](const SeekProgress&) {};

    auto result = run_seek_search(
        c, tests, /*nbtries=*/2, 1, on_prep, on_iter, on_progress);

    // 2 combos * 2 retries = 4 iterations. Pattern of is_retry:
    //   combo1 try1 = false, combo1 try2 = true,
    //   combo2 try1 = false, combo2 try2 = true.
    EXPECT_EQ(result.nbgen, 4);
    ASSERT_EQ(retry_flags.size(), 4u);
    EXPECT_FALSE(retry_flags[0]);
    EXPECT_TRUE(retry_flags[1]);
    EXPECT_FALSE(retry_flags[2]);
    EXPECT_TRUE(retry_flags[3]);
}

TEST(SeekDriver, WorksWithNoEquidistTest) {
    // No tests => nothing is selected, but iteration still happens.
    Combination c(1, 32);
    c.component(0).add_gen(*make_tgfsr(1));
    c.component(0).add_gen(*make_tgfsr(2));
    ASSERT_TRUE(c.reset());

    std::vector<SeekTestSpec> tests{};
    int selected = 0;
    auto on_iter = [&](Combination&, const SeekIterResult&) { ++selected; };
    auto on_progress = [](const SeekProgress&) {};

    auto result = run_seek_search(
        c, tests, 1, 1, nullptr, on_iter, on_progress);

    EXPECT_EQ(result.nbgen, 2);
    EXPECT_EQ(result.nb_select, 0);
    EXPECT_EQ(selected, 0);
}
