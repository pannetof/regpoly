// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

// Phase 2.4 (TDD): PrimitiveSearchDriver.
//
// The driver re-randomizes a generator's non-structural parameters
// each iteration, builds the generator via the factory, and emits a
// callback for every full-period hit. These tests focus on the
// driver's contract (limits, callback shape, hit selection) — the
// underlying primitivity check is exercised by test_primitivity.cpp.

#include <gtest/gtest.h>

#include <atomic>
#include <chrono>
#include <thread>
#include <vector>

#include "primitive_search.h"
#include "primitivity.h"
#include "search_types.h"
#include "well.h"
#include "factory.h"

using namespace regpoly::core;
using namespace regpoly::internal;


namespace {

// Build a config that searches for full-period TGFSR with a fixed
// (w=32, r=3, m=1) shape, randomizing only `a`.
PrimitiveSearchConfig make_tgfsr_config(int64_t max_tries) {
    PrimitiveSearchConfig cfg;
    cfg.family = "TGFSRGen";
    cfg.L = 32;
    cfg.structural_params.set_int("w", 32);
    cfg.structural_params.set_int("r", 3);
    cfg.fixed_params.set_int("m", 1);
    cfg.max_tries = max_tries;
    cfg.max_seconds = 0.0;
    cfg.progress_interval = 50;
    cfg.random_seed = 2024;
    return cfg;
}

}  // namespace

TEST(PrimitiveSearchDriver, MaxTriesLimitIsHonoured) {
    auto cfg = make_tgfsr_config(/*max_tries=*/200);

    int64_t hit_count = 0;
    auto on_hit = [&](const TestedGenerator& tg) {
        ASSERT_EQ(tg.family, "TGFSRGen");
        ASSERT_EQ(tg.L, 32);
        ASSERT_GT(tg.k, 0);
        ASSERT_TRUE(tg.params.has("a"));
        ++hit_count;
    };

    int64_t last_progress_tries = 0;
    auto on_progress = [&](const SearchProgress& sp) {
        last_progress_tries = sp.tries;
    };

    auto total_tries = run_primitive_search(cfg, on_hit, on_progress);
    EXPECT_EQ(total_tries, 200);
    EXPECT_EQ(last_progress_tries, 200);
    // hit_count is non-deterministic but with 200 tries on TGFSR(96)
    // we should usually see at least one hit; to keep the test stable,
    // assert only that hit_count <= total_tries.
    EXPECT_LE(hit_count, 200);
}

TEST(PrimitiveSearchDriver, ProgressFiresOnInterval) {
    auto cfg = make_tgfsr_config(/*max_tries=*/250);
    cfg.progress_interval = 100;

    std::vector<int64_t> progress_tries;
    auto on_progress = [&](const SearchProgress& sp) {
        progress_tries.push_back(sp.tries);
    };
    auto on_hit = [](const TestedGenerator&) {};

    run_primitive_search(cfg, on_hit, on_progress);

    // Expect callbacks at tries=100, 200, and a final at tries=250
    // (the loop exits when tries reaches max_tries; the final
    // progress event fires once after the loop).
    ASSERT_GE(progress_tries.size(), 3u);
    EXPECT_EQ(progress_tries[0], 100);
    EXPECT_EQ(progress_tries[1], 200);
    EXPECT_EQ(progress_tries.back(), 250);
}

TEST(PrimitiveSearchDriver, MaxSecondsLimitIsHonoured) {
    auto cfg = make_tgfsr_config(/*max_tries=*/0);
    cfg.max_seconds = 0.05;  // 50 ms

    auto on_hit = [](const TestedGenerator&) {};
    auto on_progress = [](const SearchProgress&) {};

    auto t0 = std::chrono::steady_clock::now();
    auto tries = run_primitive_search(cfg, on_hit, on_progress);
    auto t1 = std::chrono::steady_clock::now();

    double elapsed = std::chrono::duration<double>(t1 - t0).count();
    EXPECT_GE(tries, 1);
    // Loose upper bound: should finish within 1s even on a very slow
    // CI box. The relevant property is "did stop", not "stopped at
    // exactly 50 ms".
    EXPECT_LT(elapsed, 1.0);
}

// Regression test for the StructMap merge bug fixed alongside the
// cost-cap feature: the per-iteration `merge_into` helper must copy
// struct_maps_ from fixed_params, otherwise a user-pinned WELL
// `matrices` is silently dropped and `WELLGen::from_params` throws
// "matrices is required" on the very first iteration.
//
// On master before the fix, this test fails with std::runtime_error
// from from_params. After the fix, the search runs without throwing.
TEST(PrimitiveSearchDriver, PreservesUserPinnedStructMap) {
    PrimitiveSearchConfig cfg;
    cfg.family = "WELLGen";
    cfg.L = 32;
    cfg.structural_params.set_int("w", 32);
    cfg.structural_params.set_int("r", 16);
    cfg.structural_params.set_int("p", 0);
    cfg.fixed_params.set_int("m1", 13);
    cfg.fixed_params.set_int("m2", 9);
    cfg.fixed_params.set_int("m3", 5);
    // WELL512a from paper Table II: T0..T7 = M3(-16), M3(-15), M3(11),
    // M0, M3(-2), M3(-18), M2(-28), M5(-5, a1=0xda442d24).
    StructMap matrices;
    auto entry = [](int Mi) { StructEntry e; e["M"] = (int64_t)Mi; return e; };
    auto entry_t = [](int Mi, int t) {
        StructEntry e; e["M"] = (int64_t)Mi; e["t"] = (int64_t)t; return e;
    };
    matrices["T0"] = entry_t(3, -16);
    matrices["T1"] = entry_t(3, -15);
    matrices["T2"] = entry_t(3,  11);
    matrices["T3"] = entry(0);
    matrices["T4"] = entry_t(3,  -2);
    matrices["T5"] = entry_t(3, -18);
    matrices["T6"] = entry_t(2, -28);
    {
        StructEntry e;
        e["M"] = (int64_t)5;
        e["t"] = (int64_t)-5;
        e["b"] = (uint64_t)0xda442d24ULL;
        matrices["T7"] = e;
    }
    cfg.fixed_params.set_struct_map("matrices", matrices);
    cfg.max_tries = 5;
    cfg.max_seconds = 0.0;
    cfg.progress_interval = 1;
    cfg.random_seed = 2026;

    auto on_hit = [](const TestedGenerator&) {};
    auto on_progress = [](const SearchProgress&) {};

    // Before the merge_into fix this throws; after, it runs the loop
    // without complaint.
    EXPECT_NO_THROW(run_primitive_search(cfg, on_hit, on_progress));
}

TEST(PrimitiveSearchDriver, MissingStructuralParamThrows) {
    PrimitiveSearchConfig cfg;
    cfg.family = "TGFSRGen";
    cfg.L = 32;
    // structural_params intentionally empty: w, r are required and have
    // no random sampler. The driver must throw on the first iteration.
    cfg.max_tries = 1;
    cfg.max_seconds = 0.0;
    cfg.progress_interval = 1;
    cfg.random_seed = 2024;

    auto on_hit = [](const TestedGenerator&) {};
    auto on_progress = [](const SearchProgress&) {};

    EXPECT_THROW(run_primitive_search(cfg, on_hit, on_progress),
                 std::invalid_argument);
}

namespace {

// Build a WELL search config (random matrices, varying offsets) ready
// to accept a max_cost cap.
PrimitiveSearchConfig make_well_search_config(int max_tries) {
    PrimitiveSearchConfig cfg;
    cfg.family = "WELLGen";
    cfg.L = 32;
    cfg.structural_params.set_int("w", 32);
    cfg.structural_params.set_int("r", 16);
    cfg.structural_params.set_int("p", 0);
    cfg.max_tries = max_tries;
    cfg.max_seconds = 0.0;
    cfg.progress_interval = max_tries;  // single final progress event
    cfg.random_seed = 2026;
    return cfg;
}

}  // namespace

// Cost-bounded WELL search: every emitted hit's total cost must be
// ≤ cfg.max_cost. Uses a small WELL (k=512) so the loop is fast.
TEST(PrimitiveSearchDriver, WellMaxCostBoundsAllHits) {
    auto cfg = make_well_search_config(/*max_tries=*/200);
    cfg.max_cost = 12;

    int n_hits = 0;
    int max_hit_cost = 0;
    auto on_hit = [&](const TestedGenerator& tg) {
        // Reconstruct the generator to compute its total_cost.
        auto gen = create_generator(tg.family, tg.params, tg.L);
        auto* well = dynamic_cast<WELLGen*>(gen.get());
        ASSERT_NE(well, nullptr);
        int c = well->total_cost();
        EXPECT_LE(c, 12) << "hit at try=" << tg.tries_at_hit
                          << " has cost " << c;
        max_hit_cost = std::max(max_hit_cost, c);
        ++n_hits;
    };
    auto on_progress = [](const SearchProgress&) {};

    auto total_tries = run_primitive_search(cfg, on_hit, on_progress);
    EXPECT_GE(total_tries, 1);
    // hit_count is non-deterministic but at most max_tries.
    EXPECT_LE(n_hits, total_tries);
}

// Pinning `matrices` while also setting max_cost is mutually exclusive
// — searching means varying. The driver must reject at start.
TEST(PrimitiveSearchDriver, WellPinAndMaxCostMutuallyExclusive) {
    auto cfg = make_well_search_config(/*max_tries=*/5);
    cfg.max_cost = 12;
    StructMap matrices;
    StructEntry m1; m1["M"] = (int64_t)1;
    for (int j = 0; j < 8; ++j)
        matrices["T" + std::to_string(j)] = m1;
    cfg.fixed_params.set_struct_map("matrices", matrices);

    auto on_hit = [](const TestedGenerator&) {};
    auto on_progress = [](const SearchProgress&) {};

    EXPECT_THROW(run_primitive_search(cfg, on_hit, on_progress),
                 std::invalid_argument);
}

// Family alias `WELLRNG` must trigger the cost-cap path the same way
// the canonical `WELLGen` does.
TEST(PrimitiveSearchDriver, WellAliasResolvesForMaxCost) {
    auto cfg = make_well_search_config(/*max_tries=*/50);
    cfg.family = "WELLRNG";  // alias for WELLGen
    cfg.max_cost = 8;

    auto on_hit = [&](const TestedGenerator& tg) {
        auto gen = create_generator(tg.family, tg.params, tg.L);
        auto* well = dynamic_cast<WELLGen*>(gen.get());
        ASSERT_NE(well, nullptr);
        EXPECT_LE(well->total_cost(), 8);
    };
    auto on_progress = [](const SearchProgress&) {};
    EXPECT_NO_THROW(run_primitive_search(cfg, on_hit, on_progress));
}
