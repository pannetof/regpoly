// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.

// Unit tests for WELLGen::random_matrices, the cost-bounded sampler.

#include <gtest/gtest.h>

#include <random>
#include <set>
#include <stdexcept>
#include <variant>

#include "well.h"

namespace {

constexpr int kW = 32;

// Sum the per-Mi costs across the 8 slots of a sampled matrices map.
int total_cost(const StructMap& m) {
    int sum = 0;
    for (const auto& kv : m) {
        const auto& e = kv.second;
        auto it = e.find("M");
        if (it == e.end()) continue;
        int Mi = static_cast<int>(std::get<int64_t>(it->second));
        sum += WELLGen::static_cost_for_Mi(Mi);
    }
    return sum;
}

// Read the M-class index from a slot entry.
int Mi_of(const StructEntry& e) {
    return static_cast<int>(std::get<int64_t>(e.at("M")));
}

}  // namespace

// ── Cap is honoured ─────────────────────────────────────────────────────

class WellRandomMatricesCap : public testing::TestWithParam<int> {};

TEST_P(WellRandomMatricesCap, EveryDrawHonoursTheCap) {
    int cap = GetParam();
    std::mt19937_64 rng{42 + static_cast<uint64_t>(cap)};
    for (int i = 0; i < 1000; ++i) {
        StructMap m = WELLGen::random_matrices(kW, cap, rng);
        ASSERT_EQ(m.size(), 8u);
        EXPECT_LE(total_cost(m), cap)
            << "cap=" << cap << " sample=" << i;
    }
}

INSTANTIATE_TEST_SUITE_P(Caps, WellRandomMatricesCap,
                          testing::Values(1, 4, 12, 32, 64));

// ── Class-membership properties at tight caps ──────────────────────────

TEST(WellRandomMatrices, MaxCost1OnlyM0AndM1) {
    // costs: {0,1,2,3,5,4,8} for M0..M6.
    // Cap of 1 budgets at most 1 unit per slot total — only Mi with
    // cost ≤ 1 are allowed: M0 (0) and M1 (1). And the sum must be ≤ 1
    // → at most one slot is M1, the rest are M0.
    std::mt19937_64 rng{12345};
    for (int i = 0; i < 1000; ++i) {
        StructMap m = WELLGen::random_matrices(kW, 1, rng);
        for (const auto& kv : m) {
            int Mi = Mi_of(kv.second);
            ASSERT_TRUE(Mi == 0 || Mi == 1)
                << "cap=1 sample=" << i << " slot=" << kv.first
                << " unexpectedly produced M" << Mi;
        }
        EXPECT_LE(total_cost(m), 1);
    }
}

TEST(WellRandomMatrices, MaxCost4ExcludesM4M6) {
    // Cap 4 forbids M4 (cost 5) and M6 (cost 8). A single M5 (cost 4)
    // exhausts the budget; otherwise mix of M0..M3 + (optionally M5).
    std::mt19937_64 rng{67890};
    for (int i = 0; i < 1000; ++i) {
        StructMap m = WELLGen::random_matrices(kW, 4, rng);
        for (const auto& kv : m) {
            int Mi = Mi_of(kv.second);
            ASSERT_NE(Mi, 4) << "cap=4 sample=" << i;
            ASSERT_NE(Mi, 6) << "cap=4 sample=" << i;
        }
        EXPECT_LE(total_cost(m), 4);
    }
}

TEST(WellRandomMatrices, MaxCost64AllMiClassesAppearAtLeastOnce) {
    // Cap 64 = 8 × cost(M6). Rejection sampling has 100% accept rate
    // here. Across 10000 samples × 8 slots = 80000 Mi draws, every Mi
    // 0..6 must appear at least once. (Probability of missing any is
    // negligible: (6/7)^80000 ≈ 0.)
    std::mt19937_64 rng{2026};
    std::set<int> seen;
    for (int i = 0; i < 10000 && (int)seen.size() < 7; ++i) {
        StructMap m = WELLGen::random_matrices(kW, 64, rng);
        for (const auto& kv : m) seen.insert(Mi_of(kv.second));
    }
    for (int Mi = 0; Mi < 7; ++Mi)
        EXPECT_EQ(seen.count(Mi), 1u) << "M" << Mi << " never sampled";
}

// ── Per-Mi arg ranges (paper Table I) ──────────────────────────────────

TEST(WellRandomMatrices, ArgsWithinPaperRanges) {
    std::mt19937_64 rng{2027};
    for (int i = 0; i < 1000; ++i) {
        StructMap m = WELLGen::random_matrices(kW, 64, rng);
        for (const auto& kv : m) {
            const auto& e = kv.second;
            int Mi = Mi_of(e);
            switch (Mi) {
                case 0:
                case 1:
                    // No args expected beyond M.
                    EXPECT_EQ(e.size(), 1u);
                    break;
                case 2:
                case 3: {
                    EXPECT_EQ(e.size(), 2u);
                    int t = static_cast<int>(std::get<int64_t>(e.at("t")));
                    EXPECT_GE(std::abs(t), 1);
                    EXPECT_LE(std::abs(t), kW - 1);
                    break;
                }
                case 4: {
                    EXPECT_EQ(e.size(), 2u);
                    uint64_t a = std::get<uint64_t>(e.at("a"));
                    EXPECT_GE(a, 1u);
                    EXPECT_LE(a, 0xFFFFFFFFu);
                    break;
                }
                case 5: {
                    EXPECT_EQ(e.size(), 3u);
                    int t = static_cast<int>(std::get<int64_t>(e.at("t")));
                    EXPECT_GE(std::abs(t), 1);
                    uint64_t b = std::get<uint64_t>(e.at("b"));
                    EXPECT_GE(b, 1u);
                    break;
                }
                case 6: {
                    EXPECT_EQ(e.size(), 5u);
                    int q = static_cast<int>(std::get<int64_t>(e.at("q")));
                    int t = static_cast<int>(std::get<int64_t>(e.at("t")));
                    int s = static_cast<int>(std::get<int64_t>(e.at("s")));
                    uint64_t a = std::get<uint64_t>(e.at("a"));
                    EXPECT_GE(q, 0); EXPECT_LE(q, kW - 1);
                    EXPECT_GE(t, 0); EXPECT_LE(t, kW - 1);
                    EXPECT_GE(s, 0); EXPECT_LE(s, kW - 1);
                    EXPECT_GE(a, 1u);
                    break;
                }
            }
        }
    }
}

// ── Validation ─────────────────────────────────────────────────────────

TEST(WellRandomMatrices, NonPositiveCapThrows) {
    std::mt19937_64 rng{1};
    EXPECT_THROW(WELLGen::random_matrices(kW, 0, rng),
                 std::invalid_argument);
    EXPECT_THROW(WELLGen::random_matrices(kW, -1, rng),
                 std::invalid_argument);
}

TEST(WellRandomMatrices, NonStandardWordSizeThrows) {
    std::mt19937_64 rng{1};
    EXPECT_THROW(WELLGen::random_matrices(/*w=*/16, 12, rng),
                 std::invalid_argument);
}

// ── Round-trips through the construction path ─────────────────────────

// from_params() must accept every sampled config — i.e. the sampler's
// per-Mi arg ranges agree with decode_matrix_entry's validator.
TEST(WellRandomMatrices, EverySamplePassesFromParamsValidator) {
    std::mt19937_64 rng{2028};
    for (int i = 0; i < 100; ++i) {
        StructMap m = WELLGen::random_matrices(kW, 32, rng);
        Params p;
        p.set_int("w", 32);
        p.set_int("r", 16);
        p.set_int("p", 0);
        p.set_int("m1", 13);
        p.set_int("m2", 9);
        p.set_int("m3", 5);
        p.set_struct_map("matrices", m);
        // from_params throws std::out_of_range or std::runtime_error
        // if any per-Mi arg violates decode_matrix_entry's bounds.
        EXPECT_NO_THROW({
            auto gen = WELLGen::from_params(p, 32);
            EXPECT_LE(static_cast<WELLGen*>(gen.get())->total_cost(), 32);
        }) << "sample " << i << " failed validator";
    }
}
