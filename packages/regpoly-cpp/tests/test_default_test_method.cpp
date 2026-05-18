// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

// Tests for Generator::default_test_method(test_type).
//
// Coverage matrix:
//   - SIMD families (SFMTGen, DSFMTGen, MTGPGen, RMT64Gen) override and
//     return "simd_notprimitive" for "equidistribution", nullopt otherwise.
//   - CombinedGenerator overrides and returns "notprimitive" for
//     "equidistribution", nullopt otherwise.
//   - The base implementation handles non-SIMD scalar generators via the
//     full-period + k-threshold rule:
//       * full-period, k <= 100  -> "matricial"
//       * full-period, k >  100  -> "harase"
//       * not full-period        -> "notprimitive"

#include <gtest/gtest.h>

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "combined.h"
#include "dsfmt.h"
#include "equidistribution_method.h"
#include "generator.h"
#include "mtgp.h"
#include "rmt64.h"
#include "sfmt.h"
#include "tausworthe.h"

using namespace regpoly::core;


namespace {

// SFMT-19937 reference parameters (from SFMT-params19937.h).
std::unique_ptr<SFMTGen> make_sfmt19937() {
    return std::make_unique<SFMTGen>(
        /*mexp=*/19937, /*pos1=*/122,
        /*sl1=*/18, /*sl2=*/1, /*sr1=*/11, /*sr2=*/1,
        /*msk1=*/0xdfffffefu, /*msk2=*/0xddfecb7fu,
        /*msk3=*/0xbffaffffu, /*msk4=*/0xbffffff6u,
        /*L=*/32);
}

// dSFMT-19937 reference parameters (from dSFMT-params19937.h).
std::unique_ptr<DSFMTGen> make_dsfmt19937() {
    return std::make_unique<DSFMTGen>(
        /*mexp=*/19937, /*pos1=*/117, /*sl1=*/19,
        /*msk1=*/0x000ffafffffffb3full, /*msk2=*/0x000ffdfffc90fffdull,
        /*L=*/52);
}

// MTGP-23209 reference parameters (mtgp32-23209.c, idx 0). Only the
// shape matters for default_test_method; the recursion table content
// is irrelevant. tbl/tmp_tbl must be 16 entries each.
std::unique_ptr<MTGPGen> make_mtgp() {
    std::vector<uint32_t> tbl(16, 0);
    std::vector<uint32_t> tmp_tbl(16, 0);
    return std::make_unique<MTGPGen>(
        /*mexp=*/23209, /*pos=*/0, /*sh1=*/0, /*sh2=*/0,
        /*mask=*/0u, tbl, tmp_tbl, /*L=*/32);
}

// RMT64 dummy parameters; default_test_method does not consult them.
std::unique_ptr<RMT64Gen> make_rmt64() {
    return std::make_unique<RMT64Gen>(
        /*mexp=*/521, /*pos=*/8,
        /*mata=*/0ull, /*maskb=*/0ull, /*maskc=*/0ull, /*L=*/64);
}

// Full-period TauswortheGen with k=31 (Mersenne prime exponent), trinomial
// x^31 + x^3 + 1 — the canonical small full-period fixture used elsewhere
// in the test suite.
std::unique_ptr<TauswortheGen> make_taus_k31() {
    return std::make_unique<TauswortheGen>(
        31, std::vector<int>{0, 3, 31}, /*s=*/13,
        /*quicktaus=*/true, /*L=*/32);
}

// Full-period TauswortheGen with k=127 (Mersenne prime exponent), primitive
// trinomial x^127 + x + 1. quicktaus requires k <= L; with L=32 we use the
// general-purpose recurrence (quicktaus=false).
std::unique_ptr<TauswortheGen> make_taus_k127() {
    return std::make_unique<TauswortheGen>(
        127, std::vector<int>{0, 1, 127}, /*s=*/1,
        /*quicktaus=*/false, /*L=*/32);
}

// k=31 reducible polynomial: x^31 + 1 = (x+1)(x^30 + x^29 + ... + 1).
// Used to exercise the not-full-period base-class branch.
std::unique_ptr<TauswortheGen> make_taus_reducible() {
    return std::make_unique<TauswortheGen>(
        31, std::vector<int>{0, 31}, /*s=*/1,
        /*quicktaus=*/false, /*L=*/32);
}

}  // namespace

// ── SIMD families ─────────────────────────────────────────────────────────

TEST(DefaultTestMethod, SfmtReturnsSimdNotprimitive) {
    auto g = make_sfmt19937();
    EXPECT_EQ(g->default_test_method("equidistribution"),
              std::optional<std::string>{"simd_notprimitive"});
    EXPECT_FALSE(g->default_test_method("collision_free").has_value());
    EXPECT_FALSE(g->default_test_method("tuplets").has_value());
    EXPECT_FALSE(g->default_test_method("nonsense").has_value());
}

TEST(DefaultTestMethod, DsfmtReturnsSimdNotprimitive) {
    auto g = make_dsfmt19937();
    EXPECT_EQ(g->default_test_method("equidistribution"),
              std::optional<std::string>{"simd_notprimitive"});
    EXPECT_FALSE(g->default_test_method("collision_free").has_value());
}

TEST(DefaultTestMethod, MtgpReturnsSimdNotprimitive) {
    auto g = make_mtgp();
    EXPECT_EQ(g->default_test_method("equidistribution"),
              std::optional<std::string>{"simd_notprimitive"});
    EXPECT_FALSE(g->default_test_method("tuplets").has_value());
}

TEST(DefaultTestMethod, Rmt64ReturnsNotprimitive) {
    // RMT64 has reducible chi but no SIMD lane structure; the right
    // method is "notprimitive", not "simd_notprimitive".
    auto g = make_rmt64();
    EXPECT_EQ(g->default_test_method("equidistribution"),
              std::optional<std::string>{"notprimitive"});
    EXPECT_FALSE(g->default_test_method("collision_free").has_value());
}

// ── CombinedGenerator ─────────────────────────────────────────────────────

TEST(DefaultTestMethod, CombinedReturnsNotprimitive) {
    std::vector<std::unique_ptr<Generator>> comps;
    comps.push_back(make_taus_k31());
    comps.push_back(make_taus_k127());
    CombinedGenerator combined(std::move(comps), /*Lmax=*/32);
    EXPECT_EQ(combined.default_test_method("equidistribution"),
              std::optional<std::string>{"notprimitive"});
    EXPECT_FALSE(combined.default_test_method("collision_free").has_value());
    EXPECT_FALSE(combined.default_test_method("tuplets").has_value());
    EXPECT_FALSE(combined.default_test_method("foo").has_value());
}

// JEquals1 invariant: a CombinedGenerator wrapping a single component
// must be observationally equal to that component for default_test_method.
TEST(DefaultTestMethod, CombinedJEquals1MatchesComponent) {
    std::vector<std::unique_ptr<Generator>> comps;
    comps.push_back(make_taus_k31());
    CombinedGenerator combined(std::move(comps), /*Lmax=*/32);
    // Wrapped Taus k=31 is full-period and k <= 100 -> "matricial".
    EXPECT_EQ(combined.default_test_method("equidistribution"),
              std::optional<std::string>{"matricial"});
    EXPECT_FALSE(combined.default_test_method("collision_free").has_value());
}

// ── Base class — full-period, k <= 100 ───────────────────────────────────

TEST(DefaultTestMethod, FullPeriodSmallReturnsMatricial) {
    auto g = make_taus_k31();
    ASSERT_TRUE(g->is_full_period());
    ASSERT_LE(g->k(), 100);
    EXPECT_EQ(g->default_test_method("equidistribution"),
              std::optional<std::string>{"matricial"});
    EXPECT_FALSE(g->default_test_method("collision_free").has_value());
}

// ── Base class — full-period, k > 100 ────────────────────────────────────

TEST(DefaultTestMethod, FullPeriodLargeReturnsHarase) {
    auto g = make_taus_k127();
    ASSERT_TRUE(g->is_full_period());
    ASSERT_GT(g->k(), 100);
    EXPECT_EQ(g->default_test_method("equidistribution"),
              std::optional<std::string>{"harase"});
}

// ── Base class — not full-period ─────────────────────────────────────────

TEST(DefaultTestMethod, NotFullPeriodReturnsNotprimitive) {
    auto g = make_taus_reducible();
    ASSERT_FALSE(g->is_full_period());
    EXPECT_EQ(g->default_test_method("equidistribution"),
              std::optional<std::string>{"notprimitive"});
}

// ── Cache: repeated calls must return the same answer fast and not
// recompute is_full_period(). We probe this by calling many times in
// a tight loop on a generator whose runtime path goes through BM.
TEST(DefaultTestMethod, RepeatedCallsAreConsistent) {
    auto g = make_taus_k127();
    auto first = g->default_test_method("equidistribution");
    ASSERT_EQ(first, std::optional<std::string>{"harase"});
    for (int i = 0; i < 100; ++i) {
        EXPECT_EQ(g->default_test_method("equidistribution"), first);
        EXPECT_FALSE(g->default_test_method("collision_free").has_value());
    }
}

// Instrumented cache verification: a TauswortheGen subclass that
// counts every entry into compute_default_test_method. Repeated
// public default_test_method() calls must invoke the compute path
// at most ONCE per (test_type) — that's what the cache buys.
namespace {
class CountingTaus : public TauswortheGen {
public:
    using TauswortheGen::TauswortheGen;
    mutable int equid_calls = 0;
    mutable int other_calls = 0;
protected:
    std::optional<std::string>
    compute_default_test_method(const std::string& test_type) const override {
        if (test_type == "equidistribution") ++equid_calls;
        else                                  ++other_calls;
        return TauswortheGen::compute_default_test_method(test_type);
    }
};
}  // namespace

// ── R3 extensibility check: registering a new EquidistributionMethod
// without touching seek_search.cpp dispatch. ──────────────────────────

namespace {
class FakeMethod : public EquidistributionMethod {
public:
    std::string name() const override { return "fake_test_method"; }
    EquidistributionMethodResult run(
        const Generator&, int, int, int maxL,
        const std::vector<int>&, int) const override
    {
        // Distinctive marker so the test can verify *this* impl ran.
        return {std::vector<int>(maxL + 1, 42), /*se=*/99,
                /*verified=*/true};
    }
};
}

TEST(MethodRegistry, RegisteringNewMethodWorksWithoutDispatchChanges) {
    // Pre: method is unknown.
    EXPECT_FALSE(MethodRegistry::has("fake_test_method"));

    MethodRegistry::reg("fake_test_method",
        []{ return std::unique_ptr<EquidistributionMethod>(new FakeMethod); });

    // Post: registry knows the new method.
    EXPECT_TRUE(MethodRegistry::has("fake_test_method"));

    // The factory builds an instance and run() returns FakeMethod's
    // distinctive result — no edits to seek_search.cpp, no enum entry,
    // no dispatch branch.
    auto m = MethodRegistry::create("fake_test_method");
    EXPECT_EQ(m->name(), "fake_test_method");

    auto g = make_taus_k31();
    auto r = m->run(*g, /*kg=*/31, /*L=*/32, /*maxL=*/8,
                    std::vector<int>(9, 0), /*mse=*/0);
    EXPECT_EQ(r.se, 99);
    EXPECT_TRUE(r.verified);
    EXPECT_EQ(r.ecart.size(), 9u);
    for (int v : r.ecart) EXPECT_EQ(v, 42);
}

TEST(DefaultTestMethod, CacheCallsComputeOnce) {
    CountingTaus g(31, std::vector<int>{0, 3, 31}, /*s=*/13,
                   /*quicktaus=*/true, /*L=*/32);
    EXPECT_EQ(g.equid_calls, 0);
    EXPECT_EQ(g.other_calls, 0);

    // First call: compute runs.
    auto v = g.default_test_method("equidistribution");
    EXPECT_EQ(v, std::optional<std::string>{"matricial"});
    EXPECT_EQ(g.equid_calls, 1);

    // Subsequent calls for the same key: cache must short-circuit.
    for (int i = 0; i < 50; ++i) {
        EXPECT_EQ(g.default_test_method("equidistribution"), v);
    }
    EXPECT_EQ(g.equid_calls, 1) << "compute should run exactly once";

    // Different test_type goes through compute once as well.
    EXPECT_FALSE(g.default_test_method("collision_free").has_value());
    EXPECT_FALSE(g.default_test_method("collision_free").has_value());
    EXPECT_FALSE(g.default_test_method("tuplets").has_value());
    EXPECT_EQ(g.other_calls, 2)
        << "one compute per distinct test_type, regardless of repeats";
}
