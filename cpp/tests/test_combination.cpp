// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

// Phase 2.4b-pre (TDD): F2LinearSourcePool + ComboEnumerator iteration in C++.
//
// Mirrors the semantics of regpoly.core.combination. Reference
// expectations were captured from the Python product_no_repeat_ordered
// algorithm on small generator pools.

#include <gtest/gtest.h>

#include <memory>
#include <set>
#include <vector>

#include "combined_f2_linear_source.h"
#include "combined_f2_linear_source.h"
#include "combo_enumerator.h"
#include "factory.h"
#include "params.h"
#include "sobol.h"

using namespace regpoly::core;


namespace {

// Build a tiny TGFSR generator with a chosen `a` so we can produce
// distinguishable, fast-to-construct fixtures.
std::unique_ptr<Recurrence> make_tgfsr(uint32_t a) {
    ParamBag p;
    p.set_int("w", 32);
    p.set_int("r", 3);
    p.set_int("m", 1);
    p.set_int("a", static_cast<int64_t>(a));
    return create_generator("TGFSRGen", p, /*L=*/32);
}

}  // namespace

TEST(F2LinearSourcePool, AddGenAndTransGrowPools) {
    F2LinearSourcePool c;
    EXPECT_EQ(c.nb_sources(), 0);
    EXPECT_EQ(c.nb_trans(), 0);

    auto g = make_tgfsr(0xdeadbeef);
    c.add_source(*g);
    EXPECT_EQ(c.nb_sources(), 1);
}

TEST(F2LinearSourcePool, GenAtRespectsPoolOrder) {
    F2LinearSourcePool c;
    auto g1 = make_tgfsr(1);
    auto g2 = make_tgfsr(2);
    c.add_source(*g1);
    c.add_source(*g2);
    EXPECT_EQ(c.nb_sources(), 2);

    F2LinearSource& a = c.source_at(0);
    F2LinearSource& b = c.source_at(1);
    EXPECT_EQ(a.k(), b.k());
    EXPECT_NE(&a, &b);
}

TEST(F2LinearSourcePool, SharePoolGivesPointerEquality) {
    F2LinearSourcePool owner;
    F2LinearSourcePool shared;
    owner.add_source(*make_tgfsr(1));
    owner.add_source(*make_tgfsr(2));

    shared.copy_pool_from(owner);
    EXPECT_EQ(owner.pool_id(), shared.pool_id());
    EXPECT_EQ(shared.nb_sources(), 2);
    EXPECT_EQ(&owner.source_at(0), &shared.source_at(0));
}

TEST(F2LinearSourcePool, AddGenAfterShareThrows) {
    F2LinearSourcePool owner;
    F2LinearSourcePool shared;
    owner.add_source(*make_tgfsr(1));
    shared.copy_pool_from(owner);

    EXPECT_THROW(shared.add_source(*make_tgfsr(2)), std::logic_error);
}

TEST(F2LinearSourcePool, AcceptsDigitalNet) {
    // Phase 5.2-storage: storage widened to F2LinearSource so digital
    // nets can enter the pool (not just Recurrence PRNGs).
    SobolNet net(/*m=*/2, /*s_max=*/4);
    F2LinearSourcePool p;
    p.add_source(net);
    EXPECT_EQ(p.nb_sources(), 1);

    F2LinearSource& s = p.source_at(0);
    EXPECT_EQ(s.name(), net.name());
    EXPECT_EQ(s.k(), net.k());
}

TEST(ComboEnumerator, BuildProducesCombinedF2LinearSourceForRecurrencePool) {
    ComboEnumerator c(1, 32);
    c.pool(0).add_source(*make_tgfsr(1));
    ASSERT_TRUE(c.reset());
    auto built = build_combined_from_enumerator(c);
    ASSERT_NE(built.get(), nullptr);
    EXPECT_NE(dynamic_cast<CombinedF2LinearSource*>(built.get()), nullptr);
}

TEST(ComboEnumerator, BuildAcceptsDigitalNetInPool) {
    // The combined wrapper holds F2LinearSource components, so digital
    // nets compose freely. Kernels that need Recurrence components
    // (matricial χ-recovery, SIMD path) throw at call time if a
    // DigitalNet is present.
    ComboEnumerator c(1, 32);
    c.pool(0).add_source(SobolNet(/*m=*/2, /*s_max=*/4));
    ASSERT_TRUE(c.reset());
    auto built = build_combined_from_enumerator(c);
    ASSERT_NE(built.get(), nullptr);
    EXPECT_NE(dynamic_cast<CombinedF2LinearSource*>(built.get()), nullptr);
}

TEST(ComboEnumerator, EmptyCombHasNoValidCombo) {
    ComboEnumerator c(/*J=*/2, /*Lmax=*/32);
    EXPECT_FALSE(c.reset());
    EXPECT_TRUE(c.exhausted());
}

TEST(ComboEnumerator, JEquals1IteratesEveryGenerator) {
    ComboEnumerator c(1, 32);
    c.pool(0).add_source(*make_tgfsr(1));
    c.pool(0).add_source(*make_tgfsr(2));
    c.pool(0).add_source(*make_tgfsr(3));

    ASSERT_TRUE(c.reset());
    int count = 1;
    while (c.next()) ++count;
    EXPECT_EQ(count, 3);
    EXPECT_TRUE(c.exhausted());
}

TEST(ComboEnumerator, JEquals2IndependentPoolsCartesian) {
    ComboEnumerator c(2, 32);
    c.pool(0).add_source(*make_tgfsr(1));
    c.pool(0).add_source(*make_tgfsr(2));
    c.pool(1).add_source(*make_tgfsr(3));
    c.pool(1).add_source(*make_tgfsr(4));

    ASSERT_TRUE(c.reset());
    int count = 1;
    while (c.next()) ++count;
    // 2 * 2 = 4 combos; no shared-pool constraint and no identity collision.
    EXPECT_EQ(count, 4);
}

TEST(ComboEnumerator, JEquals2SharedPoolEnforcesCnk) {
    ComboEnumerator c(2, 32);
    c.pool(0).add_source(*make_tgfsr(1));
    c.pool(0).add_source(*make_tgfsr(2));
    c.pool(0).add_source(*make_tgfsr(3));
    c.pool(0).add_source(*make_tgfsr(4));
    c.pool(1).copy_pool_from(c.pool(0));

    ASSERT_TRUE(c.reset());
    int count = 1;
    while (c.next()) ++count;
    // C(4, 2) = 6 combos. No (i, i) self-pairs, and (i, j) with j > i only.
    EXPECT_EQ(count, 6);
}

TEST(ComboEnumerator, KgIsSumAndLIsMinCappedAtLmax) {
    ComboEnumerator c(2, 16);  // Lmax intentionally below the gens' L=32
    c.pool(0).add_source(*make_tgfsr(1));
    c.pool(1).add_source(*make_tgfsr(2));

    ASSERT_TRUE(c.reset());
    EXPECT_EQ(c.k_g(), c.at(0).k() + c.at(1).k());
    EXPECT_EQ(c.L(), 16);  // min(32, 32) = 32, capped at Lmax=16
}

TEST(ComboEnumerator, AtReturnsActiveGeneratorAfterAdvance) {
    ComboEnumerator c(1, 32);
    c.pool(0).add_source(*make_tgfsr(1));
    c.pool(0).add_source(*make_tgfsr(2));

    ASSERT_TRUE(c.reset());
    F2LinearSource* first = &c.at(0);
    ASSERT_TRUE(c.next());
    F2LinearSource* second = &c.at(0);
    EXPECT_NE(first, second);
}

TEST(ComboEnumerator, NextAfterExhaustionStaysFalse) {
    ComboEnumerator c(1, 32);
    c.pool(0).add_source(*make_tgfsr(1));

    ASSERT_TRUE(c.reset());
    EXPECT_FALSE(c.next());
    EXPECT_TRUE(c.exhausted());
    EXPECT_FALSE(c.next());  // still false on a second call
}

TEST(ComboEnumerator, ResetAfterExhaustionReinitializes) {
    ComboEnumerator c(1, 32);
    c.pool(0).add_source(*make_tgfsr(1));
    c.pool(0).add_source(*make_tgfsr(2));

    ASSERT_TRUE(c.reset());
    EXPECT_TRUE(c.next());
    EXPECT_FALSE(c.next());
    EXPECT_TRUE(c.exhausted());

    // Reset re-arms.
    ASSERT_TRUE(c.reset());
    EXPECT_FALSE(c.exhausted());
    EXPECT_TRUE(c.next());
}

TEST(ComboEnumerator, JEquals3MixedSharedAndIndependent) {
    ComboEnumerator c(3, 32);
    // pool A: 3 gens, slots 0 and 2 share it.
    c.pool(0).add_source(*make_tgfsr(1));
    c.pool(0).add_source(*make_tgfsr(2));
    c.pool(0).add_source(*make_tgfsr(3));
    // slot 1: independent pool of 2 gens.
    c.pool(1).add_source(*make_tgfsr(11));
    c.pool(1).add_source(*make_tgfsr(12));
    // slot 2: shares pool with slot 0.
    c.pool(2).copy_pool_from(c.pool(0));

    ASSERT_TRUE(c.reset());
    int count = 1;
    while (c.next()) ++count;
    // C(3, 2) * 2 = 3 * 2 = 6 combos.
    EXPECT_EQ(count, 6);
}
