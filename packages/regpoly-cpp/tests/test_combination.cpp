// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

// Phase 2.4b-pre (TDD): Component + Combination iteration in C++.
//
// Mirrors the semantics of regpoly.core.combination. Reference
// expectations were captured from the Python product_no_repeat_ordered
// algorithm on small generator pools.

#include <gtest/gtest.h>

#include <memory>
#include <set>
#include <vector>

#include "combination.h"
#include "factory.h"
#include "params.h"

namespace {

// Build a tiny TGFSR generator with a chosen `a` so we can produce
// distinguishable, fast-to-construct fixtures.
std::unique_ptr<Generator> make_tgfsr(uint32_t a) {
    Params p;
    p.set_int("w", 32);
    p.set_int("r", 3);
    p.set_int("m", 1);
    p.set_int("a", static_cast<int64_t>(a));
    return create_generator("TGFSRGen", p, /*L=*/32);
}

}  // namespace

TEST(Component, AddGenAndTransGrowPools) {
    Component c;
    EXPECT_EQ(c.nb_gen(), 0);
    EXPECT_EQ(c.nb_trans(), 0);

    auto g = make_tgfsr(0xdeadbeef);
    c.add_gen(*g);
    EXPECT_EQ(c.nb_gen(), 1);
}

TEST(Component, GenAtRespectsPoolOrder) {
    Component c;
    auto g1 = make_tgfsr(1);
    auto g2 = make_tgfsr(2);
    c.add_gen(*g1);
    c.add_gen(*g2);
    EXPECT_EQ(c.nb_gen(), 2);

    Generator& a = c.gen_at(0);
    Generator& b = c.gen_at(1);
    EXPECT_EQ(a.k(), b.k());
    EXPECT_NE(&a, &b);
}

TEST(Component, SharePoolGivesPointerEquality) {
    Component owner;
    Component shared;
    owner.add_gen(*make_tgfsr(1));
    owner.add_gen(*make_tgfsr(2));

    shared.share_pool_with(owner);
    EXPECT_EQ(owner.pool_id(), shared.pool_id());
    EXPECT_EQ(shared.nb_gen(), 2);
    EXPECT_EQ(&owner.gen_at(0), &shared.gen_at(0));
}

TEST(Component, AddGenAfterShareThrows) {
    Component owner;
    Component shared;
    owner.add_gen(*make_tgfsr(1));
    shared.share_pool_with(owner);

    EXPECT_THROW(shared.add_gen(*make_tgfsr(2)), std::logic_error);
}

TEST(Combination, EmptyCombHasNoValidCombo) {
    Combination c(/*J=*/2, /*Lmax=*/32);
    EXPECT_FALSE(c.reset());
    EXPECT_TRUE(c.exhausted());
}

TEST(Combination, JEquals1IteratesEveryGenerator) {
    Combination c(1, 32);
    c.component(0).add_gen(*make_tgfsr(1));
    c.component(0).add_gen(*make_tgfsr(2));
    c.component(0).add_gen(*make_tgfsr(3));

    ASSERT_TRUE(c.reset());
    int count = 1;
    while (c.next()) ++count;
    EXPECT_EQ(count, 3);
    EXPECT_TRUE(c.exhausted());
}

TEST(Combination, JEquals2IndependentPoolsCartesian) {
    Combination c(2, 32);
    c.component(0).add_gen(*make_tgfsr(1));
    c.component(0).add_gen(*make_tgfsr(2));
    c.component(1).add_gen(*make_tgfsr(3));
    c.component(1).add_gen(*make_tgfsr(4));

    ASSERT_TRUE(c.reset());
    int count = 1;
    while (c.next()) ++count;
    // 2 * 2 = 4 combos; no shared-pool constraint and no identity collision.
    EXPECT_EQ(count, 4);
}

TEST(Combination, JEquals2SharedPoolEnforcesCnk) {
    Combination c(2, 32);
    c.component(0).add_gen(*make_tgfsr(1));
    c.component(0).add_gen(*make_tgfsr(2));
    c.component(0).add_gen(*make_tgfsr(3));
    c.component(0).add_gen(*make_tgfsr(4));
    c.component(1).share_pool_with(c.component(0));

    ASSERT_TRUE(c.reset());
    int count = 1;
    while (c.next()) ++count;
    // C(4, 2) = 6 combos. No (i, i) self-pairs, and (i, j) with j > i only.
    EXPECT_EQ(count, 6);
}

TEST(Combination, KgIsSumAndLIsMinCappedAtLmax) {
    Combination c(2, 16);  // Lmax intentionally below the gens' L=32
    c.component(0).add_gen(*make_tgfsr(1));
    c.component(1).add_gen(*make_tgfsr(2));

    ASSERT_TRUE(c.reset());
    EXPECT_EQ(c.k_g(), c.at(0).k() + c.at(1).k());
    EXPECT_EQ(c.L(), 16);  // min(32, 32) = 32, capped at Lmax=16
}

TEST(Combination, AtReturnsActiveGeneratorAfterAdvance) {
    Combination c(1, 32);
    c.component(0).add_gen(*make_tgfsr(1));
    c.component(0).add_gen(*make_tgfsr(2));

    ASSERT_TRUE(c.reset());
    Generator* first = &c.at(0);
    ASSERT_TRUE(c.next());
    Generator* second = &c.at(0);
    EXPECT_NE(first, second);
}

TEST(Combination, NextAfterExhaustionStaysFalse) {
    Combination c(1, 32);
    c.component(0).add_gen(*make_tgfsr(1));

    ASSERT_TRUE(c.reset());
    EXPECT_FALSE(c.next());
    EXPECT_TRUE(c.exhausted());
    EXPECT_FALSE(c.next());  // still false on a second call
}

TEST(Combination, ResetAfterExhaustionReinitializes) {
    Combination c(1, 32);
    c.component(0).add_gen(*make_tgfsr(1));
    c.component(0).add_gen(*make_tgfsr(2));

    ASSERT_TRUE(c.reset());
    EXPECT_TRUE(c.next());
    EXPECT_FALSE(c.next());
    EXPECT_TRUE(c.exhausted());

    // Reset re-arms.
    ASSERT_TRUE(c.reset());
    EXPECT_FALSE(c.exhausted());
    EXPECT_TRUE(c.next());
}

TEST(Combination, JEquals3MixedSharedAndIndependent) {
    Combination c(3, 32);
    // pool A: 3 gens, slots 0 and 2 share it.
    c.component(0).add_gen(*make_tgfsr(1));
    c.component(0).add_gen(*make_tgfsr(2));
    c.component(0).add_gen(*make_tgfsr(3));
    // slot 1: independent pool of 2 gens.
    c.component(1).add_gen(*make_tgfsr(11));
    c.component(1).add_gen(*make_tgfsr(12));
    // slot 2: shares pool with slot 0.
    c.component(2).share_pool_with(c.component(0));

    ASSERT_TRUE(c.reset());
    int count = 1;
    while (c.next()) ++count;
    // C(3, 2) * 2 = 3 * 2 = 6 combos.
    EXPECT_EQ(count, 6);
}
