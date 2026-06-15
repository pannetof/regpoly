// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.

// Phase 1 of the hierarchy refactor (additive). Round-trip tests for
// the new TemperingChain wrapper around vector<unique_ptr<Transformation>>.
// At this phase nothing else uses it; the chain is verified in isolation
// against a known TemperMK transformation.

#include <gtest/gtest.h>

#include <cstdint>
#include <memory>
#include <utility>

#include "bitvect.h"
#include "tempering_chain.h"
#include "temper_mk.h"

using namespace regpoly::core;


namespace {

// Build the standard MT19937 tempMK2 step (b=0x9D2C5680, c=0xEFC60000).
std::unique_ptr<Transformation> make_mt_temper() {
    return std::make_unique<TemperMKTrans>(
        /*w=*/32, /*type=*/2,
        /*eta=*/7, /*mu=*/15, /*u=*/11, /*l=*/18,
        /*b=*/0x9D2C5680ULL, /*c=*/0xEFC60000ULL);
}

// Build a 32-bit BitVect carrying `val` in word 0.
BitVect make_bv32(uint64_t val) {
    BitVect bv(32);
    bv.set_word(/*idx=*/0, /*w=*/32, val);
    return bv;
}

// Compare two 32-bit BitVects by their packed word value.
bool same32(const BitVect& a, const BitVect& b) {
    return a.get_word(0, 32) == b.get_word(0, 32);
}

}  // namespace


TEST(TemperingChain, EmptyIsIdentity) {
    TemperingChain chain;
    EXPECT_TRUE(chain.empty());
    EXPECT_EQ(chain.size(), 0);

    BitVect bv     = make_bv32(0xCAFEBABE);
    BitVect before = bv.copy();

    chain.apply(bv);

    EXPECT_TRUE(same32(bv, before));
}


TEST(TemperingChain, AppliesSingleStep) {
    TemperingChain chain;
    chain.add(make_mt_temper());
    EXPECT_EQ(chain.size(), 1);
    EXPECT_FALSE(chain.empty());

    BitVect via_chain = make_bv32(0xDEADBEEF);
    chain.apply(via_chain);

    BitVect direct = make_bv32(0xDEADBEEF);
    auto t = make_mt_temper();
    t->apply(direct);

    EXPECT_TRUE(same32(via_chain, direct));
}


TEST(TemperingChain, AppliesInOrder) {
    // Two consecutive copies of the same step should compose:
    // chain(x) == step(step(x)).
    TemperingChain chain;
    chain.add(make_mt_temper());
    chain.add(make_mt_temper());

    BitVect via_chain = make_bv32(0x12345678);
    chain.apply(via_chain);

    BitVect direct = make_bv32(0x12345678);
    auto t1 = make_mt_temper();
    auto t2 = make_mt_temper();
    t1->apply(direct);
    t2->apply(direct);

    EXPECT_TRUE(same32(via_chain, direct));
}


TEST(TemperingChain, CopyCtorIsDeep) {
    TemperingChain original;
    original.add(make_mt_temper());

    TemperingChain cloned = original;  // copy ctor
    EXPECT_EQ(cloned.size(), 1);

    // The two chains' first steps should be DIFFERENT C++ objects (deep
    // clone), even though they apply identically.
    EXPECT_NE(&original.at(0), &cloned.at(0));

    // And they should produce identical output.
    BitVect a = make_bv32(0xABAD1DEA);
    BitVect b = make_bv32(0xABAD1DEA);
    original.apply(a);
    cloned.apply(b);
    EXPECT_TRUE(same32(a, b));
}


TEST(TemperingChain, CopyAssignmentReplacesAndDeepClones) {
    TemperingChain original;
    original.add(make_mt_temper());

    TemperingChain target;
    target.add(make_mt_temper());
    target.add(make_mt_temper());
    EXPECT_EQ(target.size(), 2);

    target = original;  // copy assignment replaces target's contents
    EXPECT_EQ(target.size(), 1);
    EXPECT_NE(&original.at(0), &target.at(0));

    BitVect a = make_bv32(0xC0FFEE00);
    BitVect b = make_bv32(0xC0FFEE00);
    original.apply(a);
    target.apply(b);
    EXPECT_TRUE(same32(a, b));
}


TEST(TemperingChain, MoveConstructTransfersOwnership) {
    TemperingChain src;
    src.add(make_mt_temper());
    src.add(make_mt_temper());

    Transformation* first_step_addr = &src.at(0);

    TemperingChain dst = std::move(src);
    EXPECT_EQ(dst.size(), 2);
    EXPECT_EQ(&dst.at(0), first_step_addr);  // same heap object, moved
    // After move-from, src is in a valid-but-unspecified state; size
    // is implementation-defined. Don't assert anything about src.
}


TEST(TemperingChain, ConstructFromVector) {
    std::vector<std::unique_ptr<Transformation>> steps;
    steps.push_back(make_mt_temper());
    steps.push_back(make_mt_temper());

    TemperingChain chain(std::move(steps));
    EXPECT_EQ(chain.size(), 2);
}
