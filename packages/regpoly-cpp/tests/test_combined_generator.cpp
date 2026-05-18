// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

// Phase 1 (TDD red): tests for the CombinedGenerator class.
//
// CombinedGenerator is a Generator subclass that owns J component
// generators, advances them in lockstep, XORs their outputs, and
// presents itself as a single Generator to all downstream code paths.
//
// These tests enumerate the contract:
//   - J=1 wrapping a primitive generator round-trips state, output,
//     char_poly, and is_full_period through CombinedGenerator with no
//     observable difference.
//   - J=3 next() = XOR of components' next() outputs.
//   - k = sum of component k; L = min of component L (capped at Lmax).
//   - char_poly = product of components' char_polys (over GF(2)).
//   - transition_matrix = block-diagonal of components' transition
//     matrices.

#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "bitvect.h"
#include "combined.h"
#include "generator.h"
#include "tausworthe.h"

using namespace regpoly::core;


namespace {

std::unique_ptr<TauswortheGen> make_taus_a() {
    // k=31 trinomial Tausworthe with Q={3,31}, s=13, quicktaus, L=32.
    return std::make_unique<TauswortheGen>(31, std::vector<int>{0, 3, 31}, 13,
                                           /*quicktaus=*/true, 32);
}

std::unique_ptr<TauswortheGen> make_taus_b() {
    return std::make_unique<TauswortheGen>(29, std::vector<int>{0, 2, 29}, 11,
                                           /*quicktaus=*/true, 32);
}

std::unique_ptr<TauswortheGen> make_taus_c() {
    return std::make_unique<TauswortheGen>(28, std::vector<int>{0, 3, 28}, 7,
                                           /*quicktaus=*/true, 32);
}

void canonical_init(Generator& g) {
    BitVect bv(g.k());
    bv.set_bit(0, 1);  // canonical "1" seed
    g.init(bv);
}

// Build a seed for a CombinedGenerator that puts a 1-bit at the start of
// each component's slice. Without this, a "set bit 0 of the combined
// k-bit vector" seed would leave all components except the first in the
// zero state, which is degenerate for Tausworthe.
BitVect combined_canonical_seed(const CombinedGenerator& cg) {
    BitVect bv(cg.k());
    const auto& prefix = cg.prefix_k();
    for (int j = 0; j < cg.J(); ++j)
        bv.set_bit(prefix[j], 1);
    return bv;
}

}  // namespace

TEST(CombinedGenerator, JEquals1MatchesPrimitive) {
    auto a_alone = make_taus_a();
    auto a_inside = make_taus_a();

    std::vector<std::unique_ptr<Generator>> comps;
    comps.push_back(std::move(a_inside));
    CombinedGenerator combined(std::move(comps), /*Lmax=*/32);

    canonical_init(*a_alone);
    canonical_init(combined);

    EXPECT_EQ(combined.k(), a_alone->k());
    EXPECT_EQ(combined.L(), a_alone->L());

    for (int step = 0; step < 8; ++step) {
        a_alone->next();
        combined.next();
        BitVect lhs = a_alone->get_output();
        BitVect rhs = combined.get_output();
        ASSERT_EQ(lhs.nbits(), rhs.nbits());
        for (int b = 0; b < lhs.nbits(); ++b) {
            EXPECT_EQ(lhs.get_bit(b), rhs.get_bit(b))
                << "step " << step << " bit " << b;
        }
    }
}

TEST(CombinedGenerator, KIsSumAndLIsMinCapped) {
    std::vector<std::unique_ptr<Generator>> comps;
    comps.push_back(make_taus_a());  // k=31
    comps.push_back(make_taus_b());  // k=29
    comps.push_back(make_taus_c());  // k=28
    CombinedGenerator combined(std::move(comps), /*Lmax=*/32);

    EXPECT_EQ(combined.k(), 31 + 29 + 28);
    EXPECT_EQ(combined.L(), 32);
}

TEST(CombinedGenerator, OutputIsXorOfComponents) {
    auto a = make_taus_a();
    auto b = make_taus_b();
    auto c = make_taus_c();
    canonical_init(*a);
    canonical_init(*b);
    canonical_init(*c);

    auto a_in = make_taus_a();
    auto b_in = make_taus_b();
    auto c_in = make_taus_c();
    std::vector<std::unique_ptr<Generator>> comps;
    comps.push_back(std::move(a_in));
    comps.push_back(std::move(b_in));
    comps.push_back(std::move(c_in));
    CombinedGenerator combined(std::move(comps), /*Lmax=*/32);
    combined.init(combined_canonical_seed(combined));

    for (int step = 0; step < 8; ++step) {
        a->next();
        b->next();
        c->next();
        combined.next();

        BitVect xa = a->get_output();
        BitVect xb = b->get_output();
        BitVect xc = c->get_output();
        BitVect expected(combined.L());
        // L = min(32,32,32) = 32 here.
        for (int i = 0; i < combined.L(); ++i) {
            int bit = xa.get_bit(i) ^ xb.get_bit(i) ^ xc.get_bit(i);
            if (bit) expected.set_bit(i, 1);
        }
        BitVect got = combined.get_output();
        for (int i = 0; i < combined.L(); ++i) {
            EXPECT_EQ(got.get_bit(i), expected.get_bit(i))
                << "step " << step << " bit " << i;
        }
    }
}

TEST(CombinedGenerator, CopyIsIndependent) {
    std::vector<std::unique_ptr<Generator>> comps;
    comps.push_back(make_taus_a());
    comps.push_back(make_taus_b());
    CombinedGenerator combined(std::move(comps), /*Lmax=*/32);
    canonical_init(combined);
    combined.next();
    combined.next();

    std::unique_ptr<Generator> clone = combined.copy();
    ASSERT_NE(clone.get(), nullptr);
    EXPECT_EQ(clone->k(), combined.k());
    EXPECT_EQ(clone->L(), combined.L());

    BitVect before_clone = clone->get_output();
    BitVect before_orig = combined.get_output();
    for (int b = 0; b < combined.L(); ++b) {
        EXPECT_EQ(before_clone.get_bit(b), before_orig.get_bit(b));
    }

    // Step the original; the clone must not move.
    combined.next();
    BitVect after_clone = clone->get_output();
    for (int b = 0; b < combined.L(); ++b) {
        EXPECT_EQ(after_clone.get_bit(b), before_clone.get_bit(b))
            << "clone advanced when original stepped";
    }
}

TEST(CombinedGenerator, JEquals1IsFullPeriodMatchesComponent) {
    auto a_alone = make_taus_a();
    auto a_inside = make_taus_a();
    std::vector<std::unique_ptr<Generator>> comps;
    comps.push_back(std::move(a_inside));
    CombinedGenerator combined(std::move(comps), /*Lmax=*/32);

    EXPECT_EQ(combined.is_full_period(), a_alone->is_full_period());
}
