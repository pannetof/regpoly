// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.

// Tests for CombinedF2LinearSource — the Recurrence-typed XOR wrapper
// that absorbs the role of the legacy CombinedGenerator.
//
// Contract:
//   - J=1 wrapping a primitive Recurrence round-trips state, output,
//     char_poly, and is_full_period through CombinedF2LinearSource with
//     no observable difference.
//   - J>=2 next() produces the XOR of components' next() outputs.
//   - k = sum of component k; L = min of component L (capped at Lmax).
//   - char_poly = product of components' char_polys (over GF(2)).
//   - copy()/clone_recurrence() deep-clone each component.
//   - default_test_method() unanimous-shortcut → component method;
//     disagreement or J>=2 with no unanimous → "notprimitive".
//   - flows through the lattice adapter (test_me_lat) without
//     dynamic_cast failures.

#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "bitvect.h"
#include "combined_f2_linear_source.h"
#include "equidistribution_method.h"
#include "f2_linear_source.h"
#include "generator.h"
#include "me_helpers.h"
#include "tausworthe.h"

using namespace regpoly::core;

namespace {

std::unique_ptr<TauswortheGen> make_taus_a() {
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

void canonical_init(ITestable& g) {
    BitVect bv(g.k());
    bv.set_bit(0, 1);
    g.init(bv);
}

// Seed for a combined J-component generator that puts a 1-bit at the
// start of each component's slice — otherwise components past the
// first land in the zero state, which is degenerate for Tausworthe.
BitVect combined_canonical_seed(const CombinedF2LinearSource& cg) {
    BitVect bv(cg.k());
    const auto& prefix = cg.prefix_k();
    for (int j = 0; j < cg.J(); ++j)
        bv.set_bit(prefix[j], 1);
    return bv;
}

}  // namespace

TEST(CombinedF2LinearSource, JEquals1MatchesPrimitive) {
    auto a_alone = make_taus_a();
    auto a_inside = make_taus_a();

    std::vector<std::unique_ptr<F2LinearSource>> comps;
    comps.push_back(std::move(a_inside));
    CombinedF2LinearSource combined(std::move(comps), /*Lmax=*/32);

    canonical_init(*a_alone);
    canonical_init(combined);

    EXPECT_EQ(combined.k(), a_alone->k());
    EXPECT_EQ(combined.L(), a_alone->L());
    EXPECT_EQ(combined.J(), 1);

    for (int step = 0; step < 8; ++step) {
        a_alone->next();
        combined.next();
        BitVect lhs = a_alone->get_output();
        BitVect rhs = combined.get_output();
        ASSERT_EQ(lhs.nbits(), rhs.nbits());
        for (int b = 0; b < lhs.nbits(); ++b)
            EXPECT_EQ(lhs.get_bit(b), rhs.get_bit(b))
                << "step " << step << " bit " << b;
    }
}

TEST(CombinedF2LinearSource, KIsSumAndLIsMinCapped) {
    std::vector<std::unique_ptr<F2LinearSource>> comps;
    comps.push_back(make_taus_a());  // k=31
    comps.push_back(make_taus_b());  // k=29
    comps.push_back(make_taus_c());  // k=28
    CombinedF2LinearSource combined(std::move(comps), /*Lmax=*/32);

    EXPECT_EQ(combined.k(), 31 + 29 + 28);
    EXPECT_EQ(combined.L(), 32);
}

TEST(CombinedF2LinearSource, JEquals3IsXorOfComponents) {
    auto a = make_taus_a();
    auto b = make_taus_b();
    auto c = make_taus_c();
    canonical_init(*a);
    canonical_init(*b);
    canonical_init(*c);

    std::vector<std::unique_ptr<F2LinearSource>> comps;
    comps.push_back(make_taus_a());
    comps.push_back(make_taus_b());
    comps.push_back(make_taus_c());
    CombinedF2LinearSource combined(std::move(comps), /*Lmax=*/32);
    combined.init(combined_canonical_seed(combined));

    for (int step = 0; step < 8; ++step) {
        a->next();
        b->next();
        c->next();
        combined.next();

        BitVect xa = a->get_output();
        BitVect xb = b->get_output();
        BitVect xc = c->get_output();
        BitVect got = combined.get_output();
        for (int i = 0; i < combined.L(); ++i) {
            int expected_bit = xa.get_bit(i) ^ xb.get_bit(i) ^ xc.get_bit(i);
            EXPECT_EQ(got.get_bit(i), expected_bit)
                << "step " << step << " bit " << i;
        }
    }
}

TEST(CombinedF2LinearSource, CopyIsIndependent) {
    std::vector<std::unique_ptr<F2LinearSource>> comps;
    comps.push_back(make_taus_a());
    comps.push_back(make_taus_b());
    CombinedF2LinearSource combined(std::move(comps), /*Lmax=*/32);
    combined.init(combined_canonical_seed(combined));
    combined.next();
    combined.next();

    // `copy()` returns an ITestable deep-clone. Stepping the original
    // must NOT move the clone.
    std::unique_ptr<ITestable> clone = combined.copy();
    ASSERT_NE(clone.get(), nullptr);
    EXPECT_EQ(clone->k(), combined.k());
    EXPECT_EQ(clone->L(), combined.L());

    BitVect before_clone = clone->get_output();
    BitVect before_orig = combined.get_output();
    for (int b = 0; b < combined.L(); ++b)
        EXPECT_EQ(before_clone.get_bit(b), before_orig.get_bit(b));

    combined.next();
    BitVect after_clone = clone->get_output();
    for (int b = 0; b < combined.L(); ++b)
        EXPECT_EQ(after_clone.get_bit(b), before_clone.get_bit(b))
            << "clone advanced when original stepped";
}

TEST(CombinedF2LinearSource, ComponentsExposesUnderlyingRecurrences) {
    // Composition contract: `components()` returns the Recurrence
    // pointers that kernels walk for per-component state evolution.
    std::vector<std::unique_ptr<F2LinearSource>> comps;
    comps.push_back(make_taus_a());
    comps.push_back(make_taus_b());
    F2LinearSource* a_addr = comps[0].get();
    F2LinearSource* b_addr = comps[1].get();
    CombinedF2LinearSource combined(std::move(comps), 32);

    auto cs = combined.components();
    ASSERT_EQ(cs.size(), 2u);
    EXPECT_EQ(cs[0], a_addr);
    EXPECT_EQ(cs[1], b_addr);
}

TEST(CombinedF2LinearSource, SourcesReturnsComponentPointers) {
    auto a = make_taus_a();
    auto b = make_taus_b();

    std::vector<std::unique_ptr<F2LinearSource>> comps;
    comps.push_back(std::move(a));
    comps.push_back(std::move(b));
    F2LinearSource* a_addr = comps[0].get();
    F2LinearSource* b_addr = comps[1].get();
    CombinedF2LinearSource combined(std::move(comps), 32);

    auto srcs = combined.sources();
    ASSERT_EQ(srcs.size(), 2u);
    EXPECT_EQ(srcs[0], static_cast<const F2LinearSource*>(a_addr));
    EXPECT_EQ(srcs[1], static_cast<const F2LinearSource*>(b_addr));
}

TEST(CombinedF2LinearSource, DefaultTestMethodUnanimousShortcut) {
    // Both Tausworthe components return "matricial" individually, so
    // the unanimous-shortcut honors it.
    std::vector<std::unique_ptr<F2LinearSource>> comps;
    comps.push_back(make_taus_a());
    comps.push_back(make_taus_b());
    CombinedF2LinearSource combined(std::move(comps), 32);

    auto m = combined.default_test_method("equidistribution");
    ASSERT_TRUE(m.has_value());
    EXPECT_EQ(*m, "matricial");
}

TEST(CombinedF2LinearSource, DefaultTestMethodDelegatesForJEquals1) {
    auto a = make_taus_a();
    auto expected = a->default_test_method("equidistribution");

    std::vector<std::unique_ptr<F2LinearSource>> comps;
    comps.push_back(std::move(a));
    CombinedF2LinearSource combined(std::move(comps), 32);

    auto m = combined.default_test_method("equidistribution");
    EXPECT_EQ(m.has_value(), expected.has_value());
    if (m.has_value() && expected.has_value())
        EXPECT_EQ(*m, *expected);
}

TEST(CombinedF2LinearSource, EmptyComponentsThrows) {
    std::vector<std::unique_ptr<F2LinearSource>> empty;
    EXPECT_THROW(CombinedF2LinearSource(std::move(empty), 32),
                 std::invalid_argument);
}

TEST(CombinedF2LinearSource, FlowsThroughLatticeAdapter) {
    // Lattice runners take `const ITestable&`. A CombinedF2LinearSource
    // of all-Recurrence components must run through `test_me_lat`
    // without throwing — exercises the dynamic_cast unpack inside
    // `single_gen_adapters.cpp`.
    std::vector<std::unique_ptr<F2LinearSource>> comps;
    comps.push_back(make_taus_a());
    comps.push_back(make_taus_b());
    CombinedF2LinearSource combined(std::move(comps), 32);

    BitVect seed = combined_canonical_seed(combined);
    combined.init(seed);

    int kg = combined.k();
    int L  = combined.L();
    int maxL = 4;
    std::vector<int> delta(maxL + 1, 1000);
    EXPECT_NO_THROW({
        auto r = test_me_lat(combined, kg, L, maxL, delta, /*mse=*/1000);
        EXPECT_EQ(static_cast<int>(r.ecart.size()), maxL + 1);
    });
}
