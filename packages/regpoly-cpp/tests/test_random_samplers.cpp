// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

// Phase 2.4 (TDD): generic random parameter samplers in C++.
//
// Mirrors the four generic rand_type cases of
// regpoly.core.parametric.generate_random. Each test seeds the
// thread-local RNG to a fixed value so failures are reproducible.

#include <gtest/gtest.h>

#include <algorithm>
#include <set>
#include <string>
#include <vector>

#include "param_spec.h"
#include "params.h"
#include "random_samplers.h"

using namespace regpoly::core;
using namespace regpoly::random;


namespace {

ParamSpec make_spec(const std::string& name, const std::string& rt,
                    const std::string& ra) {
    ParamSpec s;
    s.name = name;
    s.type = "int";
    s.structural = false;
    s.has_default = false;
    s.default_val = 0;
    s.rand_type = rt;
    s.rand_args = ra;
    s.optimizable = false;
    return s;
}

}  // namespace

TEST(EvalParamExpr, LiteralInt) {
    Params p;
    EXPECT_EQ(regpoly::random::eval_param_expr("42", p), 42);
    EXPECT_EQ(regpoly::random::eval_param_expr("-3", p), -3);
    EXPECT_EQ(regpoly::random::eval_param_expr("0", p), 0);
}

TEST(EvalParamExpr, LookupParam) {
    Params p;
    p.set_int("w", 32);
    p.set_int("k", 88);
    EXPECT_EQ(regpoly::random::eval_param_expr("w", p), 32);
    EXPECT_EQ(regpoly::random::eval_param_expr("k", p), 88);
    EXPECT_EQ(regpoly::random::eval_param_expr(" w ", p), 32);
}

TEST(EvalParamExpr, ParamPlusMinusN) {
    Params p;
    p.set_int("w", 32);
    EXPECT_EQ(regpoly::random::eval_param_expr("w+3", p), 35);
    EXPECT_EQ(regpoly::random::eval_param_expr("w-1", p), 31);
    EXPECT_EQ(regpoly::random::eval_param_expr("w + 3", p), 35);
}

TEST(EvalParamExpr, MissingParamThrows) {
    Params p;
    EXPECT_THROW(regpoly::random::eval_param_expr("bogus", p),
                 std::invalid_argument);
    EXPECT_THROW(regpoly::random::eval_param_expr("ghost+1", p),
                 std::invalid_argument);
}

TEST(SampleGeneric, Bitmask) {
    regpoly::random::seed(42);
    Params p;
    p.set_int("w", 8);
    auto spec = make_spec("a", "bitmask", "w");
    ASSERT_TRUE(regpoly::random::sample_generic_into(spec, p));
    int64_t v = p.get_int("a");
    EXPECT_GE(v, 0);
    EXPECT_LT(v, 1 << 8);
}

TEST(SampleGeneric, BitmaskWide) {
    regpoly::random::seed(42);
    Params p;
    p.set_int("w", 32);
    auto spec = make_spec("a", "bitmask", "w");
    bool any_high_bit = false;
    for (int i = 0; i < 200; ++i) {
        ASSERT_TRUE(regpoly::random::sample_generic_into(spec, p));
        int64_t v = p.get_int("a");
        EXPECT_GE(v, 0);
        EXPECT_LT(v, int64_t(1) << 32);
        if (v >= (int64_t(1) << 30)) any_high_bit = true;
    }
    EXPECT_TRUE(any_high_bit);
}

TEST(SampleGeneric, RangeStaysInBounds) {
    regpoly::random::seed(42);
    Params p;
    p.set_int("r", 10);
    auto spec = make_spec("m", "range", "1,r-1");
    int64_t lo_observed = INT64_MAX, hi_observed = INT64_MIN;
    for (int i = 0; i < 1000; ++i) {
        ASSERT_TRUE(regpoly::random::sample_generic_into(spec, p));
        int64_t v = p.get_int("m");
        ASSERT_GE(v, 1);
        ASSERT_LE(v, 9);
        lo_observed = std::min(lo_observed, v);
        hi_observed = std::max(hi_observed, v);
    }
    // Sanity: with 1000 draws over [1, 9] we should hit both endpoints.
    EXPECT_EQ(lo_observed, 1);
    EXPECT_EQ(hi_observed, 9);
}

TEST(SampleGeneric, RangeRejectsBadArgs) {
    Params p;
    p.set_int("r", 10);
    auto bad = make_spec("m", "range", "5");
    EXPECT_THROW(regpoly::random::sample_generic_into(bad, p),
                 std::invalid_argument);

    auto inverted = make_spec("m", "range", "10,1");
    EXPECT_THROW(regpoly::random::sample_generic_into(inverted, p),
                 std::invalid_argument);
}

TEST(SampleGeneric, PolyExponentsShape) {
    regpoly::random::seed(7);
    Params p;
    p.set_int("k", 32);
    auto spec = make_spec("nocoeff", "poly_exponents", "k");
    ASSERT_TRUE(regpoly::random::sample_generic_into(spec, p));
    auto v = p.get_int_vec("nocoeff");
    ASSERT_GE(v.size(), 2u);
    // Final element is the appended 0.
    EXPECT_EQ(v.back(), 0);
    // The leading n elements are sorted, distinct, in [1, k-1].
    for (size_t i = 0; i + 1 < v.size(); ++i) {
        EXPECT_GE(v[i], 1);
        EXPECT_LT(v[i], 32);
        if (i + 2 < v.size())
            EXPECT_LT(v[i], v[i + 1]);
    }
    // n is in [1, min(k-1, 10)].
    size_t n = v.size() - 1;
    EXPECT_GE(n, 1u);
    EXPECT_LE(n, 10u);
}

TEST(SampleGeneric, PolyExponentsRejectsSmallK) {
    Params p;
    p.set_int("k", 1);
    auto spec = make_spec("nocoeff", "poly_exponents", "k");
    EXPECT_THROW(regpoly::random::sample_generic_into(spec, p),
                 std::invalid_argument);
}

TEST(SampleGeneric, BitmaskVecMatchesLengthSourceShape) {
    regpoly::random::seed(11);
    Params p;
    p.set_int("w", 16);
    p.set_int_vec("nocoeff", {1, 5, 9, 0});  // length 4 → coeff vec length 4
    auto spec = make_spec("coeff", "bitmask_vec", "w,nocoeff");
    ASSERT_TRUE(regpoly::random::sample_generic_into(spec, p));
    auto v = p.get_uint_vec("coeff");
    ASSERT_EQ(v.size(), 4u);
    for (auto x : v) {
        EXPECT_LT(x, uint64_t(1) << 16);
    }
}

TEST(SampleGeneric, UnknownRandTypeReturnsFalse) {
    Params p;
    auto spec = make_spec("a", "tausworthe_poly", "k");
    EXPECT_FALSE(regpoly::random::sample_generic_into(spec, p));
}

TEST(SampleGeneric, NoneRandTypeReturnsFalse) {
    Params p;
    auto spec = make_spec("a", "none", "");
    EXPECT_FALSE(regpoly::random::sample_generic_into(spec, p));
    auto empty = make_spec("a", "", "");
    EXPECT_FALSE(regpoly::random::sample_generic_into(empty, p));
}

TEST(Seed, ProducesReproducibleSequence) {
    regpoly::random::seed(123);
    Params p;
    p.set_int("w", 16);
    auto spec = make_spec("a", "bitmask", "w");
    std::vector<int64_t> seq1;
    for (int i = 0; i < 10; ++i) {
        ASSERT_TRUE(regpoly::random::sample_generic_into(spec, p));
        seq1.push_back(p.get_int("a"));
    }
    regpoly::random::seed(123);
    std::vector<int64_t> seq2;
    for (int i = 0; i < 10; ++i) {
        ASSERT_TRUE(regpoly::random::sample_generic_into(spec, p));
        seq2.push_back(p.get_int("a"));
    }
    EXPECT_EQ(seq1, seq2);
}
