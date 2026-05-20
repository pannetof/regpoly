// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.
//
// Phase 4: t-value kernel — Schmid-style primal enumeration.
//
// Locked properties:
//   1) s = 1 (skipped by definition; profile is over s ≥ 2). tvals[s=1] = 0.
//   2) Trivial net (C_j = I for all j) → t(s) = m·(s-1) for s ≥ 2:
//      with all matrices equal, every composition (l_1,…,l_s) yields a
//      stacked matrix whose rank is min(l_max, m). Full rank d
//      requires d ≤ m, achieved by composition (m, 0, …, 0); but for
//      d > m the stacked rank can't exceed m. So d_max = m and
//      t(s) = 0?  Let me reconsider — the trivial-net claim in the
//      plan needs more care. See test commentary below.
//   3) Hand-derived Sobol(m=2, s=2) and Niederreiter(m=2, s=2) goldens.
//   4) Rejection short-circuit on delta[s] and max_t_sum.

#include <gtest/gtest.h>

#include "digital_net.h"
#include "sobol.h"
#include "niederreiter_f2.h"
#include "tvalue_runner.h"
#include "bitvect.h"

#include <climits>
#include <memory>
#include <vector>

using regpoly::core::BitVect;
using regpoly::core::DigitalNet;
using regpoly::core::Generator;
using regpoly::core::NiederreiterF2Gen;
using regpoly::core::SobolNet;
using regpoly::core::run_tvalue_profile_schmid;
using regpoly::core::run_tvalue_profile_dual;
using regpoly::core::TValueResult;

namespace {

// Reuse the IdentityNet from test_digital_net.cpp's pattern. Defined
// locally here too to keep the test file self-contained.
class IdentityNet : public DigitalNet {
public:
    IdentityNet(int m, int s_max) : DigitalNet(m, s_max) {
        std::vector<BitVect> rows;
        rows.reserve(m);
        for (int r = 0; r < m; ++r) {
            BitVect row(m);
            row.set_bit(r, 1);
            rows.push_back(std::move(row));
        }
        identity_rows_ = std::make_shared<const std::vector<BitVect>>(std::move(rows));
    }
    std::string name() const override { return "IdentityNet"; }
    std::string display_str() const override {
        return "IdentityNet(m=" + std::to_string(m()) +
               ", s_max=" + std::to_string(s_max()) + ")";
    }
    std::unique_ptr<Generator> copy() const override {
        return std::unique_ptr<Generator>(new IdentityNet(*this));
    }
protected:
    const std::vector<BitVect>& generating_matrix(int) const override {
        return *identity_rows_;
    }
private:
    std::shared_ptr<const std::vector<BitVect>> identity_rows_;
};

std::vector<int> uncapped_delta(int s_max) {
    return std::vector<int>(s_max + 1, INT_MAX);
}

}  // namespace

// ── s=1 boundary ──────────────────────────────────────────────────────

TEST(TValueRunner, SMaxOneReturnsTrivialProfile) {
    IdentityNet net(2, 2);
    TValueResult res = run_tvalue_profile_schmid(
        net, /*kg=*/2, /*m=*/2, /*s_max=*/1, uncapped_delta(2), INT_MAX);
    EXPECT_TRUE(res.verified);
    EXPECT_EQ(res.se, 0);
}

// ── Trivial net (C_j = I for all j) ───────────────────────────────────
//
// All matrices equal → for any composition (l_1, ..., l_s), the
// stacked matrix has rank min(max(l_j), m) = max(l_j) (since each C_j
// is rank-m and identical, adding rows from a different j contributes
// nothing new). For sum(l_j) = d ≤ m, choose (d, 0, …, 0): max=d,
// rank=d, full rank. So d_max = m and t(s) = 0 for all s ≥ 2.
// (The plan sketch — "t(s) = m·(s-1)" — was a misremembering; the
// correct degenerate behaviour for a constant net is t(s) = 0 at
// every s where d ≤ m is achievable.)

TEST(TValueRunner, IdentityNetHasZeroTValueAtEveryS) {
    const int m = 3;
    IdentityNet net(m, m);
    TValueResult res = run_tvalue_profile_schmid(
        net, /*kg=*/m, /*m=*/m, /*s_max=*/m, uncapped_delta(m), INT_MAX);
    EXPECT_TRUE(res.verified);
    for (int s = 2; s <= m; ++s) {
        EXPECT_EQ(res.tvals[s], 0) << "s=" << s;
    }
    EXPECT_EQ(res.se, 0);
}

// ── Hand-derived Sobol(m=2, s=2) ──────────────────────────────────────
//
// Sobol with m=2: C_1 = [[1,0],[0,1]], C_2 = [[1,1],[0,1]].
// For s=2, d=2: try every composition (l_1, l_2) of 2:
//   (2, 0): top-2 rows of C_1 = full C_1, rank 2. Full rank!
//   So d_max(2) = 2, t(2) = 2 - 2 = 0.
// Sobol is a (0,m,2)-net for low m — the classical result.

TEST(TValueRunner, SobolM2S2HasTValueZero) {
    SobolNet net(/*m=*/2, /*s_max=*/2);
    auto res = run_tvalue_profile_schmid(
        net, /*kg=*/2, /*m=*/2, /*s_max=*/2, uncapped_delta(2), INT_MAX);
    EXPECT_TRUE(res.verified);
    EXPECT_EQ(res.tvals[2], 0);
    EXPECT_EQ(res.se, 0);
}

// ── Hand-derived Sobol(m=3, s=2) ──────────────────────────────────────
//
// Sobol with m=3, s_max=2.
//   C_1 = I_3 = [[1,0,0],[0,1,0],[0,0,1]]
//   C_2 = [[1,1,1],[0,1,1],[0,0,1]]  (Antonov–Saleev recurrence on
//                                     p=x+1 with all 1s direction
//                                     numbers — upper-triangular all-1s.)
//
// For s=2, d=3: try (3, 0) → top-3 of C_1 = I, rank 3. Full rank.
// So d_max = 3, t(2) = 0.

TEST(TValueRunner, SobolM3S2HasTValueZero) {
    SobolNet net(/*m=*/3, /*s_max=*/3);  // s_max=3 because m=3 requires s_max>=m
    auto res = run_tvalue_profile_schmid(
        net, /*kg=*/3, /*m=*/3, /*s_max=*/2, uncapped_delta(2), INT_MAX);
    EXPECT_TRUE(res.verified);
    EXPECT_EQ(res.tvals[2], 0);
}

// ── Hand-derived Sobol(m=3, s=3) ──────────────────────────────────────
//
// Sobol with m=3, s=3.
//   C_3 comes from Joe-Kuo j=3: deg=2, a=1, m_init=[1,3].
//   Recurrence: m_{3,k} = 2·a_1·m_{3,k-1} XOR 2^2·m_{3,k-2} XOR m_{3,k-2}
//                for k > 2, with a_1 derived from bit (deg-1-i) of a=1 → a_1 = 1.
//   m_{3,1}=1, m_{3,2}=3.
//   m_{3,3} = 2·1·3 XOR 4·1 XOR 1 = 6 XOR 4 XOR 1 = 3.
//   v_{3,k} = m_{3,k} << (3-k): v_{3,1}=4 (100), v_{3,2}=6 (110), v_{3,3}=3 (011).
//   Hmm 3 doesn't fit in 3 bits MSB-first... wait it does: 011 in 3 bits = 3.
//   Bit r of v: bit r of (v_{3,k}) read MSB-first in 3 bits.
//   v_{3,1}=4=100: row r=0 has 1, row r=1 has 0, row r=2 has 0. Col 0 = (1,0,0)^T.
//   v_{3,2}=6=110: bits 1,1,0. Col 1 = (1,1,0)^T.
//   v_{3,3}=3=011: bits 0,1,1. Col 2 = (0,1,1)^T.
//   So C_3 = [[1,1,0],[0,1,1],[0,0,1]].
//
// For s=3, d=3: try (1, 1, 1) → top-1 row of each:
//   C_1[0:1] = (1,0,0). C_2[0:1] = (1,1,1). C_3[0:1] = (1,1,0). Stacked
//   matrix is 3×3 with rows (1,0,0), (1,1,1), (1,1,0). Rank: subtract
//   row 0 from row 1: (0,1,1). subtract row 0 from row 2: (0,1,0).
//   Subtract new row 1 from new row 2: (0,0,1). Full rank 3.
// So d_max=3 and t(3) = 3 - 3 = 0.

TEST(TValueRunner, SobolM3S3HasTValueZero) {
    SobolNet net(/*m=*/3, /*s_max=*/3);
    auto res = run_tvalue_profile_schmid(
        net, /*kg=*/3, /*m=*/3, /*s_max=*/3, uncapped_delta(3), INT_MAX);
    EXPECT_TRUE(res.verified);
    EXPECT_EQ(res.tvals[2], 0);
    EXPECT_EQ(res.tvals[3], 0);
    EXPECT_EQ(res.se, 0);
}

// ── Hand-derived Niederreiter(m=2, s=2) ───────────────────────────────
//
// Niederreiter with m=2: C_1 = I, C_2 = [[1,1],[0,1]] (upper-tri all-1s).
// Same structure as Sobol(m=2,s=2). d_max=2 → t(2)=0.

TEST(TValueRunner, NiederreiterM2S2HasTValueZero) {
    NiederreiterF2Gen net(/*m=*/2, /*s_max=*/2);
    auto res = run_tvalue_profile_schmid(
        net, /*kg=*/2, /*m=*/2, /*s_max=*/2, uncapped_delta(2), INT_MAX);
    EXPECT_TRUE(res.verified);
    EXPECT_EQ(res.tvals[2], 0);
}

// ── Rejection short-circuit ───────────────────────────────────────────

TEST(TValueRunner, ShortCircuitsOnDeltaCap) {
    // IdentityNet has t(s) = 0 everywhere, so to provoke a short-circuit
    // we use a *negative-cap* delta (impossible to satisfy at any t).
    // Hmm, t ≥ 0 always; negative delta makes "tvals[s] > delta[s]"
    // trip when delta[s] = -1. Easier: use a net where t > 0 somewhere.
    //
    // For now use a tighter check: when delta[2] = -1, the first
    // s = 2 must fail because 0 > -1.
    IdentityNet net(2, 2);
    std::vector<int> delta(3, INT_MAX);
    delta[2] = -1;
    auto res = run_tvalue_profile_schmid(
        net, /*kg=*/2, /*m=*/2, /*s_max=*/2, delta, INT_MAX);
    EXPECT_FALSE(res.verified);
}

TEST(TValueRunner, ShortCircuitsOnMaxTSumCap) {
    IdentityNet net(2, 2);
    auto delta = uncapped_delta(2);
    // tvals[2] = 0 for IdentityNet → cumulative se = 0. To trip
    // max_t_sum we need se > max_t_sum at *some* point. 0 > -1 trips.
    auto res = run_tvalue_profile_schmid(
        net, /*kg=*/2, /*m=*/2, /*s_max=*/2, delta, /*max_t_sum=*/-1);
    EXPECT_FALSE(res.verified);
}

// ── Dual-method stub ──────────────────────────────────────────────────

TEST(TValueRunner, DualMethodStubThrows) {
    IdentityNet net(2, 2);
    EXPECT_THROW(run_tvalue_profile_dual(
                     net, 2, 2, 2, uncapped_delta(2), INT_MAX),
                 std::runtime_error);
}
