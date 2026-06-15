// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

// Phase 2.1 (TDD): cross-check Psi_12 / Phi_4 against the published
// Python algorithm. Reference vectors below were captured from the
// existing regpoly.analyses implementation on a small parameter grid.

#include <gtest/gtest.h>

#include <vector>

#include <regpoly/resolution_sets.h>

using namespace regpoly::core;


namespace {

// Convert a vector<bool> to the list of indices where true, for easier
// assertion messages.
std::vector<int> indices(const std::vector<bool>& v) {
    std::vector<int> out;
    for (int i = 0; i < static_cast<int>(v.size()); ++i)
        if (v[i]) out.push_back(i);
    return out;
}

}  // namespace

TEST(Psi12, MatchesPythonAlgorithmKg88L32) {
    // kg = 88 (3-component Tausworthe k=31+29+28), L=32.
    // Reference computed via the Python _compute_psi12: indices set to
    // True are 1..isqrt(88)=9, then for t=88/32=2..isqrt(87)=9 mark
    // min(88/t, 32) → t=2:32, t=3:29, t=4:22, t=5:17, t=6:14, t=7:12,
    // t=8:11, t=9:9.
    auto psi = compute_psi12(88, 32);
    EXPECT_EQ(static_cast<int>(psi.size()), 33);
    auto idx = indices(psi);
    std::vector<int> expected = {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                 11, 12, 14, 17, 22, 29, 32};
    EXPECT_EQ(idx, expected);
}

TEST(Psi12, MatchesPythonAlgorithmKg256L64) {
    // kg = 256, L = 64.
    // isqrt(256) = 16: marks 1..16.
    // m = max(2, 256/64) = 4. isqrt(255) = 15.
    // t=4: 64, t=5: 51, t=6: 42, t=7: 36, t=8: 32, t=9: 28, t=10: 25,
    // t=11: 23, t=12: 21, t=13: 19, t=14: 18, t=15: 17.
    auto psi = compute_psi12(256, 64);
    auto idx = indices(psi);
    std::vector<int> expected = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                                 17, 18, 19, 21, 23, 25, 28, 32, 36, 42, 51, 64};
    EXPECT_EQ(idx, expected);
}

TEST(Psi12, KgEquals1) {
    // Edge: kg = 1, L = 32.  isqrt(1) = 1 → marks 1.  m = max(2, 1/32)=2.
    // isqrt(0) = 0 → outer loop empty.  Result: only index 1.
    auto psi = compute_psi12(1, 32);
    EXPECT_EQ(indices(psi), std::vector<int>{1});
}

TEST(Phi4, MatchesPythonAlgorithmKg88L32) {
    // kg = 88, L = 32. Reference captured from Python _compute_phi4.
    auto phi = compute_phi4(88, 32);
    EXPECT_EQ(static_cast<int>(phi.size()), 89);
    std::vector<int> expected = {3, 5, 6, 7, 9, 10, 13, 15, 18, 30};
    EXPECT_EQ(indices(phi), expected);
}

TEST(Phi4, MatchesPythonAlgorithmKg256L64) {
    auto phi = compute_phi4(256, 64);
    std::vector<int> expected = {5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 18, 19,
                                 20, 22, 24, 26, 29, 37, 43, 52, 86};
    EXPECT_EQ(indices(phi), expected);
}
