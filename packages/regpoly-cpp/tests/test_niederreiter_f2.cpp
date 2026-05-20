// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.
//
// Phase 3: NiederreiterF2Gen — irreducible enumeration and Laurent
// matrix construction.
//
// The canonical Niederreiter F_2 first few irreducibles (in
// ascending degree, lex within degree by binary value):
//
//   p_1 = x         (deg 1, "10" = 2)
//   p_2 = x+1       (deg 1, "11" = 3)
//   p_3 = x^2+x+1   (deg 2, "111" = 7)
//   p_4 = x^3+x+1   (deg 3, "1011" = 11)
//   p_5 = x^3+x^2+1 (deg 3, "1101" = 13)
//   p_6 = x^4+x+1   (deg 4, "10011" = 19)
//   p_7 = x^4+x^3+1 (deg 4, "11001" = 25)
//   p_8 = x^4+x^3+x^2+x+1 (deg 4, "11111" = 31)
//
// Hand-derived (m=2, s=2) C_1, C_2:
//   p_1 = x, e_1 = 1. For i=1: Q=0, k=0, r=1. q = x.
//     Laurent of 1/x in F_2((x^-1)): 1/x = x^-1, so a_1 = 1, a_n = 0 for n != 1.
//     Row 1, col ℓ: a^{(1)}(1, 0, ℓ) = L[1][0+ℓ+1]. ℓ=0 → L[1][1]=1; ℓ=1 → L[1][2]=0.
//     Row 1 = [1, 0].
//   For i=2: Q=1, k=0, r=2. q = x^2.
//     Laurent of 1/x^2: a_2 = 1, a_n = 0 else.
//     Row 2, col ℓ: L[2][0+ℓ+1]. ℓ=0 → L[2][1]=0; ℓ=1 → L[2][2]=1.
//     Row 2 = [0, 1].
//   So C_1 = [[1,0],[0,1]] = identity (matches the digital-sequence
//   convention where the first Niederreiter coordinate is van der Corput-like).
//
//   p_2 = x+1, e_2 = 1. For i=1: Q=0, k=0, r=1. q = x+1.
//     Laurent of 1/(x+1): coefficient of x^-n for n ≥ 1.
//     Using the recurrence a_1 = 1, a_{n} = q_0 · a_{n-1} = 1 · a_{n-1} for n ≥ 2.
//     So a_n = 1 for all n ≥ 1.
//     Row 1, col ℓ: L[1][ℓ+1] = 1 for all ℓ.
//     Row 1 = [1, 1].
//   For i=2: Q=1, k=0, r=2. q = (x+1)^2 = x^2 + 1 (in F_2).
//     Laurent of 1/(x^2+1): a_2=1, a_3 = q_0·a_1 + q_1·a_2 ... wait, q_0=1, q_1=0.
//     a_{n} for n ≥ 3 (using recurrence a_{n'} = sum_{j=0..d-1} q_j a_{n'-d+j}):
//     a_3 = q_0 · a_3-2+0 + q_1 · a_3-2+1 = 1·a_1 + 0·a_2 = 0 + 0 = 0.
//     a_4 = q_0 · a_2 + q_1 · a_3 = 1·1 + 0·0 = 1.
//     Row 2, col ℓ: L[2][ℓ+1]. ℓ=0 → L[2][1]=0; ℓ=1 → L[2][2]=1.
//     Row 2 = [0, 1].
//
//   So C_2 = [[1,1],[0,1]] — upper triangular with all 1s above
//   the diagonal. (This matches the Niederreiter literature.)

#include <gtest/gtest.h>

#include "niederreiter_f2.h"
#include "bitvect.h"
#include "equidistribution_runner.h"

#include <climits>
#include <vector>

using regpoly::core::BitVect;
using regpoly::core::NiederreiterF2Gen;
using regpoly::core::run_matricial_equidistribution;

TEST(NiederreiterF2, ConstructorRejectsExcessiveM) {
    EXPECT_THROW(NiederreiterF2Gen(/*m=*/31, /*s_max=*/31), std::invalid_argument);
}

TEST(NiederreiterF2, ConstructorRejectsSMaxLessThanM) {
    EXPECT_THROW(NiederreiterF2Gen(/*m=*/4, /*s_max=*/3), std::invalid_argument);
}

TEST(NiederreiterF2, IrreducibleEnumerationMatchesCanonical) {
    NiederreiterF2Gen net(4, 8);
    EXPECT_EQ(net.irreducible(1), 0b10u);    // x
    EXPECT_EQ(net.irreducible_degree(1), 1);
    EXPECT_EQ(net.irreducible(2), 0b11u);    // x+1
    EXPECT_EQ(net.irreducible_degree(2), 1);
    EXPECT_EQ(net.irreducible(3), 0b111u);   // x^2+x+1
    EXPECT_EQ(net.irreducible_degree(3), 2);
    EXPECT_EQ(net.irreducible(4), 0b1011u);  // x^3+x+1
    EXPECT_EQ(net.irreducible_degree(4), 3);
    EXPECT_EQ(net.irreducible(5), 0b1101u);  // x^3+x^2+1
    EXPECT_EQ(net.irreducible_degree(5), 3);
    EXPECT_EQ(net.irreducible(6), 0b10011u); // x^4+x+1
    EXPECT_EQ(net.irreducible_degree(6), 4);
    EXPECT_EQ(net.irreducible(7), 0b11001u); // x^4+x^3+1
    EXPECT_EQ(net.irreducible_degree(7), 4);
    EXPECT_EQ(net.irreducible(8), 0b11111u); // x^4+x^3+x^2+x+1
    EXPECT_EQ(net.irreducible_degree(8), 4);
}

TEST(NiederreiterF2, M2S2FirstCoordinateIsIdentity) {
    // C_1 = [[1,0],[0,1]] per the hand-derivation above.
    NiederreiterF2Gen net(/*m=*/2, /*s_max=*/2);
    BitVect e0(2); e0.set_bit(0, 1);   // e_0 = (1, 0)
    BitVect e1(2); e1.set_bit(1, 1);   // e_1 = (0, 1)

    net.init(e0);
    BitVect out = net.get_output();    // C_1 · e_0 = column 0 of C_1
    EXPECT_EQ(out.get_bit(0), 1);
    EXPECT_EQ(out.get_bit(1), 0);

    net.init(e1);
    out = net.get_output();            // C_1 · e_1 = column 1 of C_1
    EXPECT_EQ(out.get_bit(0), 0);
    EXPECT_EQ(out.get_bit(1), 1);
}

TEST(NiederreiterF2, M2S2SecondCoordinateIsUpperTriangular) {
    // C_2 = [[1,1],[0,1]] per the hand-derivation above.
    NiederreiterF2Gen net(/*m=*/2, /*s_max=*/2);
    BitVect e0(2); e0.set_bit(0, 1);
    BitVect e1(2); e1.set_bit(1, 1);

    net.init(e0);
    net.next();                        // j = 2
    BitVect out = net.get_output();    // C_2 · e_0 = column 0 of C_2 = (1, 0)
    EXPECT_EQ(out.get_bit(0), 1);
    EXPECT_EQ(out.get_bit(1), 0);

    net.init(e1);
    net.next();
    out = net.get_output();            // C_2 · e_1 = column 1 of C_2 = (1, 1)
    EXPECT_EQ(out.get_bit(0), 1);
    EXPECT_EQ(out.get_bit(1), 1);
}

TEST(NiederreiterF2, MatricialEquidistributionRunsUnchanged) {
    // Architectural-claim test: existing equidistribution kernel
    // runs on a NiederreiterF2 net unchanged.
    const int m = 4;
    NiederreiterF2Gen net(m, /*s_max=*/m);
    std::vector<int> delta(m + 1, INT_MAX);
    auto res = run_matricial_equidistribution(
        net, /*kg=*/m, /*L=*/m, /*Lmax=*/m, delta, /*mse=*/INT_MAX);
    EXPECT_TRUE(res.verified);
    for (int l = 1; l <= m; ++l) {
        EXPECT_GE(res.ecart[l], 0) << "l=" << l;
    }
    // Niederreiter at (m=4, s=4): the construction is asymptotically
    // optimal; for low m we expect a small but nonzero t-value, and
    // the equidistribution test should run without error.
}

TEST(NiederreiterF2, CopyClonesSharedMatrices) {
    NiederreiterF2Gen net(4, 4);
    auto clone = net.copy();
    EXPECT_NE(clone.get(), &net);
    EXPECT_EQ(clone->name(), "NiederreiterF2Gen");
}
