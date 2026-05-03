// Phase 2.4 (TDD): primitivity testing in C++.
//
// Mirrors the semantics of the previous Python implementation in
// packages/regpoly/src/regpoly/search/primitivity.py. Reference
// expectations were captured from the Python algorithm + standard
// number-theory sources.

#include <gtest/gtest.h>

#include <NTL/ZZ.h>

#include <algorithm>
#include <memory>
#include <sstream>
#include <vector>

#include <regpoly/primitivity.h>

TEST(MersennePrimeExponent, KnownExponents) {
    EXPECT_TRUE(is_mersenne_prime_exponent(2));
    EXPECT_TRUE(is_mersenne_prime_exponent(3));
    EXPECT_TRUE(is_mersenne_prime_exponent(31));
    EXPECT_TRUE(is_mersenne_prime_exponent(11213));
    EXPECT_TRUE(is_mersenne_prime_exponent(19937));
    EXPECT_TRUE(is_mersenne_prime_exponent(44497));
    EXPECT_TRUE(is_mersenne_prime_exponent(216091));
}

TEST(MersennePrimeExponent, NonMersenneRejected) {
    EXPECT_FALSE(is_mersenne_prime_exponent(1));
    EXPECT_FALSE(is_mersenne_prime_exponent(4));
    EXPECT_FALSE(is_mersenne_prime_exponent(6));
    EXPECT_FALSE(is_mersenne_prime_exponent(11));
    EXPECT_FALSE(is_mersenne_prime_exponent(20));
    EXPECT_FALSE(is_mersenne_prime_exponent(100));
    EXPECT_FALSE(is_mersenne_prime_exponent(521 + 1));
}

TEST(GetPrimitiveFactorsForK, MersennePrime) {
    // For a Mersenne prime exponent, the single factor is 2^k - 1.
    auto facs = get_primitive_factors_for_k(31);
    ASSERT_TRUE(facs.has_value());
    ASSERT_EQ(facs->size(), 1u);
    NTL::ZZ expected = NTL::power(NTL::ZZ(2), 31) - 1;  // 2147483647
    std::ostringstream oss;
    oss << expected;
    EXPECT_EQ((*facs)[0], oss.str());
}

TEST(GetPrimitiveFactorsForK, CompositeAggregatesDivisors) {
    // k = 12 has divisors {1,2,3,4,6,12}. Phi_d(2) factors:
    //   d=2 -> {3}, d=3 -> {7}, d=4 -> {5}, d=6 -> {3}, d=12 -> {13}.
    // Union, sorted numerically: {3, 5, 7, 13}.
    auto facs = get_primitive_factors_for_k(12);
    ASSERT_TRUE(facs.has_value());
    std::vector<std::string> got = *facs;
    std::vector<std::string> expected = {"3", "5", "7", "13"};
    EXPECT_EQ(got, expected);
}

TEST(GetPrimitiveFactorsForK, K10) {
    // k=10 divisors {1,2,5,10}: Phi_2 -> {3}, Phi_5 -> {31}, Phi_10 -> {11}.
    auto facs = get_primitive_factors_for_k(10);
    ASSERT_TRUE(facs.has_value());
    std::vector<std::string> expected = {"3", "11", "31"};
    EXPECT_EQ(*facs, expected);
}

TEST(GetPrimitiveFactorsForK, MissingKReturnsNullopt) {
    // 2_000_000 is comfortably past anything we have factor data for.
    auto facs = get_primitive_factors_for_k(2'000'000);
    EXPECT_FALSE(facs.has_value());
}

TEST(IsFullPeriod, KnownPrimitivePolynomial) {
    // x^7 + x^6 + 1 is primitive over GF(2). char_poly bit i = coeff of x^i,
    // for i in [0, k-1]; the leading x^k is implicit. So char_poly bits 0
    // and 6 are set, k = 7.
    BitVect cp(7);
    cp.set_bit(0, 1);
    cp.set_bit(6, 1);
    EXPECT_TRUE(is_full_period(cp, 7));
}

TEST(IsFullPeriod, IrreducibleButNotPrimitive) {
    // x^4 + x^3 + x^2 + x + 1 is irreducible over GF(2) (it is Phi_5),
    // but its order divides 5, not 2^4 - 1 = 15. So it must NOT be
    // reported as full-period.
    BitVect cp(4);
    cp.set_bit(0, 1);
    cp.set_bit(1, 1);
    cp.set_bit(2, 1);
    cp.set_bit(3, 1);
    EXPECT_FALSE(is_full_period(cp, 4));
}

TEST(IsFullPeriod, ReducibleRejected) {
    // x^4 + x^2 = x^2 (x^2 + 1) is reducible -> not primitive.
    BitVect cp(4);
    cp.set_bit(2, 1);
    EXPECT_FALSE(is_full_period(cp, 4));
}

TEST(IsFullPeriod, ConstantTermZeroRejected) {
    // x^3 (degree 3 with no constant term) cannot be primitive: the
    // factorization includes x as a factor.
    BitVect cp(3);  // all zeros
    EXPECT_FALSE(is_full_period(cp, 3));
}

TEST(IsFullPeriod, PrimitiveDegree10ExercisesFactorTable) {
    // x^10 + x^3 + 1 is primitive over GF(2). k=10 is non-Mersenne so
    // this path must consult the embedded factor table.
    BitVect cp(10);
    cp.set_bit(0, 1);
    cp.set_bit(3, 1);
    EXPECT_TRUE(is_full_period(cp, 10));
}

TEST(IsFullPeriod, NonPrimitiveDegree10Rejected) {
    // x^10 + 1 = (x + 1) * (x^9 + x^8 + ... + 1): obviously reducible.
    BitVect cp(10);
    cp.set_bit(0, 1);
    EXPECT_FALSE(is_full_period(cp, 10));
}
