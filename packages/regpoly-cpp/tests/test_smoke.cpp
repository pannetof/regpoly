// Phase 0 smoke test: confirms the GoogleTest harness builds and links against
// regpoly_cpp_core. Phase 1+ replaces this with the CombinedGenerator suite and
// the parameterized "iterate this property over all 18 families" tables.

#include <gtest/gtest.h>

#include <regpoly/bitvect.h>

TEST(BitVectSmoke, ConstructionRecordsRequestedSize) {
    BitVect bv(8);
    EXPECT_EQ(bv.nbits(), 8);
}

TEST(BitVectSmoke, NewBitsAreZero) {
    BitVect bv(64);
    for (int i = 0; i < 64; ++i) {
        EXPECT_EQ(bv.get_bit(i), 0) << "bit " << i << " was not zero";
    }
}
