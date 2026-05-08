// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.

// Unit tests for decoding the legacy 3-mask form of WELL paper M6 to
// the paper's (s, t) integer parameters.
//
// The legacy form held d_s as ~(1 << (31 - s)) and the test mask as
// (1 << (31 - t)). For every published M6 instance (paper Table II,
// only WELL44497a uses M6) the masks decode losslessly. These tests
// pin the decoder behaviour at single-bit / single-zero-bit shape and
// at every other shape — which must throw.

#include <gtest/gtest.h>

#include <stdexcept>

#include "well_legacy_decode.h"

using regpoly_well::decode_d_s_mask;
using regpoly_well::decode_test_mask;

TEST(WellLegacyDecode, DSMaskAcceptsAllSingleZeroBitPatterns) {
    // s ∈ [0..31] corresponds to the zero bit at MSB-position s.
    for (int s = 0; s < 32; ++s) {
        uint32_t m = ~(1u << (31 - s));
        EXPECT_EQ(decode_d_s_mask(m, "test"), s)
            << "failed for s=" << s << " mask=0x" << std::hex << m;
    }
}

TEST(WellLegacyDecode, DSMaskRejectsAllOnes) {
    EXPECT_THROW(decode_d_s_mask(0xFFFFFFFFu, "test"), std::runtime_error);
}

TEST(WellLegacyDecode, DSMaskRejectsAllZeros) {
    EXPECT_THROW(decode_d_s_mask(0u, "test"), std::runtime_error);
}

TEST(WellLegacyDecode, DSMaskRejectsTwoZeroBits) {
    // Mask with two zero bits — outside the paper's family.
    EXPECT_THROW(decode_d_s_mask(0xFFFFFFF3u, "test"), std::runtime_error);
}

TEST(WellLegacyDecode, DSMaskMatchesPublishedWELL44497a) {
    // Paper Table II for WELL44497a has T6 = M6(9, 14, 5, a7).
    // The legacy d_s mask shipped in shared/yaml/.../WELLRNG/p15_r1391_w32.yaml
    // is 0xfbffffff (zero bit at MSB-position 5).
    EXPECT_EQ(decode_d_s_mask(0xfbffffffull, "WELL44497a"), 5);
}

TEST(WellLegacyDecode, TestMaskAcceptsAllSingleSetBitPatterns) {
    for (int t = 0; t < 32; ++t) {
        uint32_t m = 1u << (31 - t);
        EXPECT_EQ(decode_test_mask(m, "test"), t)
            << "failed for t=" << t << " mask=0x" << std::hex << m;
    }
}

TEST(WellLegacyDecode, TestMaskRejectsZero) {
    EXPECT_THROW(decode_test_mask(0u, "test"), std::runtime_error);
}

TEST(WellLegacyDecode, TestMaskRejectsTwoSetBits) {
    EXPECT_THROW(decode_test_mask(0x3u, "test"), std::runtime_error);
}

TEST(WellLegacyDecode, TestMaskMatchesPublishedWELL44497a) {
    // Legacy test mask 0x00020000 from WELL44497a → t=14.
    EXPECT_EQ(decode_test_mask(0x00020000ull, "WELL44497a"), 14);
}

TEST(WellLegacyDecode, IgnoresUpper32BitsInUint64Input) {
    // legacy_reader carries masks as uint64_t; only low 32 bits matter.
    uint64_t high_garbage = (static_cast<uint64_t>(0xdeadbeef) << 32)
                            | 0xfbffffffull;
    EXPECT_EQ(decode_d_s_mask(high_garbage, "test"), 5);
}
