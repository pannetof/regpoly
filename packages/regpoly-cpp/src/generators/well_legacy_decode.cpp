// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.

#include "well_legacy_decode.h"
#include <bitset>
#include <sstream>
#include <stdexcept>

namespace regpoly_well {

namespace {
constexpr uint32_t M32 = 0xFFFFFFFFu;

// Find the position of the single set bit in `mask32`, MSB-indexed
// (so bit at position 0 is the leftmost). Returns -1 if `mask32` does
// not have exactly one bit set.
int find_unique_msb_position(uint32_t mask32) {
    if (mask32 == 0) return -1;
    if (std::bitset<32>(mask32).count() != 1) return -1;
    // Lowest-set-bit index (LSB-indexed); convert to MSB-indexed.
    int lsb = 0;
    while (((mask32 >> lsb) & 1u) == 0) ++lsb;
    return 31 - lsb;
}
}  // namespace

int decode_d_s_mask(uint64_t mask, const std::string& context) {
    uint32_t m = static_cast<uint32_t>(mask & M32);
    // d_s = ~(1 << (31 - s))  → ~m has exactly one bit set, at MSB-position s.
    uint32_t inv = (~m) & M32;
    int s = find_unique_msb_position(inv);
    if (s < 0) {
        std::ostringstream oss;
        oss << context << ": cannot decode WELL paper M6 d_s mask 0x"
            << std::hex << m
            << " — expected exactly one zero bit (i.e. ~(1 << (31 - s))), "
            << "but got " << std::dec << std::bitset<32>(inv).count()
            << " zero bits. The legacy 3-mask M6 form falls outside the "
               "paper's parametric family and must be re-encoded by hand.";
        throw std::runtime_error(oss.str());
    }
    return s;
}

int decode_test_mask(uint64_t mask, const std::string& context) {
    uint32_t m = static_cast<uint32_t>(mask & M32);
    int t = find_unique_msb_position(m);
    if (t < 0) {
        std::ostringstream oss;
        oss << context << ": cannot decode WELL paper M6 test mask 0x"
            << std::hex << m
            << " — expected exactly one set bit (i.e. (1 << (31 - t))), "
            << "but got " << std::dec << std::bitset<32>(m).count()
            << " set bits. The legacy 3-mask M6 form falls outside the "
               "paper's parametric family and must be re-encoded by hand.";
        throw std::runtime_error(oss.str());
    }
    return t;
}

}  // namespace regpoly_well
