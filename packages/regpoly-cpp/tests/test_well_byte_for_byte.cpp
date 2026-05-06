// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.

// Byte-for-byte WELL output stream test.
//
// Locks the output of WELL19937a (paper Table II, first row of
// shared/legacy_parameters/carry32_624_31_final.dat) against a frozen
// binary fixture. The fixture was captured before the M0..M6 rename;
// after the rename, the read-time remap in legacy_reader.cpp must
// produce a generator that emits the same bytes.
//
// This is the only test that catches an accidental case-body swap in
// apply_matrix(): equidistribution-based tests cannot tell two
// primitive WELLs apart.

#include <gtest/gtest.h>

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <vector>

#include "bitvect.h"
#include "legacy_reader.h"

namespace fs = std::filesystem;
using namespace regpoly_legacy;

namespace {

fs::path find_legacy_dir() {
    auto try_paths = std::vector<fs::path>{
        fs::path(__FILE__).parent_path(),
        fs::current_path(),
    };
    for (auto base : try_paths) {
        for (int up = 0; up < 8; ++up) {
            auto cand = base / "shared" / "legacy_parameters";
            if (fs::is_directory(cand)) return cand;
            if (base == base.parent_path()) break;
            base = base.parent_path();
        }
    }
    ADD_FAILURE() << "could not find shared/legacy_parameters directory";
    return {};
}

fs::path find_fixtures_dir() {
    // The fixtures live under tests/fixtures/ in the source tree.
    return fs::path(__FILE__).parent_path() / "fixtures";
}

}  // namespace

TEST(WELLByteForByte, WELL19937aMatchesFixture) {
    auto dat = (find_legacy_dir() / "carry32_624_31_final.dat").string();
    auto gens = read_generators(dat, 32);
    ASSERT_EQ(gens.size(), 1u);
    auto& gen = gens.front();
    ASSERT_EQ(gen->name(), "Carry Generator");

    // Deterministic seed: bit 0 set, all others zero.
    BitVect seed(gen->k());
    seed.set_bit(0, 1);
    gen->init(seed);

    // Pull 1024 outputs; pack each as little-endian uint32_t.
    std::vector<uint8_t> bytes;
    bytes.reserve(1024 * 4);
    for (int i = 0; i < 1024; i++) {
        gen->next();
        BitVect out = gen->get_output();
        uint32_t w = 0;
        for (int b = 0; b < 32 && b < out.nbits(); b++) {
            if (out.get_bit(b)) w |= (1u << b);
        }
        for (int by = 0; by < 4; by++)
            bytes.push_back(static_cast<uint8_t>((w >> (8 * by)) & 0xFF));
    }

    auto fixture = find_fixtures_dir() / "well19937a_1024.bin";

    // Capture mode: if the fixture is missing, write it and fail loudly.
    // Re-running the test will then verify against the captured bytes.
    if (!fs::exists(fixture)) {
        fs::create_directories(fixture.parent_path());
        std::ofstream out(fixture, std::ios::binary);
        out.write(reinterpret_cast<const char*>(bytes.data()),
                  static_cast<std::streamsize>(bytes.size()));
        FAIL() << "Captured fixture at " << fixture
               << " (" << bytes.size() << " bytes). "
               << "Re-run the test to verify the comparison path.";
    }

    std::ifstream in(fixture, std::ios::binary);
    std::vector<uint8_t> expected(
        (std::istreambuf_iterator<char>(in)),
        std::istreambuf_iterator<char>());
    ASSERT_EQ(bytes.size(), expected.size())
        << "Fixture size mismatch: fixture is " << expected.size()
        << " bytes, generator produced " << bytes.size();
    EXPECT_EQ(bytes, expected)
        << "WELL19937a byte stream diverged from fixture at "
        << fixture << ". If this divergence is intentional, regenerate "
        << "the fixture by deleting it and re-running the test.";
}
