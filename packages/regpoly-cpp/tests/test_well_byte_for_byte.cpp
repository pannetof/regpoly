// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.

// Byte-for-byte WELL output stream test.
//
// Locks the output of WELL19937a (paper Table II row 1) against a
// frozen binary fixture. This is the only test that catches an
// accidental case-body swap in apply_matrix(): equidistribution-based
// tests cannot tell two primitive WELLs apart.
//
// The WELL parameters are constructed inline below (formerly read from
// shared/legacy_parameters/carry32_624_31_final.dat through the
// pre-v2 .dat reader, which moved to the regpoly-legacy add-on).

#include <gtest/gtest.h>

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <memory>
#include <vector>

#include "bitvect.h"
#include "factory.h"
#include "params.h"

namespace fs = std::filesystem;

namespace {

fs::path find_fixtures_dir() {
    // The fixtures live under tests/fixtures/ in the source tree.
    return fs::path(__FILE__).parent_path() / "fixtures";
}

// Inline WELL19937a constructor. Mirrors paper Table II row 1 — the
// same parameters the pre-v2 fixture
// `carry32_624_31_final.dat` (now in packages/regpoly-legacy/shared/
// legacy_parameters/) would have produced via the legacy `.dat` path.
std::unique_ptr<Generator> make_well19937a(int L) {
    Params p;
    p.set_int("w", 32);
    p.set_int("r", 624);
    p.set_int("p", 31);
    p.set_int("m1", 70);
    p.set_int("m2", 179);
    p.set_int("m3", 449);

    StructMap matrices;
    auto set_M = [&](const char* slot, int64_t M) {
        StructEntry e;
        e["M"] = M;
        matrices[slot] = std::move(e);
    };
    auto set_Mt = [&](const char* slot, int64_t M, int64_t t) {
        StructEntry e;
        e["M"] = M;
        e["t"] = t;
        matrices[slot] = std::move(e);
    };

    set_Mt("T0", 3, -25);
    set_Mt("T1", 3,  27);
    set_Mt("T2", 2,   9);
    set_Mt("T3", 3,   1);
    set_M ("T4", 1);
    set_Mt("T5", 3,  -9);
    set_Mt("T6", 3, -21);
    set_Mt("T7", 3,  21);

    p.set_struct_map("matrices", std::move(matrices));
    return create_generator("WELLGen", p, L);
}

}  // namespace

TEST(WELLByteForByte, WELL19937aMatchesFixture) {
    auto gen = make_well19937a(32);
    ASSERT_NE(gen, nullptr);
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
