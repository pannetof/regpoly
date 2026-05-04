// Phase 3.3 (TDD): C++ port of regpoly.io.legacy_reader.
//
// Each test parses a real fixture from shared/legacy_parameters/,
// asserts the expected number of generators came back, and pokes at
// the first one to verify dispatch + Params construction made it
// through to a working create_generator() call.

#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <vector>

#include "legacy_reader.h"

namespace fs = std::filesystem;
using namespace regpoly_legacy;

namespace {

// Walk up from this source file (and from the CWD as a fallback) to
// find shared/legacy_parameters/. Same trick as test_catalog.cpp.
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

std::string fix(const std::string& name) {
    return (find_legacy_dir() / name).string();
}

}  // namespace

// ── PolyLCG ───────────────────────────────────────────────────────────

TEST(LegacyReader, PolyLCGSpecsFromFixture) {
    auto specs = read_generator_specs(fix("32poly.dat"), 32);
    EXPECT_EQ(specs.size(), 13u);
    EXPECT_EQ(specs.front().family, "PolyLCG");
    EXPECT_EQ(specs.front().params.get_int("k"), 32);
    auto poly = specs.front().params.get_int_vec("poly");
    ASSERT_FALSE(poly.empty());
    EXPECT_EQ(poly.back(), 0);   // sentinel preserved
}

TEST(LegacyReader, PolyLCGBuildsGenerators) {
    auto gens = read_generators(fix("32poly.dat"), 32);
    ASSERT_EQ(gens.size(), 13u);
    EXPECT_EQ(gens.front()->name(), "Polynomial LCG");
}

// ── Tausworthe ────────────────────────────────────────────────────────

TEST(LegacyReader, TauswortheSpecsFromTrinomials) {
    auto specs = read_generator_specs(fix("trinomials.dat"), 32);
    EXPECT_EQ(specs.size(), 9u);
    EXPECT_EQ(specs.front().family, "Tausworthe");
    // First record in `trinomials.dat` is `3 1 0` → reversed = {0,1,3},
    // k = Q.back() = 3.
    EXPECT_EQ(specs.front().params.get_int_vec("poly"),
              std::vector<int>({0, 1, 3}));
    EXPECT_TRUE(specs.front().params.get_bool("quicktaus"));
}

TEST(LegacyReader, TauswortheBuildsGenerators) {
    auto gens = read_generators(fix("trinomials.dat"), 32);
    ASSERT_EQ(gens.size(), 9u);
    EXPECT_EQ(gens.front()->name(), "TauswortheGen Generator");
}

// ── TGFSR ─────────────────────────────────────────────────────────────

TEST(LegacyReader, TGFSRSpecsFromFixture) {
    auto specs = read_generator_specs(fix("96_1.dt"), 32);
    EXPECT_EQ(specs.size(), 5u);
    EXPECT_EQ(specs.front().family, "TGFSR");
    EXPECT_EQ(specs.front().params.get_int("w"), 32);
    EXPECT_EQ(specs.front().params.get_int("r"), 3);
    // First row in 96_1.dt: "9965523b 1"
    EXPECT_EQ(static_cast<uint64_t>(
                 specs.front().params.get_int("a")),
              0x9965523bULL);
    EXPECT_EQ(specs.front().params.get_int("m"), 1);
}

TEST(LegacyReader, TGFSRBuildsGenerators) {
    auto gens = read_generators(fix("96_1.dt"), 32);
    ASSERT_EQ(gens.size(), 5u);
    EXPECT_EQ(gens.front()->name(), "TGFSRGen");
}

// ── Mersenne Twister ──────────────────────────────────────────────────

TEST(LegacyReader, MTSpecsFromMT19937Fixture) {
    auto specs = read_generator_specs(fix("MT19937.dat"), 32);
    ASSERT_EQ(specs.size(), 1u);
    EXPECT_EQ(specs.front().family, "MersenneTwister");
    EXPECT_EQ(specs.front().params.get_int("w"), 32);
    EXPECT_EQ(specs.front().params.get_int("r"), 624);
    EXPECT_EQ(specs.front().params.get_int("m"), 397);
    EXPECT_EQ(specs.front().params.get_int("p"), 31);
    EXPECT_EQ(static_cast<uint64_t>(
                 specs.front().params.get_int("a")),
              0x9908b0dfULL);
}

TEST(LegacyReader, MTBuildsGenerators) {
    auto gens = read_generators(fix("MT19937.dat"), 32);
    ASSERT_EQ(gens.size(), 1u);
    EXPECT_EQ(gens.front()->name(), "Mersenne Twister");
}

// ── GenF2w (PolyLCG variant) ──────────────────────────────────────────

TEST(LegacyReader, GenF2wSpecsFromFixture) {
    // genf2w_test.dat: gen_type=0 → F2wPolyLCGGen, normal_basis=false,
    // nbgen=2; per-gen: w=7 r=3 modM=0x49 step=1 nbcoeff=2
    //   coeffs: (0x3c, 2), (0x1c, 0)
    auto specs = read_generator_specs(fix("genf2w_test.dat"), 32);
    ASSERT_EQ(specs.size(), 2u);
    EXPECT_EQ(specs.front().family, "GenF2wPolyLCG");
    EXPECT_EQ(specs.front().params.get_int("w"), 7);
    EXPECT_EQ(specs.front().params.get_int("r"), 3);
    EXPECT_EQ(static_cast<uint64_t>(
                 specs.front().params.get_int("modM")),
              0x49ULL);
    EXPECT_EQ(specs.front().params.get_int("step"), 1);
    EXPECT_FALSE(specs.front().params.get_bool("normal_basis"));
    EXPECT_EQ(specs.front().params.get_string("type"), "polylcg");
    EXPECT_EQ(specs.front().params.get_uint_vec("coeff"),
              std::vector<uint64_t>({0x3c, 0x1c}));
    EXPECT_EQ(specs.front().params.get_int_vec("nocoeff"),
              std::vector<int>({2, 0}));
}

TEST(LegacyReader, GenF2wBuildsGenerators) {
    auto gens = read_generators(fix("genf2w_test.dat"), 32);
    ASSERT_EQ(gens.size(), 2u);
    // F2w generators inherit name() from F2wBaseGen.
    EXPECT_EQ(gens.front()->name(), "Generator in F_{2^w}");
}

// ── Carry (WELLGen) ───────────────────────────────────────────────────

TEST(LegacyReader, CarrySpecsFromFixture) {
    auto specs = read_generator_specs(fix("carry32_624_31_final.dat"), 32);
    ASSERT_EQ(specs.size(), 1u);
    EXPECT_EQ(specs.front().family, "WELLGen");
    EXPECT_EQ(specs.front().params.get_int("w"), 32);
    EXPECT_EQ(specs.front().params.get_int("r"), 624);
    EXPECT_EQ(specs.front().params.get_int("p"), 31);
    EXPECT_EQ(specs.front().params.get_int("m1"), 70);
    EXPECT_EQ(specs.front().params.get_int("m2"), 179);
    EXPECT_EQ(specs.front().params.get_int("m3"), 449);
    EXPECT_EQ(specs.front().params.get_int_vec("mat_types").size(), 8u);
    EXPECT_EQ(specs.front().params.get_int_vec("mat_pi").size(), 24u);
    EXPECT_EQ(specs.front().params.get_uint_vec("mat_pu").size(), 24u);
}

TEST(LegacyReader, CarryBuildsGenerators) {
    auto gens = read_generators(fix("carry32_624_31_final.dat"), 32);
    ASSERT_EQ(gens.size(), 1u);
    EXPECT_EQ(gens.front()->name(), "Carry Generator");
}

// ── Matsumoto ─────────────────────────────────────────────────────────

TEST(LegacyReader, MatsumotoSpecsFromFixture) {
    // matsumoto_test.dat: 2 gens.
    //   gen 1: type=1 n=3 m=1 nbpi=4 [5,-7,11,-15] nbpu=0
    //   gen 2: type=2 n=3 m=1 nbpi=4 [5,-7,11,-15] nbpu=1 [0x9908b0df]
    auto specs = read_generator_specs(fix("matsumoto_test.dat"), 32);
    ASSERT_EQ(specs.size(), 2u);
    EXPECT_EQ(specs.front().family, "MatsumotoGen");
    EXPECT_EQ(specs.front().params.get_int("type"), 1);
    EXPECT_EQ(specs.front().params.get_int("n"), 3);
    EXPECT_EQ(specs.front().params.get_int("m"), 1);
    EXPECT_EQ(specs.front().params.get_int_vec("paramsint"),
              std::vector<int>({5, -7, 11, -15}));
    EXPECT_TRUE(specs.front().params.get_uint_vec("paramsunsigned").empty());
    EXPECT_EQ(specs[1].params.get_int("type"), 2);
    EXPECT_EQ(specs[1].params.get_uint_vec("paramsunsigned"),
              std::vector<uint64_t>({0x9908b0dfULL}));
}

TEST(LegacyReader, MatsumotoBuildsGenerators) {
    auto gens = read_generators(fix("matsumoto_test.dat"), 32);
    ASSERT_EQ(gens.size(), 2u);
    EXPECT_EQ(gens.front()->name(), "Matsumoto Generator");
}

// ── MarsaXorshift ─────────────────────────────────────────────────────

TEST(LegacyReader, MarsaXorshiftSpecsExpandsTYPE1Variants) {
    // marsaxorshift_test.dat: header `1 1 32`, count `4`, then a
    // single TYPE1 row `1 13 -17 5`. TYPE1 produces 4 variants per row.
    auto specs = read_generator_specs(fix("marsaxorshift_test.dat"), 32);
    EXPECT_EQ(specs.size(), 4u);
    for (const auto& s : specs) {
        EXPECT_EQ(s.family, "MarsaXorshiftGen");
        EXPECT_EQ(s.params.get_int("type"), 1);
        EXPECT_EQ(s.params.get_int("w"), 32);
        EXPECT_EQ(s.params.get_int("r"), 1);
        EXPECT_EQ(s.params.get_int_vec("shifts").size(), 3u);
    }
    // Check the four variant tuples explicitly. Per the Python:
    //   (a,b,c), (c,b,a), (-a,-b,-c), (a,-c,-b)
    // with (a,b,c) = (13, -17, 5).
    EXPECT_EQ(specs[0].params.get_int_vec("shifts"),
              std::vector<int>({13, -17, 5}));
    EXPECT_EQ(specs[1].params.get_int_vec("shifts"),
              std::vector<int>({5, -17, 13}));
    EXPECT_EQ(specs[2].params.get_int_vec("shifts"),
              std::vector<int>({-13, 17, -5}));
    EXPECT_EQ(specs[3].params.get_int_vec("shifts"),
              std::vector<int>({13, -5, 17}));
}

TEST(LegacyReader, MarsaXorshiftBuildsGenerators) {
    auto gens = read_generators(fix("marsaxorshift_test.dat"), 32);
    ASSERT_EQ(gens.size(), 4u);
    EXPECT_EQ(gens.front()->name(), "Marsaglia Xor-shift");
}

// ── Transformations ───────────────────────────────────────────────────

TEST(LegacyReader, TransformationsSpecsFromTrans32) {
    // trans32.dat:
    //   2
    //   permut 32 -1 -1
    //   tempMKopt 32 7 15 -1 0 32
    auto result = read_transformation_specs(fix("trans32.dat"));
    ASSERT_EQ(result.specs.size(), 2u);
    EXPECT_TRUE(result.mk_opt) << "tempMKopt must set mk_opt=true";

    // permut: -1 sentinels are folded to 0.
    EXPECT_EQ(result.specs[0].trans_type, "permut");
    EXPECT_EQ(result.specs[0].params.get_int("w"), 32);
    EXPECT_EQ(result.specs[0].params.get_int("p"), 0);
    EXPECT_EQ(result.specs[0].params.get_int("q"), 0);

    // tempMKopt: w=32 eta=7 mu=15 nb_words=-1 (no b/c lines follow).
    EXPECT_EQ(result.specs[1].trans_type, "tempMK");
    EXPECT_EQ(result.specs[1].params.get_int("w"), 32);
    EXPECT_EQ(result.specs[1].params.get_int("eta"), 7);
    EXPECT_EQ(result.specs[1].params.get_int("mu"), 15);
    EXPECT_EQ(result.specs[1].params.get_int("u"), 0);
    EXPECT_EQ(result.specs[1].params.get_int("l"), 0);
    EXPECT_FALSE(result.specs[1].params.has("b"));
    EXPECT_FALSE(result.specs[1].params.has("c"));
}

TEST(LegacyReader, TransformationsBuiltFromTrans32) {
    auto built = read_transformations(fix("trans32.dat"));
    EXPECT_EQ(built.transformations.size(), 2u);
    EXPECT_TRUE(built.mk_opt);
    // Both transformations must be non-null.
    for (const auto& t : built.transformations) {
        ASSERT_NE(t.get(), nullptr);
    }
}

// ── Error handling ────────────────────────────────────────────────────

TEST(LegacyReader, UnknownTagThrows) {
    // Write a one-shot tempfile with a bogus tag and verify the throw.
    auto tmp = fs::temp_directory_path() / "regpoly_legacy_unknown.dat";
    {
        std::ofstream out(tmp);
        out << "notarealtag\n";
    }
    EXPECT_THROW(read_generator_specs(tmp.string(), 32), std::runtime_error);
    fs::remove(tmp);
}

TEST(LegacyReader, MissingFileThrows) {
    EXPECT_THROW(
        read_generator_specs("/nonexistent/path/legacy.dat", 32),
        std::runtime_error);
}
