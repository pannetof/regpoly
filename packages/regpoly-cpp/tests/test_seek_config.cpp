// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

// Phase 4.2: tests for the seek-style YAML config loader.
//
// Covers parse-shape verification (the SeekConfig struct populated as
// expected for two of the canonical sample configs) and the build_search
// side that translates a SeekConfig into a Combination via the existing
// factories.

#include <gtest/gtest.h>

#include "combination.h"
#include "params.h"
#include "search_types.h"
#include "seek_config.h"
#include "seek_search.h"

#include <climits>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

namespace fs = std::filesystem;
using regpoly_yaml_config::ComponentSpec;
using regpoly_yaml_config::SeekConfig;
using regpoly_yaml_config::build_search;
using regpoly_yaml_config::load_seek_config;

namespace {

// Walk up from this source file's directory to find the monorepo root
// (the directory containing both `shared/` and `docs/`).
std::string repo_root() {
    fs::path here = fs::path(__FILE__).parent_path();
    for (int i = 0; i < 8; ++i) {
        if (fs::is_directory(here / "shared")
            && fs::is_directory(here / "docs")) {
            return here.string();
        }
        if (here == here.parent_path()) break;
        here = here.parent_path();
    }
    // Fall back to relative — if this fails the test will report it.
    return ".";
}

}  // namespace

TEST(SeekConfig, RejectsLegacyFileSourceWithMigrationPointer) {
    // The pre-v2 `legacy_file:` component source is rejected at the
    // YAML loader. The error message must point at regpoly-legacy.
    auto tmpdir = fs::temp_directory_path()
        / ("regpoly_seek_cfg_legacy_reject_" + std::to_string(::getpid()));
    fs::create_directories(tmpdir);
    auto yaml_path = tmpdir / "cfg.yaml";
    {
        std::ofstream f(yaml_path);
        f << "search:\n  Lmax: 32\ncomponents:\n"
          << "  - generators:\n      legacy_file: anything.dat\n"
          << "tests:\n  - type: equidistribution\n    max_gap_sum: 100\n";
    }
    try {
        load_seek_config(yaml_path.string());
        FAIL() << "expected load_seek_config to throw";
    } catch (const std::runtime_error& exc) {
        const std::string msg = exc.what();
        EXPECT_NE(msg.find("regpoly-legacy"), std::string::npos)
            << "error message should reference regpoly-legacy add-on; got: "
            << msg;
    }
    fs::remove_all(tmpdir);
}

TEST(SeekConfig, ParsesCombinedTwoComponents) {
    std::string path = repo_root() + "/shared/yaml/equidist/combined.yaml";
    SeekConfig cfg = load_seek_config(path);

    EXPECT_EQ(cfg.Lmax, 32);
    EXPECT_EQ(cfg.nbtries, 1);
    ASSERT_EQ(cfg.components.size(), 2u);

    // First component: inline TGFSRGen with tempering (tempMK with
    // random b/c — those should be silently dropped from params).
    const auto& c0 = cfg.components[0];
    EXPECT_EQ(c0.source, ComponentSpec::Source::Inline);
    EXPECT_EQ(c0.inline_family, "TGFSRGen");
    EXPECT_EQ(c0.common_params.get_int("w"), 32);
    EXPECT_EQ(c0.common_params.get_int("r"), 3);
    ASSERT_EQ(c0.per_gen_params.size(), 3u);
    // Each entry merged the common w/r and overlaid {a, m=1}.
    EXPECT_EQ(c0.per_gen_params[0].get_int("w"), 32);
    EXPECT_EQ(c0.per_gen_params[0].get_int("r"), 3);
    EXPECT_EQ(c0.per_gen_params[0].get_int("m"), 1);
    EXPECT_EQ(static_cast<uint64_t>(c0.per_gen_params[0].get_int("a")),
              0x8ebfd028u);

    ASSERT_EQ(c0.tempering.size(), 1u);
    EXPECT_EQ(c0.tempering[0].type, "tempMK");
    EXPECT_EQ(c0.tempering[0].params.get_int("w"), 32);
    EXPECT_EQ(c0.tempering[0].params.get_int("eta"), 7);
    EXPECT_EQ(c0.tempering[0].params.get_int("mu"), 15);
    // b/c are { random: ... } maps — must be ABSENT from params so the
    // factory falls back to its default randomization behavior.
    EXPECT_FALSE(c0.tempering[0].params.has("b"));
    EXPECT_FALSE(c0.tempering[0].params.has("c"));

    // Second component: inline, empty tempering.
    const auto& c1 = cfg.components[1];
    EXPECT_EQ(c1.source, ComponentSpec::Source::Inline);
    EXPECT_EQ(c1.inline_family, "TGFSRGen");
    EXPECT_EQ(c1.common_params.get_int("r"), 5);
    ASSERT_EQ(c1.per_gen_params.size(), 3u);
    EXPECT_EQ(c1.per_gen_params[0].get_int("m"), 2);
    EXPECT_TRUE(c1.tempering.empty());

    // Test: equidistribution + delta block.
    ASSERT_EQ(cfg.tests.size(), 1u);
    EXPECT_EQ(cfg.tests[0].kind, SeekTestKind::EquidistributionMatricial);
    EXPECT_EQ(cfg.tests[0].eq_mse, 10000);
    // delta rule: from=1, to=32, max=10000 — every entry [1..32] = 10000;
    // index 0 stays at INT_MAX.
    EXPECT_EQ(cfg.tests[0].eq_delta[0], INT_MAX);
    EXPECT_EQ(cfg.tests[0].eq_delta[1], 10000);
    EXPECT_EQ(cfg.tests[0].eq_delta[32], 10000);
}

TEST(SeekConfig, BuildSearchCombinedResetsCleanly) {
    std::string path = repo_root() + "/shared/yaml/equidist/combined.yaml";
    SeekConfig cfg = load_seek_config(path);
    auto built = build_search(cfg);
    ASSERT_NE(built.combination, nullptr);
    EXPECT_EQ(built.combination->J(), 2);
    EXPECT_FALSE(built.combination->exhausted());
    // Both components hold three TGFSRGens each.
    EXPECT_EQ(built.per_component_gens[0].size(), 3u);
    EXPECT_EQ(built.per_component_gens[1].size(), 3u);
    // First component has one tempering step.
    EXPECT_EQ(built.per_component_trans[0].size(), 1u);
    EXPECT_EQ(built.per_component_trans[1].size(), 0u);
}

TEST(SeekConfig, ParsesSameSentinelAsSharedPool) {
    std::string path = repo_root() + "/shared/yaml/equidist/same.yaml";
    SeekConfig cfg = load_seek_config(path);
    ASSERT_EQ(cfg.components.size(), 2u);
    EXPECT_EQ(cfg.components[0].source, ComponentSpec::Source::Inline);
    EXPECT_EQ(cfg.components[1].source, ComponentSpec::Source::Same);

    auto built = build_search(cfg);
    ASSERT_NE(built.combination, nullptr);
    EXPECT_EQ(built.combination->J(), 2);
    // Component 1 shares pool — pool_id should match component 0's.
    EXPECT_EQ(built.combination->component(0).pool_id(),
              built.combination->component(1).pool_id());
    EXPECT_FALSE(built.combination->exhausted());
}

TEST(SeekConfig, RejectsMissingFile) {
    EXPECT_THROW(load_seek_config("/no/such/path/missing.yaml"),
                 std::runtime_error);
}

TEST(SeekConfig, RejectsUnknownTestType) {
    auto tmpdir = fs::temp_directory_path()
        / ("regpoly_seek_cfg_test_unktest_" + std::to_string(::getpid()));
    fs::create_directories(tmpdir);
    auto yaml_path = tmpdir / "cfg.yaml";
    {
        std::ofstream f(yaml_path);
        f << "search:\n  Lmax: 32\ncomponents:\n"
          << "  - generators:\n      family: TGFSRGen\n"
          << "      common: { w: 32, r: 3 }\n"
          << "      generators:\n        - { a: 0x9908b0df, m: 1 }\n"
          << "tests:\n  - type: bogus_test\n";
    }
    EXPECT_THROW(load_seek_config(yaml_path.string()), std::runtime_error);
    fs::remove_all(tmpdir);
}

TEST(SeekConfig, BuiltSearchRunsSeekDriver) {
    // Smoke test: load combined.yaml, build, and run a few iterations.
    // Verifies the wiring end-to-end.
    std::string path = repo_root() + "/shared/yaml/equidist/combined.yaml";
    SeekConfig cfg = load_seek_config(path);
    auto built = build_search(cfg);

    int64_t iters = 0;
    int64_t selections = 0;
    auto on_iter = [&](Combination&, const SeekIterResult&) { ++selections; };
    auto on_progress = [&](const SeekProgress& p) { iters = p.nbgen; };

    auto result = run_seek_search(*built.combination, cfg.tests,
                                  /*nbtries=*/1,
                                  /*progress_interval=*/1000,
                                  /*on_prep=*/nullptr,
                                  on_iter, on_progress);
    EXPECT_GT(result.nbgen, 0);
    // Selections may be zero for a tight delta, just verify the loop
    // terminated.
    EXPECT_LE(selections, result.nbgen);
    (void)iters;
}
