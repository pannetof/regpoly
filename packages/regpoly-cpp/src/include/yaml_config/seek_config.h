// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

// Phase 4.2: focused YAML loader for the `regpoly-cli search` driver.
//
// Mirrors the subset of regpoly.search.seek.Seek.from_yaml needed by
// the C++ CLI: parse a seek-style search config (search/components/
// tests), build the corresponding Combination + tempering chains via
// the existing factories, and hand off to run_seek_search().
//
// The schema supported here is a strict subset of the Python loader's
// — two component sources (inline family/common/generators and the
// literal "same" sentinel), simple tempering chains, and the three
// test types (equidistribution, collision_free, tuplets). The
// `file:` source (per-family Generator.from_yaml) is intentionally
// NOT supported in this cut; it is used by exactly one sample config
// (example3.yaml) and will land alongside the catalog integration in
// a later phase. The pre-v2 `legacy_file:` source has moved to the
// optional `regpoly-legacy` Python add-on.
//
// Tempering values that are YAML maps with a `random:` key (the
// Python placeholder for "fill in randomly") are silently dropped
// from the Params — create_transformation() will then perform its
// default randomization.
//
// All paths in the YAML are resolved relative to the YAML file's
// parent directory (matches Python's _resolve_path).

#include "combination.h"
#include "generator.h"
#include "params.h"
#include "search_types.h"
#include "seek_search.h"
#include "transformation.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace regpoly_yaml_config {

struct ComponentSpec {
    enum class Source { Inline, Same };
    Source source = Source::Inline;

    // Inline source.
    std::string inline_family;        // e.g. "TGFSRGen"
    Params common_params;             // shared by every entry
    std::vector<Params> per_gen_params;

    // Tempering chain (built directly into Transformation specs).
    struct TemperStep {
        std::string type;             // e.g. "tempMK"
        Params params;                // map-valued keys (random placeholders)
                                      // are dropped before reaching here.
    };
    std::vector<TemperStep> tempering;
};

struct SeekConfig {
    int64_t seed1 = -1;
    int64_t seed2 = 0;
    int Lmax = 32;
    int nbtries = 1;
    std::string output_dir;
    std::vector<ComponentSpec> components;
    std::vector<SeekTestSpec> tests;
};

// Throws std::runtime_error on YAML parse / schema errors with a
// message that includes the filename.
SeekConfig load_seek_config(const std::string& filename);

// Build a fully-initialised Combination from the SeekConfig:
// resolves inline/same into actual Generator objects via
// the existing factories, attaches the tempering chain to each
// component, and calls reset(). Returns the Combination plus the
// owning vectors of generators / transformations (kept alive in case
// callers want to inspect them — Combination::add_gen/add_trans
// already deep-copy into the component's internal pool, so these
// vectors are not strictly required for the search to run).
struct BuiltSearch {
    std::unique_ptr<Combination> combination;
    std::vector<std::vector<std::unique_ptr<Generator>>> per_component_gens;
    std::vector<std::vector<std::unique_ptr<Transformation>>> per_component_trans;
};
BuiltSearch build_search(const SeekConfig& cfg);

}  // namespace regpoly_yaml_config
