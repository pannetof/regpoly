// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

/**
 * @file seek_config.h
 * @brief Focused YAML loader for the regpoly-cli search driver.
 * @ingroup core
 *
 * Phase 4.2: focused YAML loader for the `regpoly-cli search`
 * driver.
 *
 * Mirrors the subset of `regpoly.search.seek.Seek.from_yaml` needed
 * by the C++ CLI: parse a seek-style search config (`search` /
 * `components` / `tests`), build the corresponding `Combination` +
 * tempering chains via the existing factories, and hand off to
 * `run_seek_search()`.
 *
 * The schema supported here is a strict subset of the Python
 * loader's: two component sources (inline `family` / `common` /
 * `generators` and the literal `"same"` sentinel), simple tempering
 * chains, and the three test types (`equidistribution`,
 * `collision_free`, `tuplets`). The `file:` source (per-family
 * `Generator.from_yaml`) is intentionally NOT supported in this
 * cut; it is used by exactly one sample config (`example3.yaml`)
 * and will land alongside the catalog integration in a later phase.
 * The pre-v2 `legacy_file:` source has moved to the optional
 * `regpoly-legacy` Python add-on.
 *
 * Tempering values that are YAML maps with a `random:` key (the
 * Python placeholder for "fill in randomly") are silently dropped
 * from the `Params` — `create_transformation()` will then perform
 * its default randomisation.
 *
 * All paths in the YAML are resolved relative to the YAML file's
 * parent directory (matches Python's `_resolve_path`).
 *
 * There is no Python counterpart for this header directly: Python
 * loads YAML via PyYAML inside `regpoly.search.seek.Seek.from_yaml`
 * rather than calling into this loader. Use the C++ loader from
 * `regpoly-cli search` or from C++-only consumers.
 */

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

namespace regpoly::yaml_config {

using regpoly::core::Combination;
using regpoly::core::Generator;
using regpoly::core::Params;
using regpoly::core::SearchProgress;
using regpoly::core::SeekIterResult;
using regpoly::core::SeekTestSpec;
using regpoly::core::Transformation;

/**
 * @brief Parsed `components:` entry of a seek-search YAML config.
 *
 * Two component sources are supported: `Inline` (the component
 * specifies its own family, common parameters, and per-generator
 * parameter overrides) and `Same` (reuse the previous slot's pool).
 *
 * @ingroup core
 */
struct ComponentSpec {
    /// Source kind for this component slot.
    enum class Source {
        Inline,  ///< Inline component with its own family / params.
        Same,    ///< Re-use the previous slot's generator pool.
    };
    Source source = Source::Inline;     ///< Source kind.

    // Inline source.
    std::string inline_family;          ///< Canonical family name (e.g. `"TGFSRGen"`).
    Params common_params;               ///< Parameters shared by every generator entry.
    std::vector<Params> per_gen_params; ///< Per-generator parameter overrides.

    /**
     * @brief One step in this component's tempering chain.
     *
     * @ingroup core
     */
    struct TemperStep {
        std::string type;             ///< Transformation type name (e.g. `"tempMK"`).
        /**
         * @brief Tempering parameters.
         *
         * Map-valued keys (Python's `random:` placeholder) are
         * dropped before reaching here; the factory then performs
         * its default randomisation.
         */
        Params params;
    };
    std::vector<TemperStep> tempering; ///< Tempering chain for this slot.
};

/**
 * @brief Top-level seek-search YAML configuration.
 *
 * Carries the global search parameters (`seed1` / `seed2` / `Lmax` /
 * `nbtries` / `output_dir`), the per-component specs, and the test
 * stack to run on each combo.
 *
 * @ingroup core
 */
struct SeekConfig {
    int64_t seed1 = -1;                       ///< Primary search seed (`-1` = OS entropy).
    int64_t seed2 = 0;                        ///< Secondary search seed.
    int Lmax = 32;                            ///< Maximum output word width across components.
    int nbtries = 1;                          ///< Number of tries per combo.
    std::string output_dir;                   ///< Output directory (relative to the YAML file).
    std::vector<ComponentSpec> components;    ///< Per-component specs.
    std::vector<SeekTestSpec> tests;          ///< Test stack run per combo.
};

/**
 * @brief Parse a seek-search YAML config from disk.
 *
 * @code{.cpp}
 *   namespace yc = regpoly::yaml_config;
 *   auto cfg   = yc::load_seek_config("shared/yaml/equidist/mt19937.yaml");
 *   auto built = yc::build_search(cfg);
 *   regpoly::core::Combination& comb = *built.combination;
 *   auto combined = regpoly::core::build_combined_from_combination(comb);
 *   // ... drive the search via run_seek_search or analyse `combined` directly.
 * @endcode
 *
 * @param filename  Path to the YAML file.
 * @return          Parsed `SeekConfig`.
 * @throws std::runtime_error  On YAML parse / schema errors. The
 *                             message includes the filename.
 */
SeekConfig load_seek_config(const std::string& filename);

/**
 * @brief Owning result of `build_search`.
 *
 * Carries the constructed `Combination` plus the owning vectors of
 * generators / transformations. `Combination::add_gen` /
 * `add_trans` already deep-copy into each component's internal
 * pool, so these vectors are not strictly required for the search
 * to run — they are kept alive in case callers want to inspect
 * them.
 *
 * @ingroup core
 */
struct BuiltSearch {
    std::unique_ptr<Combination> combination;                                       ///< Constructed iterator.
    std::vector<std::vector<std::unique_ptr<Generator>>> per_component_gens;        ///< Owning generator pools.
    std::vector<std::vector<std::unique_ptr<Transformation>>> per_component_trans;  ///< Owning tempering chains.
};

/**
 * @brief Build a fully-initialised `Combination` from a `SeekConfig`.
 *
 * Resolves `Inline` / `Same` into actual `Generator` objects via
 * the existing factories, attaches the tempering chain to each
 * component, and calls `reset()`.
 *
 * @param cfg  Parsed seek configuration.
 * @return     The constructed `Combination` plus the owning
 *             per-component generator / tempering vectors.
 */
BuiltSearch build_search(const SeekConfig& cfg);

}  // namespace regpoly::yaml_config
