// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "seek_config.h"

#include "equidistribution_method.h"
#include "factory.h"
#include "legacy_reader.h"

#include <yaml-cpp/yaml.h>

#include <cctype>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace regpoly_yaml_config {

namespace {

namespace fs = std::filesystem;

// Try parsing a scalar string as int64 — handles decimal AND 0x… hex
// (yaml-cpp leaves the latter as a string). Mirrors catalog.cpp's
// helper of the same name.
bool try_parse_int(const std::string& s, int64_t& out) {
    if (s.empty()) return false;
    try {
        if (s.size() > 2 && s[0] == '0' && (s[1] == 'x' || s[1] == 'X')) {
            uint64_t u = std::stoull(s.substr(2), nullptr, 16);
            out = static_cast<int64_t>(u);
            return true;
        }
        size_t pos = 0;
        if (s[0] == '-' || std::isdigit(static_cast<unsigned char>(s[0]))) {
            int64_t v = std::stoll(s, &pos, 10);
            if (pos == s.size()) { out = v; return true; }
        }
    } catch (...) {}
    return false;
}

// Apply one YAML-derived value into a Params, mirroring the Python
// loader's lenient int / string / list coercion. Map-valued nodes
// (e.g. `{ random: bitvect, bits: w }` for tempering) are silently
// dropped — create_transformation then falls back to its default
// randomization for that key.
void apply_yaml_to_params(Params& dst,
                          const std::string& key,
                          const YAML::Node& node) {
    if (node.IsMap()) {
        // Random-value placeholder; drop and let the factory randomize.
        return;
    }
    if (node.IsSequence()) {
        // Decide between int_vec and uint_vec by overflow detection.
        bool all_int = true;
        std::vector<int64_t> raw;
        for (const auto& el : node) {
            if (!el.IsScalar()) { all_int = false; break; }
            int64_t v = 0;
            std::string s = el.as<std::string>();
            if (!try_parse_int(s, v)) {
                try { v = el.as<int64_t>(); }
                catch (...) { all_int = false; break; }
            }
            raw.push_back(v);
        }
        if (!all_int) {
            // Stringify defensively — none of the supported tempering /
            // generator params take heterogeneous lists today.
            YAML::Emitter em; em << node;
            dst.set_string(key, em.c_str());
            return;
        }
        bool needs_uint = false;
        for (auto x : raw) {
            if (x > static_cast<int64_t>(std::numeric_limits<int>::max())
             || x < static_cast<int64_t>(std::numeric_limits<int>::min())) {
                needs_uint = true; break;
            }
        }
        if (needs_uint) {
            std::vector<uint64_t> vu;
            vu.reserve(raw.size());
            for (auto x : raw) vu.push_back(static_cast<uint64_t>(x));
            dst.set_uint_vec(key, vu);
        } else {
            std::vector<int> v32;
            v32.reserve(raw.size());
            for (auto x : raw) v32.push_back(static_cast<int>(x));
            dst.set_int_vec(key, v32);
        }
        return;
    }
    if (!node.IsScalar()) return;

    const auto& tag = node.Tag();
    std::string s = node.as<std::string>();

    // Bool detection.
    if (tag == "tag:yaml.org,2002:bool"
        || s == "true" || s == "false"
        || s == "True" || s == "False") {
        try { dst.set_bool(key, node.as<bool>()); return; }
        catch (...) {}
    }

    // Int detection (covers both decimal and 0x… hex).
    int64_t iv = 0;
    if (try_parse_int(s, iv)) {
        dst.set_int(key, iv);
        return;
    }

    // Fall back to string — preserves things like family names.
    dst.set_string(key, s);
}

std::string resolve_path(const std::string& path,
                         const std::string& base_dir) {
    if (path.empty()) return path;
    fs::path p(path);
    if (p.is_absolute()) return path;
    fs::path base(base_dir);
    return (base / p).lexically_normal().string();
}

void parse_search_block(const YAML::Node& search,
                        SeekConfig& out,
                        const std::string& filename) {
    if (search && !search.IsMap()) {
        throw std::runtime_error(filename + ": 'search' must be a mapping");
    }
    if (!search) return;

    if (search["seed"]) {
        const auto& s = search["seed"];
        if (!s.IsSequence() || s.size() < 1 || s.size() > 2) {
            throw std::runtime_error(filename
                + ": search.seed must be a list of two ints");
        }
        out.seed1 = s[0].as<int64_t>(-1);
        out.seed2 = s.size() >= 2 ? s[1].as<int64_t>(0) : 0;
    }
    if (search["Lmax"]) {
        out.Lmax = search["Lmax"].as<int>(32);
        if (out.Lmax <= 0) {
            throw std::runtime_error(filename
                + ": search.Lmax must be a positive int");
        }
    }
    if (search["nbtries"]) {
        out.nbtries = search["nbtries"].as<int>(1);
        if (out.nbtries < 1) out.nbtries = 1;
    }
    if (search["output_dir"]) {
        out.output_dir = search["output_dir"].as<std::string>("");
    }
}

ComponentSpec parse_component(const YAML::Node& comp_node,
                              const std::string& base_dir,
                              const std::string& filename,
                              int comp_index) {
    ComponentSpec out;

    if (!comp_node.IsMap()) {
        throw std::runtime_error(filename + ": components["
            + std::to_string(comp_index) + "] must be a mapping");
    }

    auto gens_node = comp_node["generators"];
    if (!gens_node) {
        throw std::runtime_error(filename + ": components["
            + std::to_string(comp_index) + "] missing 'generators'");
    }

    if (gens_node.IsScalar()) {
        std::string s = gens_node.as<std::string>("");
        if (s == "same") {
            out.source = ComponentSpec::Source::Same;
        } else {
            throw std::runtime_error(filename + ": components["
                + std::to_string(comp_index)
                + "].generators is a scalar but not 'same' (got '"
                + s + "')");
        }
    } else if (gens_node.IsMap()) {
        if (gens_node["legacy_file"]) {
            out.source = ComponentSpec::Source::LegacyFile;
            std::string raw = gens_node["legacy_file"].as<std::string>("");
            out.legacy_file_path = resolve_path(raw, base_dir);
        } else if (gens_node["family"]) {
            out.source = ComponentSpec::Source::Inline;
            out.inline_family = gens_node["family"].as<std::string>("");

            // family-level params (everything except family, common,
            // generators) feed into common alongside the explicit
            // common: block. Mirror Python's family_params / common
            // merge.
            for (auto kv : gens_node) {
                std::string k = kv.first.as<std::string>();
                if (k == "family" || k == "common" || k == "generators") {
                    continue;
                }
                apply_yaml_to_params(out.common_params, k, kv.second);
            }
            if (gens_node["common"] && gens_node["common"].IsMap()) {
                for (auto kv : gens_node["common"]) {
                    std::string k = kv.first.as<std::string>();
                    apply_yaml_to_params(out.common_params, k, kv.second);
                }
            }
            auto entries = gens_node["generators"];
            if (!entries || !entries.IsSequence() || entries.size() == 0) {
                throw std::runtime_error(filename + ": components["
                    + std::to_string(comp_index)
                    + "].generators.generators must be a non-empty list");
            }
            for (auto entry : entries) {
                if (!entry.IsMap()) {
                    throw std::runtime_error(filename + ": components["
                        + std::to_string(comp_index)
                        + "].generators.generators[i] must be a mapping");
                }
                Params merged;
                // Copy common into merged first, then overlay entry.
                for (const auto& kv : out.common_params.ints())
                    merged.set_int(kv.first, kv.second);
                for (const auto& kv : out.common_params.bools())
                    merged.set_bool(kv.first, kv.second);
                for (const auto& kv : out.common_params.strings())
                    merged.set_string(kv.first, kv.second);
                for (const auto& kv : out.common_params.int_vecs())
                    merged.set_int_vec(kv.first, kv.second);
                for (const auto& kv : out.common_params.uint_vecs())
                    merged.set_uint_vec(kv.first, kv.second);
                for (auto kv : entry) {
                    std::string k = kv.first.as<std::string>();
                    apply_yaml_to_params(merged, k, kv.second);
                }
                out.per_gen_params.push_back(std::move(merged));
            }
        } else if (gens_node["file"]) {
            throw std::runtime_error(filename + ": components["
                + std::to_string(comp_index)
                + "].generators.file is not yet supported by the C++ CLI "
                  "(use the Python `regpoly` CLI for per-family YAML loading)");
        } else {
            throw std::runtime_error(filename + ": components["
                + std::to_string(comp_index)
                + "].generators must contain one of: "
                  "legacy_file, family, file (or be the string 'same')");
        }
    } else {
        throw std::runtime_error(filename + ": components["
            + std::to_string(comp_index)
            + "].generators must be a mapping or the string 'same'");
    }

    // Tempering chain — optional list.
    if (comp_node["tempering"]) {
        const auto& tnode = comp_node["tempering"];
        if (!tnode.IsSequence()) {
            throw std::runtime_error(filename + ": components["
                + std::to_string(comp_index)
                + "].tempering must be a list");
        }
        int ti = -1;
        for (auto step_node : tnode) {
            ++ti;
            if (!step_node.IsMap() || !step_node["type"]) {
                throw std::runtime_error(filename + ": components["
                    + std::to_string(comp_index) + "].tempering["
                    + std::to_string(ti) + "] needs a 'type' key");
            }
            ComponentSpec::TemperStep step;
            step.type = step_node["type"].as<std::string>("");
            for (auto kv : step_node) {
                std::string k = kv.first.as<std::string>();
                if (k == "type") continue;
                apply_yaml_to_params(step.params, k, kv.second);
            }
            out.tempering.push_back(std::move(step));
        }
    }
    return out;
}

SeekTestSpec parse_test(const YAML::Node& tn,
                        int Lmax,
                        const std::string& filename,
                        int test_index) {
    if (!tn.IsMap() || !tn["type"]) {
        throw std::runtime_error(filename + ": tests["
            + std::to_string(test_index) + "] needs a 'type' key");
    }
    std::string type = tn["type"].as<std::string>("");
    SeekTestSpec spec;

    if (type == "equidistribution") {
        std::string method = tn["method"]
            ? tn["method"].as<std::string>("matricial")
            : "matricial";
        // Method-name validation lives in MethodRegistry; the YAML
        // parser no longer owns the list.
        if (!MethodRegistry::has(method)) {
            throw std::runtime_error(filename + ": tests["
                + std::to_string(test_index)
                + "].method '" + method + "' is unknown");
        }
        spec.kind = SeekTestKind::Equidistribution;
        spec.method_name = method;
        spec.eq_L_max_test = Lmax;
        int64_t mse = INT_MAX;
        if (tn["max_gap_sum"]) {
            int64_t v = tn["max_gap_sum"].as<int64_t>(INT_MAX);
            if (v > INT_MAX) v = INT_MAX;
            mse = v;
        }
        spec.eq_mse = static_cast<int>(mse);
        spec.eq_delta.assign(Lmax + 1, INT_MAX);
        if (tn["delta"]) {
            const auto& dn = tn["delta"];
            if (!dn.IsSequence()) {
                throw std::runtime_error(filename + ": tests["
                    + std::to_string(test_index)
                    + "].delta must be a list of {from,to,max}");
            }
            for (auto rule : dn) {
                if (!rule.IsMap()
                    || !rule["from"] || !rule["to"] || !rule["max"]) {
                    throw std::runtime_error(filename + ": tests["
                        + std::to_string(test_index)
                        + "].delta entries need from/to/max");
                }
                int lo = rule["from"].as<int>(1);
                int hi = rule["to"].as<int>(Lmax);
                int64_t mx = rule["max"].as<int64_t>(INT_MAX);
                if (mx > INT_MAX) mx = INT_MAX;
                int hi_clamped = std::min(hi, Lmax);
                for (int l = lo; l <= hi_clamped; ++l) {
                    if (l < 0 || l > Lmax) continue;
                    spec.eq_delta[l] = static_cast<int>(mx);
                }
            }
        }
        return spec;
    }

    if (type == "collision_free") {
        spec.kind = SeekTestKind::CollisionFree;
        return spec;
    }

    if (type == "tuplets") {
        spec.kind = SeekTestKind::Tuplets;
        std::vector<int> dims;
        if (tn["dimensions"] && tn["dimensions"].IsSequence()) {
            for (auto x : tn["dimensions"]) dims.push_back(x.as<int>(0));
        }
        spec.tup_d = static_cast<int>(dims.size());
        spec.tup_h.assign(1, 0);
        for (int d : dims) spec.tup_h.push_back(d);
        spec.tup_threshold = tn["threshold"]
            ? tn["threshold"].as<double>(0.0) : 0.0;
        std::string tt = tn["test_type"]
            ? tn["test_type"].as<std::string>("max") : "max";
        spec.tup_testtype = (tt == "max") ? 0 : 1;
        return spec;
    }

    throw std::runtime_error(filename + ": tests["
        + std::to_string(test_index)
        + "] unknown type '" + type + "'");
}

// Build a Generator pool for a single inline component spec — calls
// create_generator(family, merged_params, Lmax) for every entry.
std::vector<std::unique_ptr<Generator>>
build_inline_pool(const ComponentSpec& spec, int Lmax,
                  const std::string& context) {
    std::vector<std::unique_ptr<Generator>> out;
    for (size_t i = 0; i < spec.per_gen_params.size(); ++i) {
        try {
            out.push_back(create_generator(spec.inline_family,
                                           spec.per_gen_params[i], Lmax));
        } catch (const std::exception& exc) {
            throw std::runtime_error(context + " generator[" + std::to_string(i)
                + "] (" + spec.inline_family + "): " + exc.what());
        }
    }
    return out;
}

// Build a Generator pool from a legacy_file source.
std::vector<std::unique_ptr<Generator>>
build_legacy_pool(const ComponentSpec& spec, int Lmax,
                  const std::string& context) {
    std::vector<std::unique_ptr<Generator>> out;
    std::vector<regpoly_legacy::LegacyGeneratorSpec> raw;
    try {
        raw = regpoly_legacy::read_generator_specs(spec.legacy_file_path, Lmax);
    } catch (const std::exception& exc) {
        throw std::runtime_error(context + " legacy_file '"
            + spec.legacy_file_path + "': " + exc.what());
    }
    for (size_t i = 0; i < raw.size(); ++i) {
        try {
            out.push_back(create_generator(raw[i].family,
                                           raw[i].params, Lmax));
        } catch (const std::exception& exc) {
            throw std::runtime_error(context + " legacy_file generator["
                + std::to_string(i) + "] (" + raw[i].family
                + "): " + exc.what());
        }
    }
    return out;
}

// Apply default w from the inline component's common_params (if any)
// to a tempering Params block when 'w' is missing. Mirrors catalog
// .cpp's dry_build_generator.
void apply_w_default(Params& tp, const ComponentSpec& comp) {
    if (tp.has("w")) return;
    if (comp.common_params.has("w")) {
        tp.set_int("w", comp.common_params.get_int("w"));
    }
}

}  // namespace

SeekConfig load_seek_config(const std::string& filename) {
    YAML::Node root;
    try {
        root = YAML::LoadFile(filename);
    } catch (const std::exception& exc) {
        throw std::runtime_error(filename + ": YAML load failed: "
            + exc.what());
    }
    if (!root || !root.IsMap()) {
        throw std::runtime_error(filename
            + ": top-level YAML must be a mapping");
    }

    SeekConfig out;
    fs::path fp(filename);
    fs::path abs = fp.is_absolute() ? fp : fs::absolute(fp);
    std::string base_dir = abs.parent_path().string();

    parse_search_block(root["search"], out, filename);

    auto comps = root["components"];
    if (!comps || !comps.IsSequence() || comps.size() == 0) {
        throw std::runtime_error(filename
            + ": 'components' must be a non-empty list");
    }
    int idx = -1;
    for (auto cn : comps) {
        ++idx;
        out.components.push_back(parse_component(cn, base_dir, filename, idx));
    }

    auto tests = root["tests"];
    if (tests) {
        if (tests.IsMap() && tests["file"]) {
            throw std::runtime_error(filename
                + ": tests.file is not supported by the C++ CLI; inline the "
                  "tests in the search config or use the Python `regpoly` CLI");
        }
        if (!tests.IsSequence()) {
            throw std::runtime_error(filename
                + ": 'tests' must be a list");
        }
        int ti = -1;
        for (auto tn : tests) {
            ++ti;
            out.tests.push_back(parse_test(tn, out.Lmax, filename, ti));
        }
    }

    return out;
}

BuiltSearch build_search(const SeekConfig& cfg) {
    BuiltSearch built;
    const int J = static_cast<int>(cfg.components.size());
    if (J <= 0) {
        throw std::runtime_error("build_search: components list is empty");
    }
    built.combination = std::make_unique<Combination>(J, cfg.Lmax);
    built.per_component_gens.resize(J);
    built.per_component_trans.resize(J);

    for (int j = 0; j < J; ++j) {
        const auto& spec = cfg.components[j];
        std::string ctx = "components[" + std::to_string(j) + "]";

        if (spec.source == ComponentSpec::Source::Same) {
            if (j == 0) {
                throw std::runtime_error(ctx
                    + ": 'same' is invalid on the first component");
            }
            built.combination->component(j).share_pool_with(
                built.combination->component(j - 1));
        } else {
            std::vector<std::unique_ptr<Generator>> pool;
            if (spec.source == ComponentSpec::Source::LegacyFile) {
                pool = build_legacy_pool(spec, cfg.Lmax, ctx);
            } else {
                pool = build_inline_pool(spec, cfg.Lmax, ctx);
            }
            if (pool.empty()) {
                throw std::runtime_error(ctx
                    + ": resolved generator pool is empty");
            }
            for (const auto& g : pool) {
                built.combination->component(j).add_gen(*g);
            }
            built.per_component_gens[j] = std::move(pool);
        }

        // Tempering chain (always processed — even for 'same' the per-
        // component tempering applies on top of the shared pool).
        for (size_t ti = 0; ti < spec.tempering.size(); ++ti) {
            const auto& step = spec.tempering[ti];
            Params tp;
            // Copy step.params into tp.
            for (const auto& kv : step.params.ints())
                tp.set_int(kv.first, kv.second);
            for (const auto& kv : step.params.bools())
                tp.set_bool(kv.first, kv.second);
            for (const auto& kv : step.params.strings())
                tp.set_string(kv.first, kv.second);
            for (const auto& kv : step.params.int_vecs())
                tp.set_int_vec(kv.first, kv.second);
            for (const auto& kv : step.params.uint_vecs())
                tp.set_uint_vec(kv.first, kv.second);
            apply_w_default(tp, spec);
            std::unique_ptr<Transformation> t;
            try {
                t = create_transformation(step.type, tp);
            } catch (const std::exception& exc) {
                throw std::runtime_error(ctx + ".tempering["
                    + std::to_string(ti) + "] (" + step.type + "): "
                    + exc.what());
            }
            built.combination->component(j).add_trans(*t);
            built.per_component_trans[j].push_back(std::move(t));
        }
    }

    if (!built.combination->reset()) {
        throw std::runtime_error(
            "build_search: Combination::reset() returned false — "
            "no valid combo (empty pool, or shared-pool C(n,k) "
            "exhausted before a single placement?)");
    }
    return built;
}

}  // namespace regpoly_yaml_config
