// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "params.h"   // StructMap typedef

// Phase 3.2: paper-centric published-generators catalog in C++.
//
// Mirrors the data model of regpoly/library/__init__.py:
//
//   Author  → bibliographic author
//   CatalogGenerator → one generator entry inside a paper
//                      (renamed from "Generator" to avoid colliding
//                      with the runtime Generator class)
//   Paper   → one bibliographic entry holding a list of generators
//   Catalog → a directory of *.yaml files indexed by paper id
//
// The catalog reader/writer is implemented in C++ and exposed to
// Python through pybind11; the regpoly Python package keeps a thin
// shim at regpoly/library/__init__.py that just re-exports the
// pybind11 types.

namespace regpoly_catalog {

// ── Variant value type for component params + tempering steps ─────────
//
// YAML parameter values come in three shapes the catalog needs to
// preserve faithfully: ints (incl. bitmask hex literals stored as
// negative-when-MSB or wrapped to int64), strings (e.g. "0x9908b0df"
// — yaml-cpp parses hex as a string, the runtime factory converts),
// bools, and homogeneous int lists (Tausworthe `poly`, MELG `nb_lag`).
// We store them in a tagged union since pybind11 transparently
// converts std::variant<int64_t,…> to native Python.

enum class ParamKind : uint8_t {
    Int = 0,
    String = 1,
    Bool = 2,
    IntList = 3,
    StructMap = 4,
};

struct ParamValue {
    ParamKind kind = ParamKind::String;
    int64_t int_val = 0;
    std::string string_val;
    bool bool_val = false;
    std::vector<int64_t> int_list_val;
    StructMap struct_map_val;

    static ParamValue make_int(int64_t v) {
        ParamValue p;
        p.kind = ParamKind::Int;
        p.int_val = v;
        return p;
    }
    static ParamValue make_string(std::string v) {
        ParamValue p;
        p.kind = ParamKind::String;
        p.string_val = std::move(v);
        return p;
    }
    static ParamValue make_bool(bool v) {
        ParamValue p;
        p.kind = ParamKind::Bool;
        p.bool_val = v;
        return p;
    }
    static ParamValue make_int_list(std::vector<int64_t> v) {
        ParamValue p;
        p.kind = ParamKind::IntList;
        p.int_list_val = std::move(v);
        return p;
    }
    static ParamValue make_struct_map(StructMap v) {
        ParamValue p;
        p.kind = ParamKind::StructMap;
        p.struct_map_val = std::move(v);
        return p;
    }
};

using ParamMap = std::map<std::string, ParamValue>;  // ordered for determinism

struct TemperingStep {
    std::string type;            // e.g. "tempMK"
    ParamMap params;             // remaining keys (b, c, eta, mu, w, …)
};

struct Component {
    std::string family;
    int L = 0;
    ParamMap params;
    std::vector<TemperingStep> tempering;
};

struct Author {
    std::string family;
    std::string given;

    std::string display() const {
        return given.empty() ? family : (family + ", " + given);
    }
    std::string short_name() const { return family; }
};

struct CatalogGenerator {
    std::string id;              // slug, e.g. "lfsr113"
    std::string display;
    std::string family;          // C++ generator family name
    std::string target;          // "tested_generator" or "primitive_generator"
    bool combined = false;
    int Lmax = 0;
    std::vector<Component> components;
    std::string notes_md;
    bool starred = false;
    std::vector<std::string> errors;

    bool valid() const { return errors.empty(); }
};

struct Paper {
    std::string id;
    std::vector<Author> authors;
    int year = 0;
    std::string title;
    std::string venue;
    std::string volume;
    std::string issue;
    std::string pages;
    std::string doi;
    std::string pdf;
    std::string bibkey;
    std::string abstract_md;
    std::string notes_md;
    std::vector<std::string> tags;
    bool starred = false;
    bool deferred = false;
    std::vector<CatalogGenerator> generators;
    std::string source_path;     // absolute path
    double source_mtime = 0.0;   // POSIX mtime
    std::vector<std::string> errors;

    bool valid() const {
        if (!errors.empty()) return false;
        for (const auto& g : generators) {
            if (!g.valid()) return false;
        }
        return true;
    }

    std::string author_list_short() const;
    std::string display() const;
    std::string acmtrans_citation() const;
};

class Catalog {
public:
    explicit Catalog(std::string library_dir);

    // Load every *.yaml under library_dir (skipping *_params.yaml
    // cross-check fixtures). Errors at parse time are logged into
    // each Paper's errors[] field rather than aborting the load.
    void load();

    // Re-stat every *.yaml; reload only those whose mtime advanced
    // (and drop papers whose files vanished).
    void reload_if_stale();

    // Filter views.
    struct PapersFilter {
        std::optional<bool> starred;
        std::optional<std::string> tag;
        bool include_invalid = false;

        PapersFilter() = default;
    };
    std::vector<Paper> papers() const { return papers(PapersFilter()); }
    std::vector<Paper> papers(const PapersFilter& filter) const;

    // Single-paper / single-generator lookup.
    std::optional<Paper> paper(const std::string& paper_id) const;
    // Returns (paper, generator); the generator is a copy of paper.generators[idx].
    std::optional<std::pair<Paper, CatalogGenerator>>
    generator(const std::string& gen_id) const;

    std::vector<std::pair<Paper, CatalogGenerator>>
    all_generators(const std::optional<std::string>& family = {}) const;

    const std::string& library_dir() const { return library_dir_; }

private:
    std::string library_dir_;
    std::map<std::string, Paper> papers_;          // ordered by id
    // generator id -> (paper id, index within paper.generators)
    std::unordered_map<std::string, std::pair<std::string, int>> generator_index_;

    void insert(Paper paper);
    void drop(const std::string& paper_id);
};

// Stable short hash (16 hex chars) of one component config. Used by
// the web UI to group equivalent configurations. Deterministic across
// runs but not bit-equivalent to the prior Python json.dumps + sha256
// implementation; no caller currently persists the value, so the
// change is safe.
std::string config_hash(
    const std::string& family,
    const ParamMap& params,
    const std::vector<TemperingStep>& tempering);

// Parse a single paper YAML file (used internally by Catalog::load
// and exposed for unit tests).
Paper parse_paper(const std::string& path);

// Phase 4.3: textual append of one tested-generator entry into an
// existing paper's `generators:` list.
//
// Read the source file (must match the schema documented in
// regpoly/io/tested_generator.py: top-level `generator` + `tempering`
// for single-component, or `components: [{generator, tempering}, ...]`
// for multi). Produces a YAML mapping appended to the paper file's
// `generators:` list — the new entry has the keys {id, display,
// family, target, combined, Lmax, components} that the Catalog loader
// expects.
//
// Constraints (deliberate to keep the writer safe):
//   * The paper file must exist under library_dir.
//   * `generators:` must be the LAST top-level key in the paper file
//     (so a textual append at end-of-file lands inside the list).
//     If not, throws with a clear remediation message.
//   * gen_id must be unique within the paper (and the catalog).
//
// Returns the path of the modified paper file. Throws
// std::runtime_error on any precondition failure or I/O error; the
// paper file is left untouched in that case.
std::string publish_tested_generator(
    const std::string& library_dir,
    const std::string& paper_id,
    const std::string& gen_id,
    const std::string& display,
    const std::string& tested_generator_file,
    const std::string& target = "tested_generator",
    bool starred = false);

// Skip rule for the *_params.yaml MTToolBox cross-check fixtures
// living alongside paper-organised entries in docs/library/.
bool is_cross_check_yaml(const std::string& filename);

}  // namespace regpoly_catalog
