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

/**
 * @file catalog.h
 * @brief Paper-centric published-generators catalog.
 *
 * Defines the data model and on-disk reader for the catalog of
 * published parameter sets. Mirrors the Python data model surfaced
 * in `regpoly.library.__init__`:
 *
 * | C++ symbol          | Python symbol                  |
 * |---------------------|--------------------------------|
 * | `Author`            | `regpoly.library.Author`       |
 * | `CatalogGenerator`  | `regpoly.library.Generator`    |
 * | `Paper`             | `regpoly.library.Paper`        |
 * | `Catalog`           | `regpoly.library.Catalog`      |
 *
 * The reader/writer is implemented in C++ and exposed to Python
 * through pybind11; the `regpoly` Python package keeps only a thin
 * shim that re-exports the pybind11 types.
 *
 * @see [`regpoly.library.Catalog`](../../python/regpoly.library.md)
 */

namespace regpoly::library {

using regpoly::core::StructMap;

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

/**
 * @brief Registry of papers and the generators they contain.
 *
 * Walks `library_dir/*.yaml`, parses each as a `Paper`, validates
 * every generator entry against the runtime factory, and indexes the
 * result by paper id and generator id. The lookup methods return
 * value copies — papers and generators are small enough that callers
 * don't need to worry about lifetimes.
 *
 * @see [`regpoly.library.Catalog`](../../python/regpoly.library.md) — the Python wrapper.
 */
class Catalog {
public:
    /**
     * @brief Construct a catalog rooted at `library_dir`.
     *
     * The constructor does not touch the filesystem. Call `load()`
     * (or `reload_if_stale()`) explicitly.
     *
     * @param library_dir  Directory containing paper YAML files.
     */
    explicit Catalog(std::string library_dir);

    /**
     * @brief Read every `*.yaml` under `library_dir`.
     *
     * Skips files matching the `*_params.yaml` cross-check fixture
     * pattern. Errors at parse time are logged into each `Paper`'s
     * `errors[]` field rather than aborting the load — call
     * `papers(filter)` with `include_invalid = false` to filter
     * the failures back out.
     *
     * @throws std::runtime_error  If `library_dir` does not exist.
     */
    void load();

    /**
     * @brief Re-stat every YAML and reload only the ones that changed.
     *
     * Cheap to call from request handlers. Papers whose source file
     * disappeared between calls are dropped.
     */
    void reload_if_stale();

    /// Filter spec for the `papers(filter)` query.
    struct PapersFilter {
        std::optional<bool> starred;     ///< Filter by `starred:` flag.
        std::optional<std::string> tag;  ///< Filter by `tags:` membership.
        bool include_invalid = false;    ///< Include papers with `errors[]`.

        PapersFilter() = default;
    };

    /** @brief Return every (valid) paper in catalog order. */
    std::vector<Paper> papers() const { return papers(PapersFilter()); }

    /**
     * @brief Return the subset of papers matching `filter`.
     * @param filter  Filter spec (starred, tag, include_invalid).
     * @returns       Papers in catalog (load) order.
     */
    std::vector<Paper> papers(const PapersFilter& filter) const;

    /**
     * @brief Look up one paper by id.
     * @param paper_id  Paper id (YAML file stem).
     * @returns The matching `Paper`, or `std::nullopt` if unknown.
     */
    std::optional<Paper> paper(const std::string& paper_id) const;

    /**
     * @brief Look up one generator by id, plus its containing paper.
     * @param gen_id  Generator id (unique across the catalog).
     * @returns  `(paper, generator)`, or `std::nullopt` if unknown.
     */
    std::optional<std::pair<Paper, CatalogGenerator>>
    generator(const std::string& gen_id) const;

    /**
     * @brief Enumerate every generator, optionally filtered by family.
     * @param family  Canonical `-Gen` class name. `std::nullopt` →
     *                no filter.
     * @returns Catalog-order list of `(paper, generator)` pairs.
     */
    std::vector<std::pair<Paper, CatalogGenerator>>
    all_generators(const std::optional<std::string>& family = {}) const;

    /** @brief Directory the catalog was constructed with. */
    const std::string& library_dir() const { return library_dir_; }

private:
    std::string library_dir_;
    std::map<std::string, Paper> papers_;          // ordered by id
    // generator id -> (paper id, index within paper.generators)
    std::unordered_map<std::string, std::pair<std::string, int>> generator_index_;

    void insert(Paper paper);
    void drop(const std::string& paper_id);
};

/**
 * @brief Stable short hash (16 hex chars) of one component config.
 *
 * Used by the web UI to group equivalent configurations. Deterministic
 * across runs but **not** bit-equivalent to the pre-v2 Python
 * `json.dumps + sha256` output; no caller currently persists the
 * value, so the change is safe.
 *
 * @param family     Canonical generator family name.
 * @param params     Component parameter map.
 * @param tempering  Tempering step list.
 * @returns          16-character lowercase hex string.
 */
std::string config_hash(
    const std::string& family,
    const ParamMap& params,
    const std::vector<TemperingStep>& tempering);

/**
 * @brief Parse a single paper YAML file.
 *
 * Used internally by `Catalog::load()`; exposed in the public header
 * for the catalog unit tests.
 *
 * @param path  Absolute path to the YAML file.
 * @returns     A populated `Paper`. Parse errors land in `paper.errors`
 *              rather than throwing.
 */
Paper parse_paper(const std::string& path);

/**
 * @brief Append one tested-generator entry to an existing paper YAML.
 *
 * Reads the tested-generator source file (schema documented in
 * `regpoly/io/tested_generator.py`: top-level `generator` + `tempering`
 * for single-component, or `components:` list for multi), and appends
 * a YAML mapping to the paper file's `generators:` list. The new
 * entry has the keys `{id, display, family, target, combined, Lmax,
 * components}` that the catalog loader expects.
 *
 * @pre  The paper file exists under `library_dir`.
 * @pre  `generators:` is the LAST top-level key in the paper file
 *       (so a textual append at end-of-file lands inside the list).
 * @pre  `gen_id` is unique within the paper and across the catalog.
 *
 * @param library_dir              Catalog root.
 * @param paper_id                 Paper to append to.
 * @param gen_id                   Slug for the new generator entry.
 * @param display                  Human-readable display name.
 * @param tested_generator_file    Source file path.
 * @param target                   `"tested_generator"` (default) or
 *                                 `"primitive_generator"`.
 * @param starred                  Initial value of the `starred:` flag.
 *
 * @returns  Path of the modified paper file.
 * @throws std::runtime_error  On any precondition failure or I/O
 *                             error. The paper file is left untouched
 *                             in that case.
 */
std::string publish_tested_generator(
    const std::string& library_dir,
    const std::string& paper_id,
    const std::string& gen_id,
    const std::string& display,
    const std::string& tested_generator_file,
    const std::string& target = "tested_generator",
    bool starred = false);

/**
 * @brief Skip rule for `*_params.yaml` MTToolBox cross-check fixtures.
 *
 * Catalog files live next to fixture files under `docs/library/`; the
 * `_params.yaml` suffix marks fixtures. This helper centralises the
 * skip rule so the loader and the unit tests agree.
 *
 * @param filename  File name (not full path).
 * @returns         True iff this file is a cross-check fixture and
 *                  should be skipped during catalog load.
 */
bool is_cross_check_yaml(const std::string& filename);

}  // namespace regpoly::library
