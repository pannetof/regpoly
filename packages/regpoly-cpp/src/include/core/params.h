// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include <string>
#include <map>
#include <unordered_map>
#include <variant>
#include <vector>
#include <cstdint>

/**
 * @file params.h
 * @brief Heterogeneous typed parameter bag passed to factories and registries.
 *
 * Defines `Params` (the runtime parameter bag) plus three helper
 * aliases — `ParamScalar`, `StructEntry`, `StructMap` — used to
 * represent the structured "list of named-arg dicts" parameters
 * such as WELL's `matrices` field. Together these are the lowest
 * common denominator the YAML loader, registries, factories and
 * pybind11 bindings agree on; everything else flows through
 * `ParamSpec` (see param_spec.h) for the declarative half.
 *
 * @ingroup core
 */

// Structured parameter values.
//
// Used for "list of named-arg dicts" parameters such as WELL's `matrices`
// (a map keyed by slot label T0..T7, each value a map of M-class arg names
// like {M, t} or {M, q, t, s, a}). The variant covers the scalar leaves
// PyYAML / pybind11 produce; deeper nesting beyond the two levels here is
// not supported on purpose — keep the shape boring.
namespace regpoly::core {

/**
 * @brief Scalar leaf for the two-level `StructMap` type.
 *
 * Covers the four scalar types PyYAML / pybind11 surface as values
 * inside structured parameters.
 */
using ParamScalar = std::variant<int64_t, uint64_t, std::string, bool>;

/**
 * @brief Inner level of a `StructMap`: one named-arg dict.
 *
 * Maps an argument name (e.g. `"M"`, `"q"`, `"t"`) to its scalar value.
 */
using StructEntry = std::map<std::string, ParamScalar>;

/**
 * @brief Outer level of a structured parameter: map of slot label → `StructEntry`.
 *
 * E.g. WELL's `matrices` parameter is keyed by slot label `"T0"..."T7"`,
 * each value a `StructEntry` describing one M-class argument set.
 * Deeper nesting beyond these two levels is not supported on purpose —
 * keep the shape boring.
 */
using StructMap   = std::map<std::string, StructEntry>;

/**
 * @brief Heterogeneous typed parameter bag passed to factories and registries.
 *
 * `Params` is a runtime container that the YAML loader, factories,
 * registries and pybind11 bindings all agree on. It stores values
 * by string key across six type-segregated maps (ints, bools,
 * strings, int_vecs, uint_vecs, struct_maps); typed setters /
 * getters expose them, and the read-only `*_vecs()` / `ints()` /
 * `bools()` / etc. views let serialisers walk every entry.
 *
 * The structural-vs-non-structural distinction (see ParamSpec) is
 * a property of the *spec* declared by each generator subclass, not
 * of the value stored here. `Params` itself does not track that
 * distinction; it is purely a typed key/value bag.
 *
 * @see regpoly::core::ParamSpec
 *
 * @ingroup core
 */
class Params {
public:
    /**
     * @brief Store an int64 value under `key` (overwrites any prior int entry).
     * @param key  Parameter name.
     * @param val  Value to store.
     */
    void set_int(const std::string& key, int64_t val);

    /**
     * @brief Store a bool value under `key` (overwrites any prior bool entry).
     * @param key  Parameter name.
     * @param val  Value to store.
     */
    void set_bool(const std::string& key, bool val);

    /**
     * @brief Store a string value under `key` (overwrites any prior string entry).
     * @param key  Parameter name.
     * @param val  Value to store.
     */
    void set_string(const std::string& key, const std::string& val);

    /**
     * @brief Store an int-vector value under `key` (overwrites any prior int-vec entry).
     * @param key  Parameter name.
     * @param val  Vector to store (copied).
     */
    void set_int_vec(const std::string& key, const std::vector<int>& val);

    /**
     * @brief Store a uint64-vector value under `key` (overwrites any prior uint-vec entry).
     * @param key  Parameter name.
     * @param val  Vector to store (copied).
     */
    void set_uint_vec(const std::string& key, const std::vector<uint64_t>& val);

    /**
     * @brief Store a structured map under `key` (overwrites any prior struct entry).
     * @param key  Parameter name.
     * @param val  StructMap to store (moved).
     */
    void set_struct_map(const std::string& key, StructMap val);

    /**
     * @brief Retrieve an int64 value or `def` if absent.
     * @param key  Parameter name.
     * @param def  Fallback returned when `key` is not in the int table.
     * @return     The stored value, or `def`.
     */
    int64_t get_int(const std::string& key, int64_t def = 0) const;

    /**
     * @brief Retrieve a bool value or `def` if absent.
     * @param key  Parameter name.
     * @param def  Fallback returned when `key` is not in the bool table.
     * @return     The stored value, or `def`.
     */
    bool get_bool(const std::string& key, bool def = false) const;

    /**
     * @brief Retrieve a string value or `def` if absent.
     * @param key  Parameter name.
     * @param def  Fallback returned when `key` is not in the string table.
     * @return     The stored value, or `def`.
     */
    std::string get_string(const std::string& key, const std::string& def = "") const;

    /**
     * @brief Retrieve an int-vector value, falling back to the uint-vec table if needed.
     *
     * The Python binding stores list parameters as int_vec or
     * uint_vec based on a per-value heuristic that doesn't know
     * the target type; this accessor lets callers read either
     * side regardless. Returns an empty vector if `key` is absent
     * from both tables.
     *
     * @param key  Parameter name.
     * @return     Stored vector, or an empty vector.
     */
    std::vector<int> get_int_vec(const std::string& key) const;

    /**
     * @brief Retrieve a uint-vector value, falling back to the int-vec table if needed.
     *
     * Mirrored fallback: see `get_int_vec`. Returns an empty vector
     * if `key` is absent from both tables.
     *
     * @param key  Parameter name.
     * @return     Stored vector, or an empty vector.
     */
    std::vector<uint64_t> get_uint_vec(const std::string& key) const;

    /**
     * @brief Retrieve a structured map (returns an empty `StructMap` if absent).
     * @param key  Parameter name.
     * @return     Const reference to the stored `StructMap`, or an empty static map.
     */
    const StructMap& get_struct_map(const std::string& key) const;

    /**
     * @brief Probe for any entry under `key` across all typed tables.
     * @param key  Parameter name.
     * @return     True iff some typed table contains `key`.
     */
    bool has(const std::string& key) const;

    /**
     * @brief Probe specifically for a `StructMap` entry under `key`.
     * @param key  Parameter name.
     * @return     True iff the struct-map table contains `key`.
     */
    bool has_struct_map(const std::string& key) const;

    // Read-only views for serialisers (e.g. params_to_dict in bindings).

    /** @brief Read-only view of the int64 table (for serialisers). */
    const std::unordered_map<std::string, int64_t>& ints() const { return ints_; }
    /** @brief Read-only view of the bool table (for serialisers). */
    const std::unordered_map<std::string, bool>& bools() const { return bools_; }
    /** @brief Read-only view of the string table (for serialisers). */
    const std::unordered_map<std::string, std::string>& strings() const { return strings_; }
    /** @brief Read-only view of the int-vec table (for serialisers). */
    const std::unordered_map<std::string, std::vector<int>>& int_vecs() const { return int_vecs_; }
    /** @brief Read-only view of the uint-vec table (for serialisers). */
    const std::unordered_map<std::string, std::vector<uint64_t>>& uint_vecs() const { return uint_vecs_; }
    /** @brief Read-only view of the struct-map table (for serialisers). */
    const std::unordered_map<std::string, StructMap>& struct_maps() const { return struct_maps_; }

private:
    std::unordered_map<std::string, int64_t> ints_;
    std::unordered_map<std::string, bool> bools_;
    std::unordered_map<std::string, std::string> strings_;
    std::unordered_map<std::string, std::vector<int>> int_vecs_;
    std::unordered_map<std::string, std::vector<uint64_t>> uint_vecs_;
    std::unordered_map<std::string, StructMap> struct_maps_;
};

}  // namespace regpoly::core
