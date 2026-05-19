// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include <string>
#include <vector>
#include <cstdint>

/**
 * @file param_spec.h
 * @brief Declarative spec describing one parameter of a generator or transformation.
 *
 * Each generator / transformation subclass advertises a static
 * `param_specs()` returning a vector of `ParamSpec` records. The
 * specs drive three orthogonal mechanisms: the YAML loader (which
 * looks up parameter names and types), the random parameter sampler
 * in `regpoly::random` (which dispatches on `rand_type` / `rand_args`),
 * and the tempering optimiser (which only mutates entries marked
 * `optimizable`).
 *
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief Declarative spec for one parameter of a generator or transformation.
 *
 * Every subclass advertises a static `param_specs()` returning a
 * vector of these records. Three flags / fields shape how the
 * surrounding tooling treats a given parameter:
 *
 *  - `structural`: the parameter defines `k` (the size of the F_2
 *    state) or some other shape that downstream code depends on.
 *    The user MUST provide it; the random sampler never overwrites
 *    a structural entry.
 *  - `optimizable`: the tempering optimiser may mutate this value
 *    (used for tempering bitmasks `b`, `c`, etc.).
 *  - `rand_type` + `rand_args`: drive the four generic random
 *    sampling modes — `bitmask`, `range`, `poly_exponents`,
 *    `bitmask_vec` — plus `none` / `""` (not randomizable) and the
 *    family-specific dispatch path (e.g. `tausworthe_poly`).
 *
 * @see regpoly::random::sample_param_into
 *
 * @ingroup core
 */
struct ParamSpec {
    std::string name;       ///< Parameter name as it appears in YAML / Params.
    std::string type;       ///< `"int"`, `"bool"`, `"int_vec"`, or `"uint_vec"`.
    bool structural;        ///< True → defines `k`; user must provide; never random.
    bool has_default;       ///< True → `default_val` is meaningful.
    int64_t default_val;    ///< Default; meaningful only when `has_default && type=="int"|"bool"`.
    std::string rand_type;  ///< `"bitmask"`, `"range"`, `"poly_exponents"`, `"bitmask_vec"`, `"none"`, or `""`.
    /**
     * @brief Auxiliary arguments for the random sampler, format depends on `rand_type`.
     *
     * Encoded as a small string the sampler parses:
     *   - `bitmask`        : param name for the bit count (e.g. `"w"`).
     *   - `range`          : `"min,max-expr"` (e.g. `"1,r-1"`).
     *   - `poly_exponents` : param name for the polynomial degree (e.g. `"k"`).
     *   - `bitmask_vec`    : `"bits_param,length_param"` (e.g. `"w,nocoeff"`).
     *   - `none` / `""`    : not randomizable; ignored.
     */
    std::string rand_args;
    bool optimizable;       ///< True → the tempering optimiser may mutate this entry.
};

}  // namespace regpoly::core
