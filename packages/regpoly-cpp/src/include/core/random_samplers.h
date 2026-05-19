// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "param_spec.h"
#include "params.h"

#include <cstdint>
#include <random>
#include <string>
#include <vector>

/**
 * @file random_samplers.h
 * @brief Generic random parameter samplers in C++ (Phase 2.4).
 *
 * Mirrors the four "generic" `rand_type` cases of
 * `packages/regpoly/src/regpoly/core/parametric.py:generate_random`:
 *
 *  - `bitmask`         → `random_bits(eval_expr(rand_args, params))`
 *  - `range`           → `uniform_int_distribution(lo, hi)` where
 *                        `rand_args` is `"lo_expr,hi_expr"`
 *  - `poly_exponents`  → sample `n` distinct exponents in `[1, k-1]`,
 *                        return sorted list with trailing 0 appended;
 *                        `k` from `eval_expr(rand_args, params)`;
 *                        `n` drawn uniformly in `[1, min(k-1, 10)]`
 *  - `bitmask_vec`     → vector of length `params[length_param]` of
 *                        random bitmask draws of width
 *                        `eval_expr(bits_param, params)`; `rand_args`
 *                        is `"bits_param,length_param"`
 *
 * Family-specific `rand_type`s (anything not in the four above) are
 * dispatched separately by the existing `random_param` binding.
 *
 * These functions deliberately do NOT call back into Python's
 * `random` module: they own a thread-local `std::mt19937_64` so the
 * C++ search drivers can run free of the GIL. The produced sequence
 * is intentionally distinct from Python's MT for the same seed.
 *
 * @ingroup core
 */

namespace regpoly::random {

using regpoly::core::Params;
using regpoly::core::ParamSpec;

/**
 * @brief Seed the thread-local RNG.
 *
 * Pass `0` to derive a seed from the OS entropy source. Matches
 * Python's `random.seed(int)` only for the purpose of being
 * deterministic — the produced sequence WILL differ from Python's
 * MT for the same seed.
 *
 * @param s  Seed; `0` → entropy-derived.
 */
void seed(uint64_t s);

/**
 * @brief Direct access to the thread-local RNG.
 *
 * For callers that want to integrate with their own samplers
 * (e.g. the search drivers use this for non-spec-driven choices
 * like which generator-list element to pick).
 *
 * @return  Reference to the thread-local `std::mt19937_64`.
 */
std::mt19937_64& rng();

/**
 * @brief Evaluate a parameter-spec expression.
 *
 * Accepts a literal int, a parameter name, or `"name+N"` /
 * `"name-N"` with `N` a non-negative literal int. `params` must
 * contain the referenced parameter as an int.
 *
 * @param expr    Parameter-spec expression.
 * @param params  Typed parameter bag the expression resolves against.
 * @return        Evaluated value.
 * @throws std::invalid_argument  If `expr` is malformed or the referenced
 *                                parameter is missing from `params`.
 */
int64_t eval_param_expr(const std::string& expr, const Params& params);

/**
 * @brief Sample a value for one of the four generic `rand_type`s and write it into `params`.
 *
 * Returns `false` if the `rand_type` is not one of the four generic
 * types (caller must then dispatch to the family-specific
 * `random_param`).
 *
 * @param spec    Declarative spec describing the parameter to fill.
 * @param params  Bag to update under `spec.name`.
 * @return        `true` on success; `false` if `spec.rand_type` is not generic.
 * @throws std::invalid_argument  On malformed `rand_args` (e.g. missing comma,
 *                                hi < lo, missing referenced parameter).
 * @throws std::runtime_error     For `irreducible_gf2`, if the retry budget
 *                                is exhausted without finding an irreducible.
 */
bool sample_generic_into(const ParamSpec& spec, Params& params);

/**
 * @brief Unified front-end: tries generic samplers then the family-specific dispatch.
 *
 * Falls back to the family-specific dispatch
 * (`TauswortheGen::generate_random` today; future families register
 * here too). Returns `false` only when no sampler exists for the
 * spec's `rand_type`.
 *
 * `L` is needed only because the family-specific dispatch path
 * consumes it; it is included here so the caller can build a
 * uniform `sample_param(spec, params, L)` front-end.
 *
 * @param spec    Declarative spec describing the parameter to fill.
 * @param params  Bag to update under `spec.name`.
 * @param L       Output word width forwarded to family-specific dispatch.
 * @return        `true` on success; `false` if no sampler matches.
 * @throws std::invalid_argument  Propagated from the underlying sampler.
 * @throws std::runtime_error     Propagated from the underlying sampler.
 *
 * @code{.cpp}
 *   using namespace regpoly::random;
 *   seed(0);                              // 0 → OS entropy
 *   Params p;
 *   p.set_int("w", 32);
 *   ParamSpec spec{ "a", "int", false, false, 0, "bitmask", "w", false };
 *   sample_param_into(spec, p, 32);       // fills p["a"] with a random w-bit mask
 * @endcode
 */
bool sample_param_into(const ParamSpec& spec, Params& params, int L);

}  // namespace regpoly::random
