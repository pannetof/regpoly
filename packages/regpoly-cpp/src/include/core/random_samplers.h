// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "param_spec.h"
#include "params.h"

#include <cstdint>
#include <random>
#include <string>
#include <vector>

// Phase 2.4: generic random parameter samplers in C++.
//
// Mirrors the four "generic" rand_type cases of
// packages/regpoly/src/regpoly/core/parametric.py:generate_random:
//
//   bitmask         → random_bits(eval_expr(rand_args, params))
//   range           → uniform_int_distribution(lo, hi) where rand_args
//                     is "lo_expr,hi_expr"
//   poly_exponents  → sample n distinct exponents in [1, k-1], return
//                     sorted list with trailing 0 appended; k from
//                     eval_expr(rand_args, params); n drawn uniformly
//                     in [1, min(k-1, 10)]
//   bitmask_vec     → vector of length params[length_param] of random
//                     bitmask draws of width eval_expr(bits_param,
//                     params); rand_args is "bits_param,length_param"
//
// Family-specific rand_types (anything not in the four above) are
// dispatched separately by the existing `random_param` binding.
//
// These functions deliberately do NOT call back into Python's random
// module: they own a thread-local std::mt19937_64 so the C++ search
// drivers can run free of the GIL.

namespace regpoly_random {

// Seed the thread-local RNG. Pass 0 to derive a seed from the OS
// entropy source. Match Python's `random.seed(int)` only for the
// purposes of being deterministic — the produced sequence WILL differ
// from Python's MT for the same seed.
void seed(uint64_t s);

// Direct access to the RNG for callers that want to integrate with
// their own samplers (e.g. the search drivers will use this for
// non-spec-driven choices like which generator-list element to pick).
std::mt19937_64& rng();

// Evaluate a parameter-spec expression: a literal int, a parameter
// name, or "name+N"/"name-N" with N a non-negative literal int.
// `params` must contain the referenced parameter as an int.
int64_t eval_param_expr(const std::string& expr, const Params& params);

// Sample a random value for a generic rand_type and write it into
// `params` under `spec.name`. Returns true on success; returns false
// if the rand_type is not one of the four generic types (caller must
// then dispatch to the family-specific `random_param`).
//
// `L` is needed only because the family-specific dispatch path
// (which this function does NOT call) consumes it; included here so
// the caller can build a uniform `sample_param(spec, params, L)`
// front-end.
bool sample_generic_into(const ParamSpec& spec, Params& params);

// Unified front-end: tries the generic samplers first, then falls
// back to the family-specific dispatch (TauswortheGen::generate_random
// today; future families register here too). Returns false only when
// no sampler exists for the spec's rand_type.
bool sample_param_into(const ParamSpec& spec, Params& params, int L);

}  // namespace regpoly_random
