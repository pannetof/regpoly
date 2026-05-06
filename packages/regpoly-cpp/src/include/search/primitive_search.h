// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "params.h"
#include "search_types.h"

#include <cstdint>
#include <functional>
#include <string>

// Phase 2.4: full-period (primitive char-poly) search driver.
//
// Replaces the inner loop of
// packages/regpoly/src/regpoly/search/search_primitive.py:
//
//   while not (max_tries reached or max_seconds reached):
//       full_params = randomize_missing(structural + fixed, family_specs)
//       gen = create_generator(family, full_params, L)        # may throw
//       if gen.is_full_period():
//           emit TestedGenerator{family, gen.k(), L, full_params}
//       periodically emit SearchProgress{tries, elapsed_seconds}
//
// Python remains responsible for: YAML parsing, .partial recovery
// I/O, dedup against the existing-results set, the merged output
// file, header/footer printing. The driver owns the random sampling
// (via regpoly_random::sample_param_into), the factory call, and the
// is_full_period check.

struct PrimitiveSearchConfig {
    std::string family;
    int L;
    Params structural_params;  // user-pinned structural inputs (define k)
    Params fixed_params;       // user-pinned non-structural inputs
    int64_t max_tries;         // 0 = no try-count limit
    double max_seconds;        // 0.0 = no wall-time limit
    int progress_interval;     // emit on_progress every N tries; default 100
    uint64_t random_seed;      // 0 = derive from OS entropy
};

using OnHitFn = std::function<void(const TestedGenerator&)>;
using OnProgressFn = std::function<void(const SearchProgress&)>;

// Returns the total number of tries actually executed.
int64_t run_primitive_search(
    const PrimitiveSearchConfig& cfg,
    const OnHitFn& on_hit,
    const OnProgressFn& on_progress);
