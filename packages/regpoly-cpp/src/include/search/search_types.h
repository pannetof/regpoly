// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "params.h"

#include <cstdint>
#include <string>

// Phase 2.4: typed structs shared by every C++ search driver.
//
// SearchProgress is emitted periodically by the driver (every
// `progress_interval` tries by default) and at the very end. It is
// intentionally narrow: callers that need richer state (best-so-far
// SE for tempering, found count for primitive search) should track
// it themselves in the on_progress closure.
//
// TestedGenerator is what the driver hands to the on_hit callback for
// each retained candidate. It carries the params that built the
// generator plus the structural metadata (family, k, L) the caller
// needs to persist or display the result.

namespace regpoly::core {

struct SearchProgress {
    int64_t tries;
    double elapsed_seconds;
};

struct TestedGenerator {
    std::string family;
    int k;
    int L;
    Params params;
    int64_t tries_at_hit;     // iteration index at which the driver
                              // accepted this generator (1-based).
                              // 0 when not produced by an iteration loop.
};

}  // namespace regpoly::core
