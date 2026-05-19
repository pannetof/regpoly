// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "params.h"

#include <cstdint>
#include <string>

/**
 * @file search_types.h
 * @brief Typed result / progress structs shared by every C++ search driver.
 * @ingroup core
 *
 * Phase 2.4: typed structs shared by every C++ search driver.
 *
 * `SearchProgress` is emitted periodically by the driver (every
 * `progress_interval` tries by default) and at the very end. It is
 * intentionally narrow: callers that need richer state (best-so-far
 * SE for tempering, found count for primitive search) should track
 * it themselves in the `on_progress` closure.
 *
 * `TestedGenerator` is what the driver hands to the `on_hit` callback
 * for each retained candidate. It carries the params that built the
 * generator plus the structural metadata (family, k, L) the caller
 * needs to persist or display the result.
 */

namespace regpoly::core {

/**
 * @brief Periodic progress snapshot emitted by every search driver.
 *
 * Drivers fire this through their `on_progress` callback at a fixed
 * iteration cadence (default: every 100 tries) and once at the very
 * end of the loop. Keep the payload narrow — callers track richer
 * state (best-so-far SE, found-count, etc.) in their own closures.
 *
 * @ingroup core
 */
struct SearchProgress {
    int64_t tries;            ///< Number of tries executed so far (1-based).
    double elapsed_seconds;   ///< Wall-clock time since the driver started.
};

/**
 * @brief One accepted candidate handed to a primitive-search `on_hit` callback.
 *
 * Carries the parameters that materialise a runtime generator plus
 * the structural metadata (`family`, `k`, `L`) the caller needs to
 * persist or display the result without re-running the factory.
 *
 * @ingroup core
 */
struct TestedGenerator {
    std::string family;       ///< Canonical generator family name.
    int k;                    ///< State width in bits.
    int L;                    ///< Output word width in bits.
    Params params;            ///< Params used to construct the generator.
    /**
     * @brief Iteration index (1-based) at which the driver accepted this
     *        candidate. `0` when not produced by an iteration loop.
     */
    int64_t tries_at_hit;
};

}  // namespace regpoly::core
