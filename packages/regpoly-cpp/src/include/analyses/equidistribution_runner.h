// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "generator.h"
#include <vector>

/**
 * @file equidistribution_runner.h
 * @brief Top-level matricial-equidistribution and collision-free runners.
 * @ingroup core
 *
 * Phase 2.3: equidistribution test orchestration in C++.
 *
 * `run_matricial_equidistribution(gen, kg, L, Lmax, delta, mse)` packages
 * the outer loop that the Python `EquidistributionTest.run()`
 * previously owned for the matricial method:
 *
 *  1. `prepare_mat` builds a `GaussMatrix`.
 *  2. Walk `l = 1..Lmax`, calling `dimension_equid(mat.copy(), kg, l, L)`.
 *     Track `ecart[l] = floor(kg/l) - dim`. Honour the "verif" state
 *     machine (re-test missed resolutions when a non-zero gap
 *     appears at `l`). Bail out as soon as `ecart[l] > delta[l]` or
 *     `se > mse`.
 *  3. Fill remaining unmeasured `ecart` slots with `0` (verified) or
 *     `INT_MAX` (unverified beyond `maxl`) following the Python
 *     convention.
 *
 * The function consumes a single `Generator&` (typically a
 * `CombinedGenerator` per Phase 1) and returns `ecart`, `se`, and
 * `verified`.
 */

namespace regpoly::core {

/**
 * @brief Result of `run_matricial_equidistribution`.
 *
 * @ingroup core
 */
struct MatricialEquidResult {
    std::vector<int> ecart;   ///< Per-resolution gap; size `Lmax + 1`; `ecart[0]` unused.
    int se;                   ///< Cumulative equidistribution gap.
    bool verified;            ///< True iff every resolution was certified.
};

/**
 * @brief Run the matricial equidistribution analysis on a single generator.
 *
 * Used by `regpoly.analyses.pis.analyze_single_generator` (Python)
 * and `regpoly-cli search` (C++) as the canonical entry point for
 * the matricial method.
 *
 * Pass `delta` with every entry set to `INT_MAX` and `mse =
 * INT_MAX` to put the kernel in unrestricted "just measure" mode
 * (nothing is rejected; the caller wants the score).
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   std::vector<int> delta(L + 1, INT_MAX);
 *   auto res = run_matricial_equidistribution(
 *       *combined, comb.k_g(), L, L, delta, INT_MAX);
 *   std::cout << "SE = " << res.se << "  verified = " << res.verified << '\n';
 * @endcode
 *
 * @param gen    Generator to evaluate (typically a `CombinedGenerator`).
 * @param kg     Combined state size.
 * @param L      Output word width.
 * @param Lmax   Maximum resolution to test.
 * @param delta  Per-resolution gap budget; size `Lmax + 1`.
 * @param mse    Upper bound on the cumulative SE.
 * @return       Per-resolution `ecart`, cumulative `se`, and `verified` flag.
 *
 * @see :py:func:`regpoly.analyses.pis.analyze_single_generator`
 */
MatricialEquidResult run_matricial_equidistribution(
    const Generator& gen,
    int kg, int L, int Lmax,
    const std::vector<int>& delta, int mse);

/**
 * @brief Result of `run_collision_free`.
 *
 * @ingroup core
 */
struct CollisionFreeResult {
    std::vector<int> ecart_cf;  ///< Per-dimension rank gap; size `kg + 1`.
    int secf;                   ///< Sum of collision-free gaps.
    bool verified;              ///< True iff every dimension was certified.
};

/**
 * @brief Run the collision-free analysis on a single generator.
 *
 * Phase 2.3 collision-free orchestration in C++. Walks `t = kg..2`
 * (descending), computing rank deficits via `rang_cf` for every `t`
 * in `Phi_4`. Returns `ecart_cf` indexed `0..kg` and the sum of
 * gaps.
 *
 * `L_for_phi4` is the `L` value used to compute `Phi_4`. In the
 * existing Python orchestration this is `me_results.L` if available,
 * otherwise `gen.L()`. The caller passes the value it would have
 * computed.
 *
 * @param gen          Generator to evaluate.
 * @param kg           Combined state size.
 * @param L            Output word width.
 * @param L_for_phi4   `L` value used by `compute_phi4`.
 * @return             Per-dimension `ecart_cf`, sum `secf`, and `verified` flag.
 */
CollisionFreeResult run_collision_free(
    const Generator& gen,
    int kg, int L, int L_for_phi4);

}  // namespace regpoly::core
