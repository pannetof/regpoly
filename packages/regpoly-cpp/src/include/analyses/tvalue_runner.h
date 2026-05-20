// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.

#pragma once
#include "generator.h"

#include <stdexcept>
#include <string>
#include <vector>

/**
 * @file tvalue_runner.h
 * @brief Top-level t-value computation for F_2 digital nets.
 * @ingroup core
 *
 * Computes the profile `t(s)` over `s = 2..s_max` for the (t,m,s)-net
 * implied by a `Generator` (typically a `DigitalNet`, but any
 * `Generator` works — the kernel sees only the F_2-linear matrix that
 * `GaussMatrix::prepare` builds from `copy()` + `init` + `next` +
 * `get_output`).
 *
 * Definition (Niederreiter 1987; Dick–Pillichshammer 2010 §4.2):
 *   t(s) = m − d_max(s)
 * where d_max(s) is the largest non-negative integer d such that
 * some composition (l_1, …, l_s) of d with `l_j ∈ [0, m]` yields
 * full rank d in the stacked top-l_j rows of the s generating
 * matrices C_1, …, C_s.
 *
 * **v1 algorithm**: naive primal Schmid-style enumeration of every
 * composition with full Gaussian elimination on each — correct but
 * O(C(d+s−1, s−1) · m^3) per s. Designed for the small-`(m, s)`
 * regime described in the plan (m ≤ 30, s ≤ 8). The Niederreiter–
 * Pirsic dual method registered alongside throws — Phase-4 sets up
 * the dispatch site, the dual implementation arrives later.
 */

namespace regpoly::core {

/**
 * @brief Result of `run_tvalue_profile_*`.
 *
 * @ingroup core
 */
struct TValueResult {
    std::vector<int> tvals;   ///< Per-dimension t-value; size `s_max + 1`;
                              ///< `tvals[0]` unused, `tvals[1] = 0` by convention.
    int se;                   ///< Σ_{s=2..s_max} tvals[s] (plain sum, mirrors
                              ///< equidistribution `se`).
    bool verified;            ///< True iff every dimension finished without
                              ///< tripping the per-s or aggregate cap.
};

/**
 * @brief Schmid-style primal t-value profile over s = 2..s_max.
 *
 * Uses `GaussMatrix::prepare` to build the master matrix once, then
 * enumerates compositions for each `s`. Short-circuits on the same
 * pattern as `run_matricial_equidistribution`: bail when any
 * `tvals[s] > delta[s]` or running `se > max_t_sum`.
 *
 * @param gen        Generator (typically a `DigitalNet` or a
 *                   `CombinedGenerator` wrapping one).
 * @param kg         Combined state size.
 * @param m          Per-coordinate output width (must equal `gen.L()`
 *                   passed to `GaussMatrix::prepare`).
 * @param s_max      Number of dimensions to profile (must be ≥ 2 and
 *                   ≤ the underlying net's coordinate count when
 *                   applicable).
 * @param delta      Per-s cap on t(s), indexed 0..s_max; `delta[s] = INT_MAX`
 *                   for "no cap".
 * @param max_t_sum  Cap on cumulative `se`; INT_MAX for no cap.
 * @return           Per-s `tvals`, cumulative `se`, and `verified` flag.
 */
TValueResult run_tvalue_profile_schmid(
    const Generator& gen,
    int kg, int m, int s_max,
    const std::vector<int>& delta, int max_t_sum);

/**
 * @brief Niederreiter–Pirsic dual t-value profile (NOT YET IMPLEMENTED).
 *
 * Registered for dispatch by name (matches the precedent of
 * `METHOD_DUALLATTICE` in `equidistribution_test.py`). Throws
 * `std::runtime_error` when called.
 */
TValueResult run_tvalue_profile_dual(
    const Generator& gen,
    int kg, int m, int s_max,
    const std::vector<int>& delta, int max_t_sum);

}  // namespace regpoly::core
