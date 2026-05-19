// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "generator.h"
#include <vector>

/**
 * @file tuplets_runner.h
 * @brief Top-level tuplets uniformity-test runner.
 * @ingroup core
 *
 * Phase 2.3: tuplets test orchestration in C++.
 *
 * Computes `Delta(t_1, ..., t_d)` for a combined generator's current
 * state space. Mirrors `regpoly.analyses.tuplets_test.TupletsTest.run()`.
 *
 * Conventions match the Python algorithm exactly:
 *  - `tuph` is 1-indexed: `tuph[0]` is unused, `tuph[1..d]` holds
 *    the `t` values.
 *  - `testtype`: `1` = MAX (matches `_MAX_TYPE`), `0` = SUM
 *    (`_SUM_TYPE`).
 *  - DOUBLE_INF stand-in: `std::numeric_limits<double>::max()` is
 *    used in places where the Python code uses `sys.maxsize` for
 *    double values (`firstpart_max`, `secondpart_max`, `DELTA`).
 *    The Python sentinel is the int `sys.maxsize` cast to float
 *    (~9.22e18); the C++ counterpart uses the same numeric value
 *    via the double of `LLONG_MAX`, so result comparisons across
 *    languages stay byte-identical for the typical test parameter
 *    regimes.
 */

namespace regpoly::core {

/// SUM-type tuplets test (`_SUM_TYPE` in the Python wrapper).
constexpr int TUPLETS_TYPE_SUM = 0;
/// MAX-type tuplets test (`_MAX_TYPE` in the Python wrapper).
constexpr int TUPLETS_TYPE_MAX = 1;

/**
 * @brief Result of `run_tuplets`.
 *
 * @ingroup core
 */
struct TupletsRunResult {
    int tupd;                        ///< Tuplets dimension echoed for the caller.
    std::vector<int> tuph;           ///< Tuplets shape parameters (1-indexed; size `tupd + 1`).
    std::vector<double> gap;         ///< Per-shape gap (1-indexed; size `tuph[1] + 1`).
    std::vector<double> DELTA;       ///< Per-dimension delta (1-indexed; size `tupd + 1`).
    std::vector<double> pourcentage; ///< Per-dimension percentage (1-indexed; size `tupd + 1`).
    double firstpart_max;            ///< Max gap on the first part.
    double firstpart_sum;            ///< Sum of gaps on the first part.
    double secondpart_max;           ///< Max gap on the second part.
    double secondpart_sum;           ///< Sum of gaps on the second part.
};

/**
 * @brief Run the tuplets uniformity test on a single generator.
 *
 * @param gen        Generator to evaluate (typically a `CombinedGenerator`).
 * @param kg         Combined state size.
 * @param L          Output word width.
 * @param tupd       Tuplets dimension.
 * @param tuph       Tuplets shape parameters (1-indexed; size `tupd + 1`).
 * @param threshold  Rejection threshold passed to the kernel.
 * @param testtype   `TUPLETS_TYPE_SUM` (0) or `TUPLETS_TYPE_MAX` (1).
 * @return           Per-shape `gap` plus per-dimension `DELTA` /
 *                   `pourcentage` and the first/second-part summary.
 */
TupletsRunResult run_tuplets(
    const Generator& gen,
    int kg, int L,
    int tupd,
    const std::vector<int>& tuph,
    double threshold,
    int testtype);

}  // namespace regpoly::core
