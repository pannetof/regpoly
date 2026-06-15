// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "bitvect.h"
#include "combined_f2_linear_source.h"
#include "generator.h"
#include "transformation.h"
#include "me_helpers.h"   // EquidistributionResult
#include <vector>

/**
 * @file me_notprimitive_simd.h
 * @brief SIMD-aware not-primitive equidistribution test (Saito-Matsumoto 2008 §3.2).
 * @ingroup core
 *
 * Selected via `EquidistributionMethodRegistry` by the name `"simd_notprimitive"`.
 * Generic across generators via the abstract `Recurrence::simd_lane_count()`
 * hint. For `simd_lane_count() == 1`, collapses to behaviour identical
 * to plain `test_me_notprimitive`.
 */

namespace regpoly::core {

/**
 * @brief Run the SIMD-aware not-primitive equidistribution test.
 *
 * @param cs     Combined source. Throws `std::invalid_argument` if
 *               any component is not a Recurrence (SIMD path needs
 *               recurrence-state evolution and `simd_*` virtuals).
 * @param kg     Combined state size.
 * @param L      Output word width.
 * @param maxL   Maximum resolution to test.
 * @param delta  Per-resolution gap budget (size `maxL + 1`).
 * @param mse    Upper bound on the cumulative gap.
 * @return       Per-resolution `ecart` and cumulative `se`.
 */
EquidistributionResult test_me_notprimitive_simd(
    const CombinedF2LinearSource& cs,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse);

}  // namespace regpoly::core
