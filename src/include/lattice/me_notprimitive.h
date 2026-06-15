// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "bitvect.h"
#include "combined_f2_linear_source.h"
#include "generator.h"
#include "transformation.h"
#include "me_helpers.h"  // EquidistributionResult
#include <vector>

/**
 * @file me_notprimitive.h
 * @brief Equidistribution test for generators with non-primitive characteristic polynomial.
 * @ingroup core
 *
 * Equidistribution test that does NOT assume the combined generator
 * is full-period. Selected via `EquidistributionMethodRegistry` by
 * the name `"notprimitive"`.
 *
 * Pipeline (per `docs/theory/equidistribution-spec.md`):
 *  1. Recover `chi_f` via Berlekamp-Massey on a scalar functional of
 *     the combined output.
 *  2. Factor `chi_f` over `F_2` (NTL `CanZass` on `GF2X`).
 *  3. Pick the largest-degree irreducible factor `phi` with
 *     certifiable maximal period.
 *  4. Build a basis `B` of the invariant subspace `V = Ker phi(f)`.
 *  5. Run the matricial DE core on `V`: maintain `p` virtual-register
 *     clones, step them in lockstep, read `L`-bit outputs, insert
 *     dual rows into an `F_2` echelon and read off `k(v)`.
 */

namespace regpoly::core {

/**
 * @brief Run the not-primitive equidistribution test.
 *
 * @param cs     Combined source. Throws `std::invalid_argument` if
 *               any component is not a Recurrence (the Krylov
 *               χ-recovery + per-component projection step both need
 *               recurrence-state evolution).
 * @param kg     Combined state size.
 * @param L      Output word width.
 * @param maxL   Maximum resolution to test.
 * @param delta  Per-resolution gap budget (size `maxL + 1`).
 * @param mse    Upper bound on the cumulative gap.
 * @return       Per-resolution `ecart` and cumulative `se`.
 */
EquidistributionResult test_me_notprimitive(
    const CombinedF2LinearSource& cs,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse);

}  // namespace regpoly::core
