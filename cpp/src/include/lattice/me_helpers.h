// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "bitvect.h"
#include "combined_f2_linear_source.h"
#include "equidistribution_runner.h"
#include "generator.h"
#include "transformation.h"
#include <vector>
#include <memory>

/**
 * @file me_helpers.h
 * @brief Shared helpers for the lattice-family equidistribution methods.
 * @ingroup core
 *
 * Common building blocks consumed by `lattice` / `harase` /
 * `notprimitive` / `simd_notprimitive`: combined characteristic
 * polynomial computation (`polychar_comb`), raw generating-polynomial
 * assembly (`find_polys`), polynomial normalisation
 * (`normalize_polys`), and the canonical Couture-L'Ecuyer dual-lattice
 * test (`test_me_lat`).
 *
 * Every kernel takes `const CombinedF2LinearSource& cs` and internally
 * calls `cs.recurrence_components("kernel_name")` — which throws
 * `std::invalid_argument` if any component is a `DigitalNet`.
 */

namespace regpoly::core {

/**
 * @brief Product of the individual characteristic polynomials.
 *
 * @param cs  Combined source. Throws if any component is not a Recurrence.
 * @return    `BitVect` of `(K_total + 1)` bits where bit `i` =
 *            coefficient of `z^i`.
 */
BitVect polychar_comb(const CombinedF2LinearSource& cs);

/**
 * @brief Build the raw generating polynomials `g_i(z)` for `i = 0..resolution-1`.
 *
 * Each polynomial is stored in `polys[i]` as a `BitVect` of `(K + 1)`
 * bits, bit `j` = coefficient of `z^j`. Tempering chains live on each
 * generator's intrinsic `F2LinearSource::tempering_`, so the kernel
 * reads chain-composed bits via `get_output()` directly.
 *
 * @param cs          Combined source. Throws if any component is not a Recurrence.
 * @param K           Combined state size.
 * @param M           Combined characteristic polynomial.
 * @param polys       Output polynomials (resized by the routine).
 * @param resolution  Number of polynomials to compute.
 */
void find_polys(
    const CombinedF2LinearSource& cs,
    int K, const BitVect& M,
    std::vector<BitVect>& polys, int resolution);

/**
 * @brief Normalise the polynomials by multiplying each by `g_1^{-1} mod M(z)`.
 *
 * If `g_0` shares a common factor with `M`, divides it out from `M`
 * and all polys. Updates `M` and returns the effective degree.
 */
int normalize_polys(
    std::vector<BitVect>& polys, BitVect& M, int K, int resolution);

/**
 * @brief Couture-L'Ecuyer dual-lattice equidistribution test.
 *
 * Selected via `EquidistributionMethodRegistry` by the name `"lattice"`.
 *
 * @param cs     Combined source. Throws if any component is not a Recurrence.
 * @param kg     Combined state size.
 * @param L      Output word width.
 * @param maxL   Maximum resolution to test.
 * @param delta  Per-resolution gap budget (size `maxL + 1`).
 * @param mse    Upper bound on the cumulative gap.
 * @return       Per-resolution `ecart` and cumulative `se`.
 */
EquidistributionResult test_me_lat(
    const CombinedF2LinearSource& cs,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse);

}  // namespace regpoly::core
