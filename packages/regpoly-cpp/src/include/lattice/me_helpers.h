// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "bitvect.h"
#include "generator.h"
#include "transformation.h"
#include <vector>
#include <memory>

/**
 * @file me_helpers.h
 * @brief Shared helpers for the lattice-family equidistribution methods.
 * @ingroup core
 *
 * Common result type and building blocks consumed by `lattice` /
 * `harase` / `notprimitive` / `simd_notprimitive`: combined
 * characteristic polynomial computation (`polychar_comb`), raw
 * generating-polynomial assembly (`find_polys`), polynomial
 * normalisation (`normalize_polys`), and the canonical Couture-
 * L'Ecuyer dual-lattice test (`test_me_lat`).
 */

namespace regpoly::core {

/**
 * @brief Result of the lattice-family equidistribution kernels.
 *
 * @ingroup core
 */
struct MeLatResult {
    std::vector<int> ecart;  ///< Per-resolution gap; indexed `0..maxL`; `ecart[0]` unused.
    int se;                  ///< Cumulative equidistribution gap.
};

/**
 * @brief Compute the product of the individual characteristic polynomials.
 *
 * @param gens  Component generators.
 * @return      `BitVect` of `(K_total + 1)` bits where bit `i` = coefficient of `z^i`.
 */
BitVect polychar_comb(const std::vector<Generator*>& gens);

/**
 * @brief Build the raw generating polynomials `g_i(z)` for `i = 0..resolution-1`.
 *
 * Each polynomial is stored in `polys[i]` as a `BitVect` of `(K + 1)`
 * bits, bit `j` = coefficient of `z^j`.
 *
 * @param gens        Component generators.
 * @param trans       Per-component tempering chains.
 * @param K           Combined state size.
 * @param M           Combined characteristic polynomial.
 * @param polys       Output polynomials (resized by the routine).
 * @param resolution  Number of polynomials to compute.
 */
void find_polys(
    const std::vector<Generator*>& gens,
    const std::vector<std::vector<Transformation*>>& trans,
    int K, const BitVect& M,
    std::vector<BitVect>& polys, int resolution);

/**
 * @brief Normalise the polynomials by multiplying each by `g_1^{-1} mod M(z)`.
 *
 * If `g_0` shares a common factor with `M`, divides it out from `M`
 * and all polys. Updates `M` and returns the effective degree (may
 * be less than `K`).
 *
 * @param polys       Polynomials to normalise in place.
 * @param M           Combined characteristic polynomial (mutated).
 * @param K           Combined state size.
 * @param resolution  Number of polynomials.
 * @return            Effective degree after any GCD reduction.
 */
int normalize_polys(
    std::vector<BitVect>& polys, BitVect& M, int K, int resolution);

/**
 * @brief Couture-L'Ecuyer dual-lattice equidistribution test (vector form).
 *
 * Selected via `MethodRegistry` by the name string `"lattice"`.
 *
 * @param gens   Component generators.
 * @param trans  Per-component tempering chains.
 * @param kg     Combined state size.
 * @param L      Output word width.
 * @param maxL   Maximum resolution to test.
 * @param delta  Per-resolution gap budget (size `maxL + 1`).
 * @param mse    Upper bound on the cumulative gap.
 * @return       Per-resolution `ecart` and cumulative `se`.
 */
MeLatResult test_me_lat(
    const std::vector<Generator*>& gens,
    const std::vector<std::vector<Transformation*>>& trans,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse);

/**
 * @brief Couture-L'Ecuyer dual-lattice test (single-`Generator&` overload).
 *
 * Phase 1 forwarding adapter — if `gen` is a `CombinedGenerator`,
 * its components and per-component tempering chains are unpacked
 * and dispatched to the vector overload. Any other `Generator` is
 * treated as a 1-component combination with no tempering. The new
 * overload is the public API; future phases migrate the underlying
 * implementation to consume `Generator&` natively.
 *
 * @param gen    Generator (typically a `CombinedGenerator`).
 * @param kg     Combined state size.
 * @param L      Output word width.
 * @param maxL   Maximum resolution to test.
 * @param delta  Per-resolution gap budget (size `maxL + 1`).
 * @param mse    Upper bound on the cumulative gap.
 * @return       Per-resolution `ecart` and cumulative `se`.
 */
MeLatResult test_me_lat(
    const Generator& gen,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse);

}  // namespace regpoly::core
