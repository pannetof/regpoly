// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "bitvect.h"
#include "generator.h"

#include <optional>
#include <string>
#include <vector>

/**
 * @file primitivity.h
 * @brief Primitivity testing for characteristic polynomials over F_2.
 * @ingroup core
 *
 * Phase 2.4: primitivity testing for characteristic polynomials.
 *
 * Mirrors the algorithm previously in
 * `packages/regpoly/src/regpoly/search/primitivity.py`:
 *
 *   - For Mersenne-prime exponents `k`, irreducibility implies
 *     primitivity (since `2^k - 1` is itself prime).
 *   - Otherwise, every prime factor `p` of `2^k - 1` must satisfy
 *     `x^((2^k - 1)/p) != 1 (mod char_poly)`.
 *
 * Prime factors of `Phi_n(2)` for every divisor `n` of `k` are
 * sourced from the embedded Cunningham-style table (see
 * `src/algebra/primitive_factors_data.cpp`, generated from
 * `packages/regpoly/src/regpoly/data/primitive_factors.json`).
 *
 * @see :py:mod:`regpoly.search.primitivity`
 */

namespace regpoly::core {

/**
 * @brief True iff `k` is a Mersenne prime exponent.
 *
 * @param k  Exponent to test.
 * @return   True iff `2^k - 1` is itself a (Mersenne) prime.
 */
bool is_mersenne_prime_exponent(int k);

/**
 * @brief Sorted prime factors of `2^k - 1` as decimal strings.
 *
 * Factors can exceed 64 bits, so they are returned as decimal
 * strings rather than integers.
 *
 * @param k  Exponent.
 * @return   Sorted prime factors, or `std::nullopt` when the
 *           factorisation cannot be assembled (a divisor's entry is
 *           missing or flagged as incomplete in the table).
 */
std::optional<std::vector<std::string>> get_primitive_factors_for_k(int k);

/**
 * @brief Check whether a generator's characteristic polynomial is primitive.
 *
 * Equivalent to asking whether the generator has full period
 * `2^k - 1` with `k = gen.k()`.
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   if (is_full_period(gen)) {
 *       // gen.char_poly() is primitive over GF(2).
 *   }
 * @endcode
 *
 * @param gen  Generator to test.
 * @return     True iff the characteristic polynomial is primitive.
 * @throws std::runtime_error  When the factorisation of `2^k - 1` is unavailable.
 */
bool is_full_period(const Generator& gen);

/**
 * @brief Primitivity check on a precomputed characteristic polynomial.
 *
 * Convenience overload for callers that already have `char_poly` and
 * don't want to recompute it.
 *
 * @param char_poly  Characteristic polynomial.
 * @param k          Polynomial degree.
 * @return           True iff the polynomial is primitive over GF(2).
 */
bool is_full_period(const BitVect& char_poly, int k);

/**
 * @brief Test whether `z^w + lower_w_bits` is irreducible over GF(2).
 *
 * Here `lower_w_bits = modM_value & ((1 << w) - 1)`. Used by the
 * `irreducible_gf2` random-parameter sampler to validate a candidate
 * `modM` for the F2w generator family (which packs `modM` LSB-first
 * in its lower `w` bits, with the leading `z^w` term implicit).
 *
 * @param modM_value  Candidate modulus packed LSB-first.
 * @param w           Field width.
 * @return            True iff the polynomial is irreducible over GF(2).
 */
bool is_irreducible_gf2w_modM(uint64_t modM_value, int w);

}  // namespace regpoly::core

/**
 * @brief Internal data-access helpers (exposed for GoogleTest fixtures).
 *
 * Lives under `regpoly::internal` rather than `regpoly::core` to
 * mark these as implementation details that may be reshaped without
 * notice.
 */
namespace regpoly::internal {

/**
 * @brief Look up the prime factors of `Phi_n(2)` for a single `n`.
 *
 * @param n         Cyclotomic index.
 * @param complete  Output: set to the table entry's `"complete"` flag.
 * @return          Pointer to the sorted factor list, or `nullptr` if
 *                  `n` is not in the table.
 */
const std::vector<std::string>* lookup_factors(int n, bool& complete);

}  // namespace regpoly::internal
