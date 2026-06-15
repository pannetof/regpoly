// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include <cstdint>
#include <vector>

/**
 * @file gf2w_arith.h
 * @brief GF(2^w) polynomial / normal-basis arithmetic primitives.
 * @ingroup core
 *
 * Helpers used by the F2w generator family to multiply elements of
 * GF(2^w) given a degree-`w` irreducible modulus polynomial.
 */

namespace regpoly::core {

/**
 * @brief Multiply `a` by `z` in GF(2^w) modulo `modM`.
 * @param a     Element packed LSB-first in the low `w` bits.
 * @param k     Modulus degree (= `w`).
 * @param modM  Modulus polynomial (lower `w` bits, leading `z^w` implicit).
 * @return      `a * z mod modM`.
 */
uint64_t gf2w_multiply_z(uint64_t a, int k, uint64_t modM);

/**
 * @brief Multiply two GF(2^w) elements via polynomial multiplication.
 * @param a     First element.
 * @param b     Second element.
 * @param w     Field width.
 * @param modM  Modulus polynomial.
 * @return      `a * b mod modM`.
 */
uint64_t gf2w_multiply_poly(uint64_t a, uint64_t b, int w, uint64_t modM);

/**
 * @brief Multiply two GF(2^w) elements in normal-basis form via a precomputed table.
 * @param a      First element.
 * @param b      Second element.
 * @param w      Field width.
 * @param table  Normal-basis multiplication table from `gf2w_make_table`.
 * @return       `a * b` in the normal-basis representation.
 */
uint64_t gf2w_multiply_normal(uint64_t a, uint64_t b, int w,
                               const uint64_t* table);

/**
 * @brief Build the normal-basis multiplication table for `gf2w_multiply_normal`.
 * @param w     Field width.
 * @param modM  Modulus polynomial.
 * @return      Multiplication table consumed by `gf2w_multiply_normal`.
 */
std::vector<uint64_t> gf2w_make_table(int w, uint64_t modM);

}  // namespace regpoly::core
