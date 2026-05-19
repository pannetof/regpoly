// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "bitvect.h"
#include "generator.h"

/**
 * @file bm.h
 * @brief Packed Berlekamp-Massey over F_2.
 * @ingroup core
 *
 * Compute the minimal polynomial of an F_2 sequence via
 * Berlekamp-Massey. Runs the generator for `2 * K` steps after
 * `init(init_state)`, observing one output bit per step.
 */

namespace regpoly::core {

/**
 * @brief Packed Berlekamp-Massey on a `Generator`'s output bit stream.
 *
 * Runs `gen` (after `init(init_state)`) for `2 * K` steps, observing
 * bit `bit_idx` of `get_output()` at each step. Returns the linear
 * complexity `L` of the resulting binary sequence. If `out_min_poly`
 * is non-null, writes the recovered minimal polynomial into a
 * freshly-sized `BitVect` of `K` bits, MSB-first: bit `j` holds the
 * coefficient of `z^(K - j)`. (Same convention as
 * `Generator::char_poly()` returns.)
 *
 * `K` should be an upper bound on the actual minimal polynomial
 * degree (typically `gen.k()` — but for non-full-period generators
 * the true minimal polynomial may be smaller).
 *
 * @param gen           Generator to probe.
 * @param init_state    Initial state for `gen.init`.
 * @param K             Upper bound on the minimal polynomial degree.
 * @param out_min_poly  Output minimal polynomial (`nullptr` to skip).
 * @param bit_idx       Which output bit to feed Berlekamp-Massey (default `0`).
 * @return              Linear complexity of the observed sequence.
 */
int packed_bm(const Generator& gen,
              const BitVect& init_state,
              int K,
              BitVect* out_min_poly,
              int bit_idx = 0);

}  // namespace regpoly::core
