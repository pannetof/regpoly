// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "bitvect.h"
#include "generator.h"

#include <optional>
#include <string>
#include <vector>

// Phase 2.4: primitivity testing for characteristic polynomials.
//
// Mirrors the algorithm previously in
// packages/regpoly/src/regpoly/search/primitivity.py:
//
//   * For Mersenne prime exponents k, irreducibility implies
//     primitivity (since 2^k - 1 is itself prime).
//   * Otherwise, every prime factor p of 2^k - 1 must satisfy
//     x^((2^k - 1)/p) != 1 (mod char_poly).
//
// Prime factors of Phi_n(2) for every divisor n of k are sourced from
// the embedded Cunningham-style table (see
// src/algebra/primitive_factors_data.cpp, generated from
// packages/regpoly/src/regpoly/data/primitive_factors.json).

namespace regpoly::core {

bool is_mersenne_prime_exponent(int k);

// Returns the sorted list of prime factors of 2^k - 1 as decimal
// strings (factors can exceed 64 bits). Returns nullopt when the
// factorisation cannot be assembled (a divisor's entry is missing or
// flagged as incomplete in the table).
std::optional<std::vector<std::string>> get_primitive_factors_for_k(int k);

// True iff the generator's characteristic polynomial is primitive
// (period 2^k - 1 with k = gen.k()). Throws std::runtime_error when
// the factorisation of 2^k - 1 is unavailable.
bool is_full_period(const Generator& gen);

// Convenience: same as above but takes the precomputed char_poly so
// callers that already have it don't recompute. `k` is the polynomial
// degree.
bool is_full_period(const BitVect& char_poly, int k);

// Tests whether `z^w + lower_w_bits` is irreducible over GF(2), where
// `lower_w_bits = modM_value & ((1<<w) - 1)`. Used by the
// `irreducible_gf2` random-parameter sampler to validate a candidate
// `modM` for the F2w generator family (which packs `modM` LSB-first
// in its lower `w` bits, with the leading `z^w` term implicit).
bool is_irreducible_gf2w_modM(uint64_t modM_value, int w);

}  // namespace regpoly::core

// Internal data-access functions (used by the implementation; exposed
// for the GoogleTest fixtures).
namespace regpoly::internal {

// Lookup factors of Phi_n(2) for a single n. `complete` is set to the
// table entry's "complete" flag. Returns nullptr if n is not in the
// table.
const std::vector<std::string>* lookup_factors(int n, bool& complete);

}  // namespace regpoly::internal
