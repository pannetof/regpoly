// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.

#pragma once
#include <cstdint>
#include <string>

// Helpers for decoding the legacy 3-mask form of WELL paper M6 back to
// the paper's (q, t, s, a) signature.
//
// Background: the original WELLGen implementation accepted three full
// 32-bit masks per M6 instance (paramsulong[0..2]) instead of the
// paper's parametric (s, t) integers + single 32-bit constant `a`. For
// every published WELL generator (Panneton et al. 2006, Table II),
// those three masks decompose losslessly:
//   pu[0] = a                 (XOR constant)
//   pu[1] = ~(1u << (31 - s)) (the d_s mask: zero bit at MSB-position s)
//   pu[2] =  (1u << (31 - t)) (the test mask: one bit at MSB-position t)
//
// These helpers extract s and t from the masks. They throw if a mask is
// not in single-bit / single-zero-bit form — in that case the data
// fell outside the paper's M6 family and must be reviewed by hand.

namespace regpoly_well {

// Returns s ∈ [0..31] such that mask == ~(1u << (31 - s)) & 0xFFFFFFFFu.
// Throws std::runtime_error if no such s exists. `context` is included
// in the error message for diagnostics.
int decode_d_s_mask(uint64_t mask, const std::string& context);

// Returns t ∈ [0..31] such that mask == (1u << (31 - t)).
// Throws std::runtime_error if no such t exists.
int decode_test_mask(uint64_t mask, const std::string& context);

}  // namespace regpoly_well
