// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include <vector>

/**
 * @file resolution_sets.h
 * @brief Psi-12 / Phi-4 resolution-set helpers.
 * @ingroup core
 *
 * Resolution-set helpers used by the equidistribution and
 * collision-free test orchestration. Phase 2 migrates these from
 * `regpoly.analyses` (Python) to C++ so a C++-only user can drive
 * the same tests without going through the Python wrapper.
 *
 * Convention: both functions return a `vector<bool>` indexed by
 * resolution (size `kg + 1` for `Phi_4`, size `L + 1` for `Psi_12`).
 * `out[r]` is true iff resolution `r` belongs to the set.
 */

namespace regpoly::core {

/**
 * @brief `Psi_12`: the resolutions l in {1..L} whose
 *        dimension-equidistribution gap must be tested.
 *
 * Mirrors `regpoly.analyses.equidistribution_test._compute_psi12`:
 *  1. Mark `l = 1..min(isqrt(kg), L)`.
 *  2. For `t = max(2, kg/L) .. isqrt(kg-1)`, mark `min(kg/t, L)`.
 *
 * @param kg  Combined state size.
 * @param L   Output word width.
 * @return    Vector of size `L + 1`; `out[0]` is always false;
 *            `out[l]` is true iff `l` is in `Psi_12`.
 */
std::vector<bool> compute_psi12(int kg, int L);

/**
 * @brief `Phi_4`: the dimensions t in {2..kg} whose collision-free
 *        rank must be checked.
 *
 * Mirrors `regpoly.analyses.collision_free_test._compute_phi4`: for
 * each `t in {2..kg}` where `kg % t != 0`, let `lt = kg / t` and
 * require `lt < L && kg % (lt + 1) != 0 && kg / (t - 1) > lt`.
 *
 * @param kg  Combined state size.
 * @param L   Output word width.
 * @return    Vector of size `kg + 1`; `out[0] = out[1] = false`.
 */
std::vector<bool> compute_phi4(int kg, int L);

}  // namespace regpoly::core
