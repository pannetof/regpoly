// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "bitvect.h"
#include "generator.h"
#include "transformation.h"
#include "me_helpers.h"   // MeLatResult
#include <vector>

/**
 * @file me_notprimitive_simd.h
 * @brief SIMD-aware not-primitive equidistribution test (Saito-Matsumoto 2008 §3.2).
 * @ingroup core
 *
 * SIMD-aware non-primitive equidistribution test. Selected via
 * `MethodRegistry` by the name string `"simd_notprimitive"`.
 *
 * Pipeline:
 *  1-3. Recover `chi_f`, factor it, select the largest Mersenne
 *       primitive `phi` (shared with `test_me_notprimitive` via
 *       `chi_recovery.h`).
 *  4.   Compute `V`-projected seeds via `apply_polynomial(psi, .)`.
 *
 * Fast path: if `view_proto.simd_lane_count() == 1` AND `p == kg`,
 * dispatch to `test_me_lat`. Otherwise fall through.
 *
 *  5. For each `(sm, wm)` in `{0..elementNo - 1} x {1..elementNo}`:
 *       - advance generator clones by `sm` `next()` calls
 *         (= `start_mode`);
 *       - run a SIMD-PIS reduction (port of Saito's
 *         `AlgorithmSIMDEquidistribution.hpp`) with the
 *         `V`-projected seed as the PIS "rand" vector, where each
 *         PIS step reads `elementNo` consecutive `next()` outputs
 *         and assembles them as a 128-bit super-word with
 *         `weight_mode` blending against a cached `previous`
 *         super-word, then extracts `v` MSBs of each lane
 *         interleaved into the basis "next" `BitVect`;
 *       - rescale `e = min_count * elementNo - (elementNo - wm)`.
 *     Aggregate per-`sm` MAX over `wm`, then MIN over `sm`.
 *  6. Emit `ecart[v] = floor(p/v) - k(v)`.
 *
 * Generic across generators via the abstract
 * `Generator::simd_lane_count()` hint. For `simd_lane_count() == 1`,
 * the algorithm collapses to behaviour identical to plain
 * `test_me_notprimitive`.
 */

namespace regpoly::core {

/**
 * @brief Run the SIMD-aware not-primitive equidistribution test (vector form).
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
MeLatResult test_me_notprimitive_simd(
    const std::vector<Generator*>& gens,
    const std::vector<std::vector<Transformation*>>& trans,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse);

/**
 * @brief Run the SIMD-aware not-primitive test (single-`Generator&` overload).
 *
 * Phase 1 forwarding adapter — see `me_helpers.h` for the unpacking
 * convention.
 *
 * @param gen    Generator (typically a `CombinedGenerator`).
 * @param kg     Combined state size.
 * @param L      Output word width.
 * @param maxL   Maximum resolution to test.
 * @param delta  Per-resolution gap budget (size `maxL + 1`).
 * @param mse    Upper bound on the cumulative gap.
 * @return       Per-resolution `ecart` and cumulative `se`.
 */
MeLatResult test_me_notprimitive_simd(
    const Generator& gen,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse);

}  // namespace regpoly::core
