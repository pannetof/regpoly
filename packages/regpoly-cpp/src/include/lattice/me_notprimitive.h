// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "bitvect.h"
#include "generator.h"
#include "transformation.h"
#include "me_helpers.h"  // MeLatResult
#include <vector>

/**
 * @file me_notprimitive.h
 * @brief Equidistribution test for generators with non-primitive characteristic polynomial.
 * @ingroup core
 *
 * Equidistribution test that does NOT assume the combined generator
 * is full-period (i.e. does not assume the characteristic polynomial
 * `chi_f` is primitive of degree `kg`). Selected via
 * `MethodRegistry` by the name string `"notprimitive"`.
 *
 * Pipeline (per `docs/theory/equidistribution-spec.md`):
 *  1. Recover `chi_f` via Berlekamp-Massey on a scalar functional of
 *     the combined output; LCM across several initial states /
 *     functionals.
 *  2. Factor `chi_f` over `F_2` (NTL `CanZass` on `GF2X`).
 *  3. Pick the largest-degree irreducible factor `phi` whose period
 *     `2^p - 1` can be certified maximal (Mersenne fast path; full
 *     Knuth check otherwise; fall back to "highest certifiable" if
 *     intractable). Set `V = Ker phi(f)`, `p = deg phi`.
 *  4. Build a basis `B` of `V` (Route B: pick random `s`, compute
 *     `s' = psi(f) * s` via Horner where `psi = chi_f / phi`,
 *     collect orbit `{f^v * s' : v = 0..p-1}`; Route A fallback if
 *     `rank < p`).
 *  5. Run the matricial DE core on `V`: maintain `p` "virtual
 *     register" clones of the combined generator, each initialised
 *     to one `B[v]`; step them in lockstep; read `L`-bit outputs;
 *     insert dual rows into a single `F_2` echelon (one per
 *     output-bit phase) and read off `k(v)` for `v = 1..maxL`.
 *
 * The reported `k(v)` is on the invariant subspace `V`, so the
 * corresponding equidistribution gap is `ecart[v] = floor(p/v) -
 * k(v)` (NOT `kg/v - k(v)`); the upper bound is `p`, not `kg`.
 * This is the "honest" `k(v)` for any `F_2`-linear generator
 * regardless of full-period status.
 */

namespace regpoly::core {

/**
 * @brief Run the not-primitive equidistribution test (vector form).
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
MeLatResult test_me_notprimitive(
    const std::vector<Generator*>& gens,
    const std::vector<std::vector<Transformation*>>& trans,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse);

/**
 * @brief Run the not-primitive equidistribution test (single-`Generator&` overload).
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
MeLatResult test_me_notprimitive(
    const Generator& gen,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse);

}  // namespace regpoly::core
