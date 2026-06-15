// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "bitvect.h"
#include "combined_f2_linear_source.h"
#include "generator.h"
#include "transformation.h"
#include "me_helpers.h"  // EquidistributionResult
#include <vector>

/**
 * @file me_harase.h
 * @brief Harase-Matsumoto-Saito (2011) primal-lattice equidistribution method.
 * @ingroup core
 *
 * Harase-Matsumoto-Saito (2011) fast lattice reduction for
 * equidistribution. Works with the PRIMAL lattice (not dual). Uses
 * Mulders-Storjohann weak reduction.
 *
 * Selected via `EquidistributionMethodRegistry` by the name `"harase"`.
 *
 * Reference:
 *   Harase, Matsumoto, Saito (2011). "Fast Lattice Reduction for
 *   F_2-Linear Pseudorandom Number Generators." Math. Comp. 80(273).
 */

namespace regpoly::core {

/**
 * @brief Run the Harase primal-lattice equidistribution test.
 *
 * @param cs     Combined source. Throws `std::invalid_argument` if
 *               any component is not a Recurrence.
 * @param kg     Combined state size.
 * @param L      Output word width.
 * @param maxL   Maximum resolution to test.
 * @param delta  Per-resolution gap budget (size `maxL + 1`).
 * @param mse    Upper bound on the cumulative gap.
 * @return       Per-resolution `ecart` and cumulative `se`.
 */
EquidistributionResult test_me_harase(
    const CombinedF2LinearSource& cs,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse);

/**
 * @brief Compute `k(v)` for a single resolution `v` via the PIS method.
 *
 * @param cs  Combined source. Throws if any component is not a Recurrence.
 * @param kg  Combined state size.
 * @param v   Target resolution.
 * @return    The equidistribution dimension `k(v)`.
 */
int compute_kv(const CombinedF2LinearSource& cs, int kg, int v);

// ── PIS basis with StackBase caching for tempering optimization ─────────

/**
 * @brief PIS basis cache for incremental tempering optimisation.
 *
 * Computes all `k(v)` from `v = L` down to `1`, caching the basis at
 * each `v`. After a bitmask perturbation, `restore_and_reduce(v)`
 * restores the cached basis at `v` and re-reduces to get the new
 * `k(v)`.
 *
 * @ingroup core
 */
class HaraseRankCache {
public:
    /**
     * @brief Construct a cache from a CombinedF2LinearSource.
     *
     * @param cs  Combined source. Throws if any component is not a Recurrence.
     * @param kg  Combined state size.
     * @param L   Output word width / maximum resolution.
     */
    HaraseRankCache(const CombinedF2LinearSource& cs, int kg, int L);

    /**
     * @brief Compute every `k(v)` for `v = L..1` and populate the cache.
     *
     * @return  `ecart[v] = kg/v - k(v)` for `v = 1..L`.
     */
    std::vector<int> compute_all();

    /**
     * @brief Recompute `k(v)` at a single resolution after a perturbation.
     */
    int restore_and_reduce(int v);

    /** @brief Combined state size. */
    int kg() const { return kg_; }
    /** @brief Maximum resolution. */
    int L() const { return L_; }

private:
    std::vector<Recurrence*> gens_;
    int kg_;
    int L_;

    struct Impl;
    std::shared_ptr<Impl> impl_;
};

}  // namespace regpoly::core
