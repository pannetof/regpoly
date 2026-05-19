// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "bitvect.h"
#include "generator.h"
#include "transformation.h"
#include "dual_lattice.h"
#include "me_helpers.h"
#include <NTL/GF2X.h>
#include <vector>
#include <memory>

/**
 * @file temper_optimizer.h
 * @brief StackBase-cached dual lattice used by the tempering optimiser.
 * @ingroup core
 *
 * `TemperOptCache` — dual lattice with StackBase for tempering
 * optimisation.
 *
 * Mirrors the C code's strategy in `temperMK.c`:
 *  - `compute_all()`: run Lenstra for `v = 1..L`, caching the basis
 *    at each `v`.
 *  - `compute_gap(v)`: recompute ONLY resolution `v` after a bitmask
 *    change:
 *      1. Restore the cached basis at `v - 1`.
 *      2. Recompute the generating polynomial `g_v` for the current
 *         bitmask.
 *      3. Normalise: `h_v = g_v * g0_inv mod M`.
 *      4. `DualBaseIncrease` + Lenstra at resolution `v`.
 *      5. Return the gap.
 *
 * Cost per `compute_gap(v)`: O(k) for `FindPolys` at one resolution
 * plus one Lenstra.
 */

namespace regpoly::core {

/**
 * @brief Stateful dual-lattice cache backing the tempering optimiser.
 *
 * Owns the combined characteristic polynomial, the normalised
 * generating polynomials, the StackBase of dual-lattice snapshots
 * (one per resolution), and the incremental snapshot used by the
 * step-driven optimisation loop. After each bitmask perturbation,
 * `compute_gap(v)` recomputes only resolution `v` by restoring the
 * cached basis at `v - 1`.
 *
 * @ingroup core
 */
class TemperOptCache {
public:
    /**
     * @brief Construct a cache from component generators + tempering chains.
     *
     * @param gens   Component generators.
     * @param trans  Per-component tempering chains.
     * @param kg     Combined state size.
     * @param L      Output word width / maximum resolution.
     */
    TemperOptCache(
        const std::vector<Generator*>& gens,
        const std::vector<std::vector<Transformation*>>& trans,
        int kg, int L);

    /**
     * @brief Single-`Generator&` forwarding constructor.
     *
     * @param gen  Generator (typically a `CombinedGenerator`).
     * @param kg   Combined state size.
     * @param L    Output word width.
     */
    TemperOptCache(const Generator& gen, int kg, int L);

    /**
     * @brief Compute all gaps `v = 1..L` and populate the StackBase cache.
     *
     * @return  `ecart[v]` for `v = 1..L` (`ecart[0]` unused).
     */
    std::vector<int> compute_all();

    /**
     * @brief Recompute the gap at resolution `v` only.
     *
     * Restores the cached basis at `v - 1`, recomputes `g_v` for the
     * current bitmask, and runs Lenstra at `v`. Cost: O(k) for
     * sequence collection plus one Lenstra.
     *
     * @param v  Resolution to recompute.
     * @return   The new gap at resolution `v`.
     */
    int compute_gap(int v);

    /**
     * @brief Recompute `g_0` and `inv_g0` from the current bitmask.
     *
     * Must be called after any perturbation that changes output bit 0.
     */
    void refresh_inv_g0();

    /**
     * @brief Recompute all polys, `inv_g0`, and snapshots from the current bitmask.
     *
     * @return  `ecart[v]` for `v = 1..L`.
     */
    std::vector<int> rebuild();

    // ── Incremental optimization (mirrors temperMK.c OptimizeTemper) ──

    /**
     * @brief Reset the incremental optimisation state.
     *
     * Call once before starting the optimisation loop that drives
     * `step()`.
     */
    void reset_step();

    /**
     * @brief Compute the gap at resolution `v` with incremental StackBase.
     *
     * Behaviour:
     *  - `v == 1`: rebuild `inv_g0` + DualBase from scratch.
     *  - `v > previous_v`: save snapshot, `DualBaseIncrease`,
     *    Lenstra (forward).
     *  - `v <= previous_v`: restore snapshot, `DualBaseIncrease`,
     *    Lenstra (retry).
     *
     * @param v  Resolution to step to.
     * @return   The gap at resolution `v`.
     */
    int step(int v);

    /** @brief Effective degree (after any GCD reduction). */
    int kg() const { return Deg_; }
    /** @brief Maximum resolution. */
    int L() const { return L_; }

private:
    std::vector<Generator*> gens_;
    std::vector<std::vector<Transformation*>> trans_;
    int kg_orig_;  // original degree (before GCD reduction)
    int Deg_;      // effective degree (after GCD reduction)
    int L_;
    int RES_;

    // Combined char poly and normalization
    BitVect M_orig_;  // original M (before GCD reduction)
    BitVect M_;
    NTL::GF2X M_ntl_;
    NTL::GF2X inv_g0_;   // g_0^{-1} mod M

    // All normalized generating polynomials
    std::vector<BitVect> polys_;

    // StackBase for compute_all / compute_gap
    std::vector<std::unique_ptr<DualLatticeBase>> stack_;

    // Incremental state for step()
    std::unique_ptr<DualLatticeBase> current_base_;
    std::vector<std::unique_ptr<DualLatticeBase>> opt_stack_;
    int previous_v_ = 0;

    // Compute generating polynomial for a SINGLE resolution index
    void find_single_poly(int idx, BitVect& out_poly);

    // Normalize a single polynomial: h = g * inv_g0 mod M
    void normalize_single(BitVect& poly);
};

}  // namespace regpoly::core
