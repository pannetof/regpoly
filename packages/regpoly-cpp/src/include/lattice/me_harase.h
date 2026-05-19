// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "bitvect.h"
#include "generator.h"
#include "transformation.h"
#include "me_helpers.h"  // MeLatResult
#include <vector>

/**
 * @file me_harase.h
 * @brief Harase-Matsumoto-Saito (2011) primal-lattice equidistribution method.
 * @ingroup core
 *
 * Harase-Matsumoto-Saito (2011) fast lattice reduction for
 * equidistribution. Works with the PRIMAL lattice (not dual). Uses
 * Mulders-Storjohann weak reduction: simpler inner loop than
 * Lenstra, no linear-system solve. Avoids the polynomial inversion
 * step (no NTL `InvMod` needed for the lattice reduction itself).
 *
 * The generating polynomials are still computed the same way as in
 * the Couture-L'Ecuyer method (run the generator from canonical
 * state for `k` steps, collect output bits, form `g_i(z)`).
 *
 * Selected via `MethodRegistry` by the name string `"harase"`.
 *
 * Reference:
 *   Harase, Matsumoto, Saito (2011). "Fast Lattice Reduction for
 *   F_2-Linear Pseudorandom Number Generators." Math. Comp. 80(273).
 */

namespace regpoly::core {

/**
 * @brief Run the Harase primal-lattice equidistribution test (vector form).
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
MeLatResult test_me_harase(
    const std::vector<Generator*>& gens,
    const std::vector<std::vector<Transformation*>>& trans,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse);

/**
 * @brief Run the Harase primal-lattice equidistribution test (single-`Generator&` form).
 *
 * Phase 1 forwarding adapter — unpacks a `CombinedGenerator`'s
 * components and tempering chains and dispatches to the vector form.
 *
 * @param gen    Generator (typically a `CombinedGenerator`).
 * @param kg     Combined state size.
 * @param L      Output word width.
 * @param maxL   Maximum resolution to test.
 * @param delta  Per-resolution gap budget (size `maxL + 1`).
 * @param mse    Upper bound on the cumulative gap.
 * @return       Per-resolution `ecart` and cumulative `se`.
 */
MeLatResult test_me_harase(
    const Generator& gen,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse);

/**
 * @brief Compute `k(v)` for a single resolution `v` via the PIS method.
 *
 * Returns the equidistribution dimension (not the gap). Cost: O(k/v)
 * generator steps — much cheaper than full `test_me_harase`.
 *
 * @param gens   Component generators.
 * @param trans  Per-component tempering chains.
 * @param kg     Combined state size.
 * @param v      Target resolution.
 * @return       The equidistribution dimension `k(v)`.
 */
int compute_kv(
    const std::vector<Generator*>& gens,
    const std::vector<std::vector<Transformation*>>& trans,
    int kg, int v);

/**
 * @brief Compute `k(v)` for a single resolution `v` (single-`Generator&` overload).
 *
 * @param gen  Generator (typically a `CombinedGenerator`).
 * @param kg   Combined state size.
 * @param v    Target resolution.
 * @return     The equidistribution dimension `k(v)`.
 */
int compute_kv(const Generator& gen, int kg, int v);

// ── PIS basis with StackBase caching for tempering optimization ─────────

/**
 * @brief PIS basis cache for incremental tempering optimisation.
 *
 * Computes all `k(v)` from `v = L` down to `1`, caching the basis at
 * each `v`. After a bitmask perturbation, `restore_and_reduce(v)`
 * restores the cached basis at `v` and re-reduces to get the new
 * `k(v)`. This avoids recomputing resolutions `v+1..L`.
 *
 * @ingroup core
 */
class PISCache {
public:
    /**
     * @brief Construct a cache from component generators + tempering chains.
     *
     * @param gens   Component generators.
     * @param trans  Per-component tempering chains.
     * @param kg     Combined state size.
     * @param L      Output word width / maximum resolution.
     */
    PISCache(const std::vector<Generator*>& gens,
             const std::vector<std::vector<Transformation*>>& trans,
             int kg, int L);

    /**
     * @brief Single-`Generator&` forwarding constructor.
     *
     * @param gen  Generator (typically a `CombinedGenerator`).
     * @param kg   Combined state size.
     * @param L    Output word width.
     */
    PISCache(const Generator& gen, int kg, int L);

    /**
     * @brief Compute every `k(v)` for `v = L..1` and populate the cache.
     *
     * @return  `ecart[v] = kg/v - k(v)` for `v = 1..L`.
     */
    std::vector<int> compute_all();

    /**
     * @brief Recompute `k(v)` at a single resolution after a perturbation.
     *
     * Restores the cached basis at resolution `v`, rebuilds the
     * generator vector for the current bitmask, reduces, and returns
     * the NEW `k(v)`. Resolutions `v+1..L` are unchanged (guaranteed
     * by safe masks). Cost: O(k/v) — same as a single PIS reduction.
     *
     * @param v  Resolution to recompute.
     * @return   The new `k(v)`.
     */
    int restore_and_reduce(int v);

    /** @brief Combined state size. */
    int kg() const { return kg_; }
    /** @brief Maximum resolution. */
    int L() const { return L_; }

private:
    std::vector<Generator*> gens_;
    std::vector<std::vector<Transformation*>> trans_;
    int kg_;
    int L_;

    // The implementation details are in the .cpp file (opaque structs).
    struct Impl;
    std::shared_ptr<Impl> impl_;
};

}  // namespace regpoly::core
