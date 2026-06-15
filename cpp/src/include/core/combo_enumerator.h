// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "combined_f2_linear_source.h"
#include "generator.h"
#include "transformation.h"

#include <cstdint>
#include <memory>
#include <vector>

/**
 * @file combination.h
 * @brief Cartesian-product iterator over per-slot source pools.
 * @defgroup core Core
 *
 * Defines the core types used by the search loop to enumerate
 * candidate combined generators: `F2LinearSourcePool` (one slot's
 * owned or shared pool of candidate sources plus a per-slot tempering
 * chain) and `ComboEnumerator` (a stateful iterator over the
 * cartesian product of J such slots).
 *
 * The iterator enforces two semantic constraints that mirror the
 * original Python implementation in `regpoly.core.{component,combination}`:
 *
 *  1. **Identity-uniqueness.** The same `Recurrence` object never
 *     appears twice in a single combo — XORing a component with
 *     itself produces the zero sequence, which is never useful.
 *  2. **Shared-pool C(n,k) selection.** When two pool slots
 *     reference the same underlying pool object (set up via
 *     `F2LinearSourcePool::copy_pool_from`), the indices in those
 *     slots must be strictly increasing. A pool of size `n` shared
 *     across `k` slots therefore yields `C(n,k)` combos instead of
 *     `n^k` permutations.
 *
 * On each successful `reset()` / `next()`, `ComboEnumerator` recomputes
 * `k_g` (sum of active generator `k`'s) and `L` (min active `L`,
 * capped at `Lmax`) so callers can read the active combo's effective
 * state size and output width without re-walking the slots.
 *
 * Phase 5.2 rename (2026-05): `Component` → `F2LinearSourcePool`;
 * accessor methods now say `source` instead of `gen`. Storage type
 * is still `unique_ptr<Recurrence>` (`Recurrence`) for now — widening
 * to `unique_ptr<F2LinearSource>` is a follow-up that requires
 * teaching `build_combined_from_enumerator` to choose between
 * `CombinedF2LinearSource` and `CombinedF2LinearSource` depending on
 * whether any pooled source is a `DigitalNet`.
 *
 * @see :py:class:`regpoly.core.combination.ComboEnumerator`
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief One slot of a `ComboEnumerator`: a pool of candidate sources plus a tempering chain.
 *
 * An `F2LinearSourcePool` owns (or shares) a pool of candidate
 * `Recurrence` (currently `Recurrence`; future widening to
 * `F2LinearSource`) instances and an ordered list of `Transformation`s
 * applied on top of the active source's output by the search loop.
 * Two ownership modes are supported:
 *
 *  - **Own pool.** Default after construction. Use `add_source()` to
 *    append deep copies of candidate sources.
 *  - **Shared pool.** After `copy_pool_from(other)`, this pool
 *    references the same underlying `SourcePool` (via `shared_ptr`) as
 *    `other`. In this mode, `add_source()` is forbidden — only `other`
 *    (the pool owner) may mutate the shared pool. This is what
 *    enables the `C(n,k)` shared-pool selection enforced by
 *    `ComboEnumerator`.
 *
 * Sources and transformations stored in an `F2LinearSourcePool` are
 * owned via `unique_ptr`; copies into the pool are deep copies made
 * through `Recurrence::copy()` / `Transformation::copy()`, so the
 * caller retains ownership of the originals.
 *
 * @code{.cpp}
 *   using regpoly::core::F2LinearSourcePool;
 *   F2LinearSourcePool p;
 *   p.add_source(some_lfsr_a);   // deep copy into pool
 *   p.add_source(some_lfsr_b);
 *   p.add_trans(temper_step);    // tempering chain on top of active source
 *
 *   F2LinearSourcePool p_share;
 *   p_share.copy_pool_from(p);   // now references p's pool
 *   // p_share.add_source(...);  // would throw — shared pool
 * @endcode
 *
 * @ingroup core
 */
class F2LinearSourcePool {
public:
    /// Owned pool of sources (each `unique_ptr` is the sole owner).
    /// Storage widened to `F2LinearSource` in Phase 5.2-storage so the
    /// pool can hold any single F_2-linear bit source — `Recurrence`
    /// PRNGs or `DigitalNet` instances. `build_combined_from_enumerator`
    /// rejects pools that mix Recurrence + non-Recurrence sources for
    /// now (the search loop's runners are still typed on `Recurrence`);
    /// the wider type lets non-search callers manipulate digital-net
    /// pools too.
    using SourcePool = std::vector<std::unique_ptr<F2LinearSource>>;
    /// Ordered tempering chain applied on top of the active source.
    using TransChain = std::vector<std::unique_ptr<Transformation>>;

    /** @brief Construct an empty pool owning a fresh, empty source pool. */
    F2LinearSourcePool();
    F2LinearSourcePool(const F2LinearSourcePool&) = delete;
    F2LinearSourcePool& operator=(const F2LinearSourcePool&) = delete;
    F2LinearSourcePool(F2LinearSourcePool&&) = default;
    F2LinearSourcePool& operator=(F2LinearSourcePool&&) = default;

    /**
     * @brief Append a deep copy of `src` to this pool.
     *
     * `src` may be any `F2LinearSource` subclass — a `Recurrence` PRNG
     * or a `DigitalNet`. The deep copy goes through
     * `F2LinearSource::clone_source()`.
     *
     * @param src  Source to deep-copy.
     * @throws std::runtime_error  If this pool currently shares its
     *                             underlying storage with another (only
     *                             the pool owner may mutate it).
     */
    void add_source(const F2LinearSource& src);

    /**
     * @brief Append a deep copy of `t` to this pool's tempering chain.
     *
     * @param t  Transformation to deep-copy via `Transformation::copy()`.
     */
    void add_trans(const Transformation& t);

    /**
     * @brief Adopt `other`'s pool (currently via shared_ptr aliasing).
     *
     * Phase 2 step 9 (Q1): renamed from the legacy `share_pool_with` to
     * `copy_pool_from` to reflect the *intended* end-state semantics —
     * each slot's pool is a deep copy of the donor's, so per-slot
     * tempering chains (now living on each `F2LinearSource::tempering_`
     * intrinsically after step 6a) stay independent across slots. The
     * actual deep-copy flip is deferred until the search-enumeration
     * logic (`compute_min_index`, `already_used`, both keyed on
     * `pool_id()`) migrates from "same-pool" to "same-pool-of-origin"
     * identity. Until that flip, this still aliases the pool via
     * `shared_ptr` — same observable behaviour as the old
     * `share_pool_with` minus the misleading name.
     *
     * @param other  Pool whose underlying storage to adopt.
     */
    void copy_pool_from(const F2LinearSourcePool& other);

    /** @brief Number of sources in the (possibly shared) pool. */
    int nb_sources() const;
    /** @brief Number of transformations in the tempering chain. */
    int nb_trans() const;
    /** @brief Current pool index (`-1` until `ComboEnumerator` places it). */
    int current_index() const { return current_index_; }
    /**
     * @brief Set the current pool index without validation.
     * @param i  Pool index to make active.
     */
    void set_current_index(int i) { current_index_ = i; }

    /**
     * @brief Access the source at pool index `i` without advancing state.
     * @param i  Pool index in `[0, nb_sources())`.
     * @return   Reference to the pooled source (lifetime == this pool).
     */
    F2LinearSource& source_at(int i) const;

    /**
     * @brief Access the currently active source (`source_at(current_index())`).
     * @return Reference to the active source.
     */
    F2LinearSource& active_source() const;

    /** @brief Read-only view of the tempering chain. */
    const TransChain& trans() const { return trans_; }
    /**
     * @brief Access the transformation at chain index `i`.
     * @param i  Index in `[0, nb_trans())`.
     * @return   Reference to the chain element.
     */
    Transformation& trans_at(int i) const;

    /**
     * @brief Pool identity for shared-pool detection.
     *
     * Two pools share storage iff their `pool_id()` pointers compare
     * equal. `ComboEnumerator` uses this to apply the shared-pool
     * C(n,k) constraint.
     *
     * @return  Raw pointer to the underlying `SourcePool` (do not
     *          dereference for ownership; use only for identity
     *          comparison).
     */
    const SourcePool* pool_id() const { return pool_.get(); }

    /**
     * @brief Concatenate `display_str()` of every tempering step.
     *
     * @return  Newline-separated tempering chain description.
     */
    std::string display() const;

private:
    std::shared_ptr<SourcePool> pool_;     // never null after construction
    TransChain trans_;
    int current_index_;
    bool owns_pool_;
};


/**
 * @brief Stateful iterator over the cartesian product of J component pools.
 *
 * The search loop drives a `ComboEnumerator` to enumerate every legal
 * tuple of `(active generator at slot 0, ..., active generator at
 * slot J-1)` subject to the two semantic constraints captured at the
 * top of this header:
 *
 *  1. Identity-uniqueness (no two slots reference the same
 *     `Recurrence` object).
 *  2. Shared-pool C(n,k) selection (slots that share a pool see
 *     strictly increasing indices).
 *
 * Iteration is driven through three operations:
 *
 *  - `reset()` places the iterator at the first legal combo and
 *    returns false if no combo exists (e.g. an empty component pool).
 *  - `next()` advances to the next legal combo, returning false on
 *    exhaustion. Subsequent calls keep returning false (idempotent).
 *  - `at(j)` returns the active generator at slot `j` once placed.
 *
 * After every successful placement, `ComboEnumerator` recomputes `k_g`
 * (sum of active `Recurrence::k()`) and `L` (the minimum active
 * `Recurrence::L()`, capped at the per-combination ceiling `Lmax`).
 *
 * @code{.cpp}
 *   using regpoly::core::ComboEnumerator;
 *   ComboEnumerator comb(2, 32);   // J = 2 slots, Lmax = 32
 *   // ... configure comb.pool(0) / comb.pool(1) with sources ...
 *   if (comb.reset()) {
 *       do {
 *           regpoly::core::F2LinearSource& g0 = comb.at(0);
 *           regpoly::core::F2LinearSource& g1 = comb.at(1);
 *           // ... evaluate the combo, e.g. via build_combined_from_enumerator ...
 *       } while (comb.next());
 *   }
 * @endcode
 *
 * @see :py:class:`regpoly.core.combination.ComboEnumerator`
 *
 * @ingroup core
 */
class ComboEnumerator {
public:
    /**
     * @brief Construct an unplaced combination of `J` slots.
     *
     * After construction, configure each slot via `pool(j)`
     * (adding sources / transformations or wiring up shared
     * pools) before calling `reset()`.
     *
     * @param J     Number of component slots.
     * @param Lmax  Ceiling for the recomputed `L` of any combo.
     */
    ComboEnumerator(int J, int Lmax);
    ComboEnumerator(const ComboEnumerator&) = delete;
    ComboEnumerator& operator=(const ComboEnumerator&) = delete;
    ComboEnumerator(ComboEnumerator&&) = default;
    ComboEnumerator& operator=(ComboEnumerator&&) = default;

    /** @brief Number of component slots. */
    int J() const { return J_; }
    /** @brief Ceiling on the recomputed `L`. */
    int Lmax() const { return Lmax_; }
    /** @brief Sum of active `Recurrence::k()` across slots (after placement). */
    int k_g() const { return k_g_; }
    /** @brief Minimum active `Recurrence::L()`, capped at `Lmax` (after placement). */
    int L() const { return L_; }

    /**
     * @brief Mutable access to slot `j`'s source pool.
     * @param j  Slot index in `[0, J())`.
     * @return   Reference to the slot's `F2LinearSourcePool`.
     */
    F2LinearSourcePool& pool(int j);
    /**
     * @brief Read-only access to slot `j`'s source pool.
     * @param j  Slot index in `[0, J())`.
     * @return   Const reference to the slot's `F2LinearSourcePool`.
     */
    const F2LinearSourcePool& pool(int j) const;

    /**
     * @brief Active source at slot `j` (equivalent to Python's `comb[j]`).
     *
     * @pre   `reset()` has returned true (or a prior `next()` did).
     * @param j  Slot index in `[0, J())`.
     * @return   Reference to the active source at slot `j`.
     */
    F2LinearSource& at(int j) const;

    /**
     * @brief Place the iterator at the first legal combination.
     *
     * @return  True if a legal combo exists; false otherwise (e.g.
     *          some component has an empty pool, or the shared-pool
     *          constraint admits no placement).
     */
    bool reset();

    /**
     * @brief Advance to the next legal combination.
     *
     * After exhaustion, subsequent calls keep returning false (the
     * iterator is sticky on exhaustion).
     *
     * @return  True on success; false on exhaustion.
     */
    bool next();

    /** @brief True iff the iterator has been exhausted. */
    bool exhausted() const { return exhausted_; }

private:
    int J_;
    int Lmax_;
    int k_g_;
    int L_;
    std::vector<std::shared_ptr<F2LinearSourcePool>> pools_;
    std::vector<int> indices_;   // indices_[j] = current index in
                                 // components_[j]'s pool; -1 = unset
    bool exhausted_;

    // Recompute k_g and L from the current active generators.
    void update_stats();

    // Place valid indices at slots j..J-1 starting from a fresh state.
    // On success, indices_[j..J-1] are populated and components_'
    // current_gen are set.
    bool place_from(int j);

    // Try to advance the index at slot j (and reset slots j+1..J-1).
    // On exhaustion at slot j, recursively try slot j-1.
    bool advance_from(int j);

    // Compute the minimum allowed start index at slot j based on
    // shared-pool constraints with slots 0..j-1.
    int compute_min_index(int j) const;

    // Check whether a source pointer at slot j collides with any
    // already-placed slot < j.
    bool already_used(const F2LinearSource* g, int j) const;
};


/**
 * @brief Build an independent combined source from `comb`'s current state.
 *
 * Each active component is deep-cloned via `clone_source()` /
 * `clone_recurrence()` and each per-slot tempering chain is cloned
 * via `Transformation::copy()`. The returned source is fully
 * independent of `comb` — subsequent `next()` / `reset()` calls on
 * `comb` do not affect the returned object, and conversely driving
 * the returned source does not perturb the iterator.
 *
 * Dispatch by pool content (Phase 5.5):
 *
 * - If every active source is a `Recurrence`, returns an
 *   `ITestable`-typed `CombinedF2LinearSource` — preserves the legacy
 *   inheritance shape and the per-slot tempering-chain plumbing
 *   (`F2LinearSourcePool::trans_` installed onto each component).
 * - Otherwise (at least one `DigitalNet` in some slot), returns a
 *   `CombinedF2LinearSource` — the composition wrapper that accepts
 *   heterogeneous F_2-linear sources. Per-slot tempering chains are
 *   installed onto each component's intrinsic `F2LinearSource::tempering_`
 *   so they still take effect; no top-level chain.
 *
 * Callers should consume the result as `ITestable&` since the
 * concrete shape depends on the pool composition. Lattice-family
 * predicates downstream throw if any component is a non-Recurrence
 * (the kernels can't operate on digital-net state evolution); the
 * matricial / t-value / tuplets / collision-free predicates handle
 * either shape polymorphically.
 *
 * @pre   `comb.reset()` (or a prior `next()`) returned true.
 * @param comb  ComboEnumerator iterator positioned on a legal combo.
 * @return      A heap-allocated, independent combined source.
 */
std::unique_ptr<ITestable>
build_combined_from_enumerator(const ComboEnumerator& comb);

}  // namespace regpoly::core
