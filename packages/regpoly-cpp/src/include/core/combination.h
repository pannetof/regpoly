// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "combined.h"
#include "generator.h"
#include "transformation.h"

#include <cstdint>
#include <memory>
#include <vector>

/**
 * @file combination.h
 * @brief Cartesian-product iterator over component generator pools.
 * @defgroup core Core
 *
 * Defines the core types used by the search loop to enumerate
 * candidate combined generators: `Component` (one slot's owned or
 * shared pool of generators plus a tempering chain) and `Combination`
 * (a stateful iterator over the cartesian product of J such slots).
 *
 * The iterator enforces two semantic constraints that mirror the
 * original Python implementation in `regpoly.core.{component,combination}`:
 *
 *  1. **Identity-uniqueness.** The same `Generator` object never
 *     appears twice in a single combo — XORing a component with
 *     itself produces the zero sequence, which is never useful.
 *  2. **Shared-pool C(n,k) selection.** When two component slots
 *     reference the same underlying pool object (set up via
 *     `Component::share_pool_with`), the indices in those slots
 *     must be strictly increasing. A pool of size `n` shared across
 *     `k` slots therefore yields `C(n,k)` combos instead of `n^k`
 *     permutations.
 *
 * On each successful `reset()` / `next()`, `Combination` recomputes
 * `k_g` (sum of active generator `k`'s) and `L` (min active `L`,
 * capped at `Lmax`) so callers can read the active combo's effective
 * state size and output width without re-walking the components.
 *
 * @see :py:class:`regpoly.core.combination.Combination`
 * @see :py:class:`regpoly.core.component.Component`
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief One slot of a `Combination`: a pool of generators plus a tempering chain.
 *
 * A `Component` owns (or shares) a pool of candidate `Generator`
 * instances and an ordered list of `Transformation`s that the search
 * loop applies on top of the active generator's output. Two ownership
 * modes are supported:
 *
 *  - **Own pool.** Default after construction. Use `add_gen()` to
 *    append deep copies of candidate generators.
 *  - **Shared pool.** After `share_pool_with(other)`, this component
 *    references the same underlying `GenPool` (via `shared_ptr`) as
 *    `other`. In this mode, `add_gen()` is forbidden — only `other`
 *    (the pool owner) may mutate the shared pool. This is what
 *    enables the `C(n,k)` shared-pool selection enforced by
 *    `Combination`.
 *
 * Generators and transformations stored in a `Component` are owned
 * via `unique_ptr`; copies into the component are deep copies made
 * through `Generator::copy()` / `Transformation::copy()`, so the
 * caller retains ownership of the originals.
 *
 * @code{.cpp}
 *   using regpoly::core::Component;
 *   Component c;
 *   c.add_gen(some_lfsr_a);   // deep copy into pool
 *   c.add_gen(some_lfsr_b);
 *   c.add_trans(temper_step); // tempering chain on top of active gen
 *
 *   Component c_share;
 *   c_share.share_pool_with(c);  // now references c's pool
 *   // c_share.add_gen(...);     // would throw — shared pool
 * @endcode
 *
 * @see :py:class:`regpoly.core.component.Component`
 *
 * @ingroup core
 */
class Component {
public:
    /// Owned pool of generators (each `unique_ptr` is the sole owner).
    using GenPool = std::vector<std::unique_ptr<Generator>>;
    /// Ordered tempering chain applied on top of the active generator.
    using TransChain = std::vector<std::unique_ptr<Transformation>>;

    /** @brief Construct an empty component owning a fresh, empty pool. */
    Component();
    Component(const Component&) = delete;
    Component& operator=(const Component&) = delete;
    Component(Component&&) = default;
    Component& operator=(Component&&) = default;

    /**
     * @brief Append a deep copy of `gen` to this component's own pool.
     *
     * @param gen  Generator to deep-copy via `Generator::copy()`.
     * @throws std::runtime_error  If this component currently shares
     *                             its pool with another (only the
     *                             pool owner may mutate it).
     */
    void add_gen(const Generator& gen);

    /**
     * @brief Append a deep copy of `t` to this component's tempering chain.
     *
     * @param t  Transformation to deep-copy via `Transformation::copy()`.
     */
    void add_trans(const Transformation& t);

    /**
     * @brief Make this component reference `other`'s pool.
     *
     * After this call, the two components share the same underlying
     * `GenPool` via `shared_ptr` semantics. `add_gen()` on this
     * component is forbidden from then on — only the original owner
     * (`other`) may grow the pool.
     *
     * @param other  Component to share the pool of.
     */
    void share_pool_with(const Component& other);

    /** @brief Number of generators in the (possibly shared) pool. */
    int nb_gen() const;
    /** @brief Number of transformations in the tempering chain. */
    int nb_trans() const;
    /** @brief Current pool index (`-1` until `Combination` places it). */
    int current_gen() const { return current_gen_; }
    /**
     * @brief Set the current pool index without validation.
     * @param i  Pool index to make active.
     */
    void set_current_gen(int i) { current_gen_ = i; }

    /**
     * @brief Access the generator at pool index `i` without advancing state.
     * @param i  Pool index in `[0, nb_gen())`.
     * @return   Reference to the pooled generator (lifetime == this component).
     */
    Generator& gen_at(int i) const;

    /**
     * @brief Access the currently active generator (`gen_at(current_gen())`).
     * @return Reference to the active generator.
     */
    Generator& active_gen() const;

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
     * Two components share a pool iff their `pool_id()` pointers
     * compare equal. `Combination` uses this to apply the
     * shared-pool C(n,k) constraint.
     *
     * @return  Raw pointer to the underlying `GenPool` (do not
     *          dereference for ownership; use only for identity
     *          comparison).
     */
    const GenPool* pool_id() const { return pool_.get(); }

    /**
     * @brief Concatenate `display_str()` of every tempering step.
     *
     * Mirrors Python's `Component.display()`. Useful for debug
     * logging and inclusion in result tables.
     *
     * @return  Newline-separated tempering chain description.
     */
    std::string display() const;

private:
    std::shared_ptr<GenPool> pool_;     // never null after construction
    TransChain trans_;
    int current_gen_;
    bool owns_pool_;
};


/**
 * @brief Stateful iterator over the cartesian product of J component pools.
 *
 * The search loop drives a `Combination` to enumerate every legal
 * tuple of `(active generator at slot 0, ..., active generator at
 * slot J-1)` subject to the two semantic constraints captured at the
 * top of this header:
 *
 *  1. Identity-uniqueness (no two slots reference the same
 *     `Generator` object).
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
 * After every successful placement, `Combination` recomputes `k_g`
 * (sum of active `Generator::k()`) and `L` (the minimum active
 * `Generator::L()`, capped at the per-combination ceiling `Lmax`).
 *
 * @code{.cpp}
 *   using regpoly::core::Combination;
 *   Combination comb(2, 32);   // J = 2 slots, Lmax = 32
 *   // ... configure comb.component(0) / comb.component(1) with pools ...
 *   if (comb.reset()) {
 *       do {
 *           regpoly::core::Generator& g0 = comb.at(0);
 *           regpoly::core::Generator& g1 = comb.at(1);
 *           // ... evaluate the combo, e.g. via build_combined_from_combination ...
 *       } while (comb.next());
 *   }
 * @endcode
 *
 * @see :py:class:`regpoly.core.combination.Combination`
 *
 * @ingroup core
 */
class Combination {
public:
    /**
     * @brief Construct an unplaced combination of `J` slots.
     *
     * After construction, configure each slot via `component(j)`
     * (adding generators / transformations or wiring up shared
     * pools) before calling `reset()`.
     *
     * @param J     Number of component slots.
     * @param Lmax  Ceiling for the recomputed `L` of any combo.
     */
    Combination(int J, int Lmax);
    Combination(const Combination&) = delete;
    Combination& operator=(const Combination&) = delete;
    Combination(Combination&&) = default;
    Combination& operator=(Combination&&) = default;

    /** @brief Number of component slots. */
    int J() const { return J_; }
    /** @brief Ceiling on the recomputed `L`. */
    int Lmax() const { return Lmax_; }
    /** @brief Sum of active `Generator::k()` across slots (after placement). */
    int k_g() const { return k_g_; }
    /** @brief Minimum active `Generator::L()`, capped at `Lmax` (after placement). */
    int L() const { return L_; }

    /**
     * @brief Mutable access to component slot `j`.
     * @param j  Slot index in `[0, J())`.
     * @return   Reference to the slot's `Component`.
     */
    Component& component(int j);
    /**
     * @brief Read-only access to component slot `j`.
     * @param j  Slot index in `[0, J())`.
     * @return   Const reference to the slot's `Component`.
     */
    const Component& component(int j) const;

    /**
     * @brief Active generator at slot `j` (equivalent to Python's `comb[j]`).
     *
     * @pre   `reset()` has returned true (or a prior `next()` did).
     * @param j  Slot index in `[0, J())`.
     * @return   Reference to the active generator at slot `j`.
     */
    Generator& at(int j) const;

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
    std::vector<std::shared_ptr<Component>> components_;
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

    // Check whether a generator pointer at slot j collides with any
    // already-placed slot < j.
    bool already_used(const Generator* g, int j) const;
};


/**
 * @brief Build a `CombinedGenerator` that snapshots the current state of `comb`.
 *
 * Each active component generator is deep-copied via
 * `Generator::copy()` and each component's tempering chain is cloned
 * via `Transformation::copy()`. The returned `CombinedGenerator` is
 * fully independent of `comb` — subsequent `next()` / `reset()` calls
 * on `comb` do not affect the returned object, and conversely
 * driving the returned generator does not perturb the iterator.
 *
 * @pre   `comb.reset()` (or a prior `next()`) returned true.
 * @param comb  Combination iterator positioned on a legal combo.
 * @return      A heap-allocated, independent `CombinedGenerator`.
 */
std::unique_ptr<CombinedGenerator>
build_combined_from_combination(const Combination& comb);

}  // namespace regpoly::core
