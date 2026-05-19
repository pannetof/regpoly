// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "generator.h"
#include "transformation.h"
#include <memory>
#include <utility>
#include <vector>

/**
 * @file combined.h
 * @brief Meta-generator that XORs J component generators together.
 *
 * Defines `CombinedGenerator`, the polymorphic `Generator` subclass
 * the rest of the search loop and the equidistribution / lattice
 * machinery operates on. A `CombinedGenerator` owns `J` component
 * `Generator`s plus a per-component tempering chain; `next()`
 * advances every component in lockstep and `get_output()` returns
 * the bitwise XOR of each component's (post-tempering) `L`-bit
 * output.
 *
 * Mathematically, the combined characteristic polynomial is the
 * product over GF(2) of the components' characteristic polynomials,
 * and the combined transition matrix is block-diagonal across
 * components. This is the API uniform across primitive and combined
 * generators that lets every kernel consume a single
 * `const Generator&` — whether the underlying generator is a single
 * MT19937 or a J-component combined construction.
 *
 * In Phase 1, tempering chains and the full SIMD-aware story (lane
 * counts, output phases) are stored but `simd_lane_count()` /
 * `output_phases()` fall through to the base defaults; the search
 * loop materialises a `CombinedGenerator` from a `Combination`
 * iterator via `regpoly::core::build_combined_from_combination`.
 *
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief Meta-generator: XORs J component generators (each with its own tempering chain).
 *
 * `CombinedGenerator` is the polymorphic outer object the search loop
 * and the equidistribution / lattice analyses operate on. It owns:
 *
 *  - `J >= 1` component `Generator`s (each themselves a `Generator` —
 *    primitive families like MT, WELL, SFMT, or even nested combined
 *    generators).
 *  - A per-component tempering chain (`vector<unique_ptr<Transformation>>`).
 *  - The associated `Lmax` ceiling for the published output width.
 *
 * `next()` advances every component in lockstep; `get_output()`
 * XORs each component's post-tempering output, truncating to `L()`.
 *
 * The `components()` and `tempering_chains()` overrides expose the
 * inner state to the kernel adapter layer (the base `Generator`
 * class returns the trivial single-component view; combined
 * generators override both to surface their J slots).
 *
 * Relationship to `regpoly::core::Combination`:
 *   - `Combination` is the iterator over the cartesian product of J
 *     candidate pools.
 *   - `build_combined_from_combination(comb)` materialises the
 *     iterator's current snapshot into a fully independent
 *     `CombinedGenerator` (deep-copies every active component and
 *     every tempering step).
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   Combination comb(2, 32);
 *   // ... configure component pools and tempering chains ...
 *   if (comb.reset()) {
 *       auto combined = build_combined_from_combination(comb);
 *       combined->init(seed_bv);
 *       combined->next();
 *       BitVect out = combined->get_output();   // XOR of both components
 *   }
 * @endcode
 *
 * @see :py:class:`regpoly.core.combination.Combination` — the Python
 *      counterpart wraps both `Combination` and `CombinedGenerator`
 *      under a single class.
 * @ingroup core
 */
class CombinedGenerator : public Generator {
public:
    /// Per-component tempering chain (ordered list of `Transformation`s).
    using ComponentTempering =
        std::vector<std::unique_ptr<Transformation>>;

    /**
     * @brief Construct a `CombinedGenerator` with no per-component tempering.
     *
     * Each component is taken into the wrapper with its own
     * tempering already applied (or no tempering at all). Use the
     * other constructor when callers want the wrapper to own the
     * tempering chains.
     *
     * @param components  Owned component generators (transferred in).
     * @param Lmax        Ceiling on the published output width.
     */
    CombinedGenerator(
        std::vector<std::unique_ptr<Generator>> components,
        int Lmax);

    /**
     * @brief Construct a `CombinedGenerator` with per-component tempering chains.
     *
     * The chains are applied to each component's output before the
     * XOR step in `get_output()`.
     *
     * @param components         Owned component generators (transferred in).
     * @param tempering_chains   Per-component owned tempering chains
     *                           (one chain per component, possibly empty).
     * @param Lmax               Ceiling on the published output width.
     */
    CombinedGenerator(
        std::vector<std::unique_ptr<Generator>> components,
        std::vector<ComponentTempering> tempering_chains,
        int Lmax);

    /** @brief Family display name — "CombinedGenerator" plus J. */
    std::string name() const override;
    /** @brief Human-readable summary including every component's `display_str()`. */
    std::string display_str() const override;

    /**
     * @brief Seed every component from the leading bits of `init_bv`.
     *
     * `init_bv` must carry at least `k()` bits. The leading `k_0`
     * bits go to component 0, the next `k_1` to component 1, and so
     * on. If `init_bv` is too short it is zero-padded; if too long
     * the trailing bits are ignored.
     */
    void init(const BitVect& init_bv) override;

    /** @brief Advance every component generator by one step (in lockstep). */
    void next() override;

    /**
     * @brief Return the XOR of components' tempered L-bit outputs, truncated to `L()`.
     */
    BitVect get_output() const override;

    /** @brief Deep copy the wrapper, every component, and every tempering chain. */
    std::unique_ptr<Generator> copy() const override;

    /**
     * @brief Combined characteristic polynomial — product over GF(2) of the components'.
     *
     * Returns `k()` bits (no leading `z^k` term, matching the
     * `Generator::char_poly()` convention).
     */
    BitVect char_poly() const override;

    /** @brief Number of component generators (`J`). */
    int J() const { return static_cast<int>(components_.size()); }
    /** @brief Read-only access to component `j`. */
    const Generator& component(int j) const { return *components_[j]; }

    /**
     * @brief Cumulative `k` partition across components.
     *
     * `prefix_k_[j] = sum_{i<j} k_i`, and `prefix_k_[J] == k()`.
     * Useful for slicing seeds and outputs to the right component.
     */
    const std::vector<int>& prefix_k() const { return prefix_k_; }

    /**
     * @brief Raw component pointers (lattice / equidistribution adapter helper).
     *
     * Internal accessor used by the kernel adapter overloads. Phase 1
     * keeps the kernel implementations on their original
     * `vector<Generator*>` + per-component tempering API; the new
     * `const Generator&` overloads dispatch through these accessors
     * to recover the vector form. Phase 2+ migrates the
     * implementations to consume `Generator&` natively, at which
     * point these can become private again.
     */
    std::vector<Generator*> raw_component_pointers() const;
    /** @brief Raw per-component tempering chain pointers (adapter helper). */
    std::vector<std::vector<Transformation*>>
        raw_tempering_pointers() const;

    /**
     * @brief Expose component decomposition for polymorphic dispatch.
     *
     * Overrides the base `Generator::components()` to surface the J
     * inner generators rather than treating `*this` as a single
     * primitive — polymorphic replacement for the `dynamic_cast`
     * pattern formerly used at every kernel adapter site.
     */
    std::vector<Generator*> components() const override {
        return raw_component_pointers();
    }
    /** @brief Per-component tempering chains for polymorphic dispatch. */
    std::vector<std::vector<Transformation*>> tempering_chains() const override {
        return raw_tempering_pointers();
    }

protected:
    /**
     * @brief Default test-method override: `notprimitive` for J >= 2.
     *
     * The combined chi is the product of components' chi over GF(2),
     * hence reducible whenever there are >= 2 components of positive
     * degree — `notprimitive` is the only safe choice. The J=1
     * degenerate case defers to the wrapped component to preserve
     * the JEquals1MatchesPrimitive invariant.
     */
    std::optional<std::string> compute_default_test_method(const std::string& test_type) const override;

private:
    std::vector<std::unique_ptr<Generator>> components_;
    std::vector<ComponentTempering> tempering_;
    std::vector<int> prefix_k_;

    int compute_k(
        const std::vector<std::unique_ptr<Generator>>& comps);
    int compute_L(
        const std::vector<std::unique_ptr<Generator>>& comps,
        int Lmax);

    void refresh_concatenated_state();
};

}  // namespace regpoly::core
