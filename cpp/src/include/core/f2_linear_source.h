// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.

#pragma once

#include "bitvect.h"
#include "i_testable.h"
#include "tempering_chain.h"
#include <cassert>
#include <memory>
#include <vector>

/**
 * @file f2_linear_source.h
 * @brief Abstract base for everything that IS a single F_2-linear bit
 *        source — a PRNG (`Recurrence`) or a digital net (`DigitalNet`).
 *
 * Implements the `ITestable` slot common to leaf sources: owns a
 * `TemperingChain` member, fulfils `k()` / `L()` / `sources()` once,
 * and exposes a tempered `get_output()` (raw_output composed with the
 * source's intrinsic chain). Subclasses implement `raw_output()` plus
 * the remaining `ITestable` virtuals (`name`, `display_str`, `init`,
 * `next`, `copy`, `default_test_method`).
 *
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief Abstract base for single F_2-linear bit sources.
 *
 * Owns the tempering chain by value (copyable). Subclasses implement
 * `raw_output()`; the base exposes both `get_output()` (raw during
 * the Phase 2 transition; will flip to tempered at the end) and
 * `tempered_output()` (always tempered, kernels migrate to this).
 *
 * @ingroup core
 */
class F2LinearSource : public ITestable {
public:
    F2LinearSource(int k, int L) : k_(k), L_(L) {}

    // Copyable + movable. The implicit defaults are well-formed now
    // that `TemperingChain` is copyable. Leaf-class `copy()` overrides
    // (`return std::make_unique<MTGen>(*this);`) rely on this chain.
    F2LinearSource(const F2LinearSource&) = default;
    F2LinearSource& operator=(const F2LinearSource&) = default;
    F2LinearSource(F2LinearSource&&) noexcept = default;
    F2LinearSource& operator=(F2LinearSource&&) noexcept = default;

    // ── ITestable fulfilment ─────────────────────────────────────────

    int k() const final { return k_; }
    int L() const final { return L_; }

    /**
     * @brief Return the tempered output: chain ∘ `raw_output()`.
     *
     * The architectural target of the Phase 2 redesign — every kernel
     * reads tempered bits via this single method. Subclasses implement
     * `raw_output()` (the un-tempered output); the chain composition is
     * handled here using the source's intrinsic `tempering_` chain.
     */
    BitVect get_output() const final {
        BitVect out = raw_output();
        tempering_.apply(out);
        return out;
    }

    // Phase 3: dropped `final` so `CombinedF2LinearSource` (which inherits
    // F2LinearSource via Recurrence) can override to expose its J
    // components instead of `{this}`. The pure leaf F2LinearSource
    // (PRNGs, DigitalNets) keeps the default.
    std::vector<const F2LinearSource*> sources() const override {
        return {this};
    }

    // ── New pure-virtual: subclasses implement the un-tempered read ──

    /**
     * @brief Return the current un-tempered output word.
     *
     * Concrete subclasses produce the family-defined L-bit slice of
     * their state (e.g. top L bits of the MT state, or
     * `C_j * input_index` for a digital net coordinate). The base
     * class wraps this with the tempering chain via `tempered_output()`.
     *
     * Hot path: called once per Gaussian-elimination row × per basis
     * vector × per coordinate. Avoid unnecessary allocations.
     */
    virtual BitVect raw_output() const = 0;

    // ── Typed-clone forwarder (non-virtual) ──────────────────────────

    /**
     * @brief Deep clone, typed as `unique_ptr<F2LinearSource>`.
     *
     * Wraps the polymorphic `copy()` virtual with a `static_cast` that
     * recovers the `F2LinearSource` static type. Any `F2LinearSource`
     * subclass's `copy()` must return an `F2LinearSource`-derived
     * pointer (enforced by the architecture; debug builds assert via
     * `dynamic_cast`). The `Recurrence::clone_recurrence()` tier
     * introduced in Phase 2 Step 3 follows the same pattern.
     *
     * Used by `F2LinearSourcePool::add_source` and
     * `CombinedF2LinearSource` construction, which need
     * `F2LinearSource`-typed storage without dynamic-casting at every
     * callsite.
     */
    std::unique_ptr<F2LinearSource> clone_source() const {
        auto p = copy();
        assert(dynamic_cast<F2LinearSource*>(p.get()) != nullptr
               && "F2LinearSource subclass copy() must return an "
                  "F2LinearSource-derived pointer");
        return std::unique_ptr<F2LinearSource>(
            static_cast<F2LinearSource*>(p.release()));
    }

    // ── TemperingChain ownership ─────────────────────────────────────

    /** @brief Read-only view of the source's tempering chain. */
    const TemperingChain& tempering() const { return tempering_; }
    /** @brief Mutable view (for in-place randomisation by the optimiser). */
    TemperingChain& tempering() { return tempering_; }
    /** @brief Replace the chain wholesale (move-in). */
    void set_tempering(TemperingChain chain) { tempering_ = std::move(chain); }

protected:
    int k_;
    int L_;
    TemperingChain tempering_;
};

}  // namespace regpoly::core
