// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.

#pragma once

#include "bitvect.h"
#include "f2_linear_source.h"
#include "generator.h"      // Recurrence (for recurrence_components())
#include "i_testable.h"
#include "tempering_chain.h"
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

/**
 * @file combined_f2_linear_source.h
 * @brief Composition wrapper: XOR of J component `F2LinearSource`s.
 *
 * `CombinedF2LinearSource` is a **composition** of J F_2-linear bit
 * sources. It implements `ITestable` directly — it is NOT a
 * `Recurrence` itself, by design. Components can be any
 * `F2LinearSource` (Recurrence-driven PRNGs *and* DigitalNet point
 * sources). Kernels that genuinely need recurrence-typed components
 * (matricial χ-recovery, lattice spectral test, SIMD path) walk
 * `components()` and `dynamic_cast` each one — they throw
 * `std::invalid_argument` if any component is not a Recurrence.
 *
 * Each component carries its own intrinsic tempering chain on
 * `F2LinearSource::tempering_`. `next()` advances every component
 * in lockstep; `get_output()` is the XOR of each component's
 * already-tempered output, truncated to `L()`.
 *
 * @ingroup core
 */

namespace regpoly::core {

class CombinedF2LinearSource : public ITestable {
public:
    /**
     * @brief Construct from owned components.
     *
     * @param components  Owned F2LinearSource components (transferred
     *                    in). Each may be a Recurrence or DigitalNet.
     *                    Must be non-empty.
     * @param Lmax        Upper bound on the published output width.
     */
    CombinedF2LinearSource(
        std::vector<std::unique_ptr<F2LinearSource>> components,
        int Lmax);

    /**
     * @brief Construct with owned components plus per-component tempering chains.
     *
     * Each chain is moved onto its component's intrinsic
     * `F2LinearSource::tempering_`.
     *
     * @param components         Owned F2LinearSource components (transferred in).
     * @param tempering_chains   Per-component `TemperingChain`s.
     * @param Lmax               Upper bound on the published output width.
     */
    CombinedF2LinearSource(
        std::vector<std::unique_ptr<F2LinearSource>> components,
        std::vector<TemperingChain> tempering_chains,
        int Lmax);

    // ── ITestable fulfilment ─────────────────────────────────────────

    int k() const final { return k_; }
    int L() const final { return L_; }
    std::string name() const final;
    std::string display_str() const final;

    void init(const BitVect& init_bv) final;
    void next() final;
    BitVect get_output() const final;

    std::unique_ptr<ITestable> copy() const final;

    /**
     * @brief Per-component F2LinearSource views (each pointer aliases an
     *        owned component).
     */
    std::vector<const F2LinearSource*> sources() const final;

    std::optional<std::string>
        default_test_method(const std::string& test_type) const final;

    // ── Composition-specific accessors (kernels walk these) ──────────

    /** @brief Number of component sources (`J`). */
    int J() const { return static_cast<int>(components_.size()); }

    /** @brief Read-only access to component `j`. */
    const F2LinearSource& component(int j) const { return *components_[j]; }
    /** @brief Mutable access to component `j` (tempering search). */
    F2LinearSource& component(int j) { return *components_[j]; }

    /**
     * @brief Per-component F2LinearSource pointers — kernels walk
     *        these. Kernels that require Recurrence components
     *        should call `recurrence_components()` instead.
     */
    std::vector<F2LinearSource*> components() const;

    /**
     * @brief Per-component Recurrence pointers — for kernels that
     *        genuinely need recurrence-state evolution (matricial
     *        χ-recovery, lattice spectral test, SIMD path).
     *
     * `dynamic_cast`s each component to `Recurrence*` and throws
     * `std::invalid_argument` if any component is not a Recurrence
     * (e.g. a DigitalNet — which has no F_2-linear state to evolve).
     *
     * @param caller  Caller name embedded in the exception message
     *                (so the user can see which kernel rejected the
     *                composition). Defaults to a generic message.
     */
    std::vector<Recurrence*>
    recurrence_components(const char* caller = "kernel") const;

    /**
     * @brief Cumulative `k` partition across components.
     *
     * `prefix_k_[j] = sum_{i<j} k_i`, and `prefix_k_[J] == k()`.
     */
    const std::vector<int>& prefix_k() const { return prefix_k_; }

private:
    std::vector<std::unique_ptr<F2LinearSource>> components_;
    std::vector<int> prefix_k_;
    int k_;
    int L_;
};

}  // namespace regpoly::core
