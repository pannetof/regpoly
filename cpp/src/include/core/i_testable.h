// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.

#pragma once

#include "bitvect.h"
#include <memory>
#include <optional>
#include <string>
#include <vector>

/**
 * @file i_testable.h
 * @brief Pure-abstract interface that every type consumed by an analysis
 *        kernel implements.
 *
 * `ITestable` is the kernel-facing minimum surface — once the
 * hierarchy refactor lands, every C++ analysis kernel (matricial
 * equidistribution, dual-lattice, harase, notprimitive, t-value,
 * tuplets, collision-free) takes `const ITestable&` rather than a
 * concrete `Recurrence&`. Two distinct kinds of objects satisfy it:
 *
 * - `F2LinearSource` (a single F_2-linear bit source: any
 *   `Recurrence` PRNG or any `DigitalNet`). Implements `sources()`
 *   to return `{this}`.
 * - `CombinedF2LinearSource` (composition: an XOR of J
 *   `F2LinearSource`s, each carrying its own tempering chain).
 *   Implements `sources()` to return its J components.
 *
 * Phase 1 introduces the type and leaves it unused. Phase 3 will flip
 * every kernel signature to consume `const ITestable&`.
 *
 * @ingroup core
 */

namespace regpoly::core {

class F2LinearSource;  // forward declaration

/**
 * @brief Pure-abstract interface for everything analysis kernels
 *        accept as input.
 *
 * No state, no members; just the minimum polymorphic surface. Lifetime
 * of the returned `F2LinearSource*` from `sources()` is bound to the
 * `ITestable` instance — never store them beyond the lifetime of the
 * containing testable.
 *
 * @see regpoly::core::F2LinearSource
 * @see regpoly::core::CombinedF2LinearSource
 *
 * @ingroup core
 */
class ITestable {
public:
    virtual ~ITestable() = default;

    /** @brief Total F_2-linear state width in bits (sum of components'). */
    virtual int k() const = 0;
    /** @brief Output word width in bits. */
    virtual int L() const = 0;
    /** @brief Canonical name (family for a leaf, "CombinedF2LinearSource"
     *         for a composition). */
    virtual std::string name() const = 0;
    /** @brief Human-readable parametrised string for diagnostics. */
    virtual std::string display_str() const = 0;

    /** @brief Initialise from the given seed. Composition implementations
     *         slice the seed across components by per-component k. */
    virtual void init(const BitVect& seed) = 0;
    /** @brief Advance by one fundamental step. */
    virtual void next() = 0;
    /** @brief Return the current (tempered) output. For
     *         `CombinedF2LinearSource`, this is the XOR of each
     *         component's already-tempered output. */
    virtual BitVect get_output() const = 0;

    /** @brief Deep clone. Returns an independent `ITestable` that
     *         produces identical outputs given identical seeds.
     *
     *  Note: C++ does not allow covariant `std::unique_ptr` return
     *  types, so the typed-clone variants live on the subclass layers
     *  (`F2LinearSource::clone_source()`, `Recurrence::clone_recurrence()`)
     *  introduced in Phase 2.
     */
    virtual std::unique_ptr<ITestable> copy() const = 0;

    /**
     * @brief Decompose into a list of `F2LinearSource` views.
     *
     * - For a leaf `F2LinearSource`: returns `{this}`.
     * - For `CombinedF2LinearSource`: returns the J components.
     *
     * The returned pointers are non-owning and remain valid only for
     * the lifetime of `*this`. Used by the kernel `unpack_for_kernel`
     * to walk per-component bits.
     */
    virtual std::vector<const F2LinearSource*> sources() const = 0;

    /**
     * @brief Recommend a default test method for the given test type.
     *
     * Returns the canonical method name (e.g. `"matricial"`,
     * `"harase"`, `"notprimitive"`, `"schmid"`) the family prefers,
     * or `std::nullopt` if the test type is not supported. Phase 1
     * leaves the existing `Recurrence::default_test_method` plumbing
     * in place; Phase 3 lifts the signature here.
     */
    virtual std::optional<std::string>
        default_test_method(const std::string& test_type) const = 0;
};

}  // namespace regpoly::core
