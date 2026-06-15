// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.

#pragma once

#include "bitvect.h"
#include "transformation.h"
#include <memory>
#include <vector>

/**
 * @file tempering_chain.h
 * @brief First-class wrapper around an ordered list of `Transformation`s.
 *
 * Replaces the bare `std::vector<std::unique_ptr<Transformation>>` that
 * `Component` and `CombinedF2LinearSource` used to carry as raw members.
 * The chain owns its steps; `apply(BitVect&)` composes every step in
 * order against the supplied output word (mutate-in-place — no
 * heap allocations on the hot path).
 *
 * Phase 1 introduces the type and leaves it unused; Phase 2 migrates
 * tempering ownership from `Component` and the kernel-side
 * `vector<Transformation*>` arg-lists into this class.
 *
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief Ordered chain of `Transformation` steps applied on top of an
 *        `F2LinearSource`'s raw output.
 *
 * Copyable (deep-clones every step) and movable. Copy semantics are
 * load-bearing — `F2LinearSource` holds a `TemperingChain` by value,
 * and leaf-class `copy()` overrides rely on the implicit copy-ctor
 * chain (`std::make_unique<MTGen>(*this)`) walking through this type.
 *
 * @ingroup core
 */
class TemperingChain {
public:
    /** @brief Empty chain (identity transformation). */
    TemperingChain() = default;

    /**
     * @brief Construct from an existing vector of steps (takes ownership).
     *
     * @param steps  Vector of owned `Transformation`s; moved in.
     */
    explicit TemperingChain(
        std::vector<std::unique_ptr<Transformation>> steps)
        : steps_(std::move(steps)) {}

    // Copyable (deep-clones every step) + movable (defaulted).
    TemperingChain(const TemperingChain& other) {
        steps_.reserve(other.steps_.size());
        for (const auto& s : other.steps_) {
            steps_.push_back(s->copy());
        }
    }
    TemperingChain& operator=(const TemperingChain& other) {
        if (this != &other) {
            steps_.clear();
            steps_.reserve(other.steps_.size());
            for (const auto& s : other.steps_) {
                steps_.push_back(s->copy());
            }
        }
        return *this;
    }
    TemperingChain(TemperingChain&&) noexcept = default;
    TemperingChain& operator=(TemperingChain&&) noexcept = default;

    /** @brief Append a step (takes ownership). */
    void add(std::unique_ptr<Transformation> step) {
        steps_.push_back(std::move(step));
    }

    /** @brief Number of steps in the chain. */
    int size() const { return static_cast<int>(steps_.size()); }
    /** @brief True iff the chain has no steps (identity). */
    bool empty() const { return steps_.empty(); }

    /**
     * @brief Access the i-th step.
     * @param i  Index in `[0, size())`.
     * @throws std::out_of_range  If `i` is outside the range.
     */
    Transformation& at(int i) const {
        return *steps_.at(static_cast<std::size_t>(i));
    }

    /** @brief Read-only view of the underlying step vector. */
    const std::vector<std::unique_ptr<Transformation>>& steps() const {
        return steps_;
    }

    /**
     * @brief Apply every step in order against `state` (mutate-in-place).
     *
     * No heap allocations: each step's `Transformation::apply(BitVect&)`
     * mutates `state` in place, matching the hot-path call site in the
     * Gaussian-elimination row builder.
     *
     * @param state  Output word; mutated in place.
     */
    void apply(BitVect& state) const {
        for (const auto& step : steps_) {
            step->apply(state);
        }
    }

private:
    std::vector<std::unique_ptr<Transformation>> steps_;
};

}  // namespace regpoly::core
