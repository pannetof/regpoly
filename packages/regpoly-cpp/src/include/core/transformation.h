// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "bitvect.h"
#include "params.h"
#include <string>
#include <memory>

/**
 * @file transformation.h
 * @brief Abstract base for every output-tempering transformation.
 *
 * Concrete subclasses implement the F_2-linear output-tempering steps
 * (shift / mask / XOR families) that REGPOLY chains on top of a
 * `Generator` output. Each step transforms an output `BitVect` in
 * place via `apply()`; the tempering optimiser mutates the step's
 * bitmask parameters between trials via `update()`; the search loop
 * deep-clones via `copy()` so candidate transformations can be
 * evaluated independently.
 *
 * @see :py:class:`regpoly.core.transformation.Transformation`
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief Abstract base for every output-tempering transformation.
 *
 * Subclasses implement a single F_2-linear transformation applied
 * on top of a `Generator`'s output word. They are owned by
 * `Component` (via `unique_ptr`) and assembled into per-component
 * tempering chains by `Combination`. The optimiser explores the
 * parameter space of a given chain by calling `update(params)`
 * between trials.
 *
 * @see :py:class:`regpoly.core.transformation.Transformation`
 *
 * @ingroup core
 */
class Transformation {
public:
    /** @brief Construct an unparameterised transformation (width set by `update()`). */
    Transformation() : w_(0) {}
    virtual ~Transformation() = default;

    /** @brief Canonical step name (e.g. `"tempMK"`, `"shiftL"`). */
    virtual std::string name() const = 0;

    /**
     * @brief Human-readable parametrised string for diagnostics.
     *
     * Typically includes the step name and the current mask /
     * shift parameter values.
     *
     * @return  Display string.
     */
    virtual std::string display_str() const = 0;

    /**
     * @brief Apply the F_2-linear transformation to `state` in place.
     *
     * @param state  Output word to transform (mutated in place).
     */
    virtual void apply(BitVect& state) const = 0;

    /**
     * @brief Return an independent deep clone of this transformation.
     *
     * Used by `Component::add_trans()` and
     * `build_combined_from_combination()` to keep candidate chains
     * decoupled from their source.
     *
     * @return  A `unique_ptr` to a freshly heap-allocated clone.
     */
    virtual std::unique_ptr<Transformation> copy() const = 0;

    /**
     * @brief Re-parameterise the step from a `Params` bundle.
     *
     * Called by the tempering optimiser between trials to mutate the
     * step's bitmask / shift parameters in place. Implementations
     * must update `w_` to the new effective width when applicable.
     *
     * @param params  New parameter values keyed by name.
     */
    virtual void update(const Params& params) = 0;

    /** @brief Effective bit-width of the step (set by `update()`). */
    int w() const { return w_; }

protected:
    int w_;   ///< Effective bit-width (subclass-maintained).
};

}  // namespace regpoly::core
