// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "transformation.h"
#include "param_spec.h"
#include <memory>
#include <string>

/**
 * @file permutation.h
 * @brief Bit-pair permutation tempering (`permut`).
 *
 * Registered with the transformation factory under the canonical
 * name `permut`. Swaps the two bits at positions `p` and `q` of the
 * output word — a discrete, non-bitmask tempering used by some
 * lattice-method variants where a single bit-pair exchange
 * minimises the dimension defect.
 *
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief Single bit-pair swap on the output word.
 *
 * Implements the linear permutation `y[p] ↔ y[q]` on the `w`-bit
 * output. Both `p` and `q` are structural — they are not bitmasks
 * and the tempering optimiser does not mutate them between trials.
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   Params p;
 *   p.set_int("w", 32);
 *   p.set_int("p", 5);
 *   p.set_int("q", 17);
 *   auto t = create_transformation("permut", p);
 * @endcode
 *
 * @see :py:class:`regpoly.core.transformation.Transformation`  — the Python wrapper.
 * @ingroup core
 */
class PermutationTrans : public Transformation {
public:
    /**
     * @brief Construct a PermutationTrans with explicit parameters.
     * @param w  Word width in bits.
     * @param p  First bit position (`0 <= p < w`).
     * @param q  Second bit position (`0 <= q < w`, distinct from `p`).
     */
    PermutationTrans(int w, int p, int q);

    /**
     * @brief Factory hook — build from a Params dict.
     * @param params  Parameter dict with keys `w, p, q`.
     * @return        Owning pointer to the new PermutationTrans.
     * @throws std::runtime_error  If a required parameter is missing.
     */
    static std::unique_ptr<Transformation> from_params(const Params& params);

    /** @brief Parameter specs declared by this transformation. */
    static std::vector<ParamSpec> param_specs();

    /** @brief Canonical type name — returns `"permut"`. */
    std::string name() const override;

    /** @brief Human-readable parameter dump. */
    std::string display_str() const override;

    /** @brief Apply the bit-pair swap in place to the output word. */
    void apply(BitVect& state) const override;

    /** @brief Deep clone — independent of `*this`. */
    std::unique_ptr<Transformation> copy() const override;

    /**
     * @brief No-op for PermutationTrans — both `p` and `q` are structural.
     *
     * Present to satisfy the `Transformation` virtual interface; the
     * tempering optimiser will never mutate this kind of transform.
     */
    void update(const Params& params) override;

    /** @brief First bit position swapped. */
    int p() const { return p_; }
    /** @brief Second bit position swapped. */
    int q() const { return q_; }

private:
    int p_, q_;
};

}  // namespace regpoly::core
