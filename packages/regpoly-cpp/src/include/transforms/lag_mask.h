// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "transformation.h"
#include "param_spec.h"
#include <cstdint>

/**
 * @file lag_mask.h
 * @brief Two-reference lagged-mask output tempering (MELG-style `lag` chain).
 *
 * The tempering canonical name is `lag` in the transformation registry.
 * Designed for the MELG generator family (Harase–Kimoto 2018): it
 * combines a one-word left-shift fold with an AND-mask draw from a
 * second word at a fixed lag in the state array. Works with any
 * generator whose F₂-linear state can be viewed as an array of
 * `w`-bit words.
 *
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief Two-reference output tempering used by MELG.
 *
 * Applies, in order, to the output word `y`:
 *
 *     y  = y ^ (y << sigma)              // T1: left-shift XOR
 *     y ^= state_word_at_lag_L & b       // T2: lag-mask AND
 *
 * where `state_word_at_lag_L` is the `L`-th `w`-bit word of the
 * generator's state array (bit offset `L * w`). The structural fields
 * `w`, `sigma`, `L` are fixed at construction; the bitmask `b` is
 * the tempering-optimisable knob.
 *
 * Construction goes through `create_transformation("lag", params)`;
 * direct constructor calls are for testing and for the catalog loader.
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   Params p;
 *   p.set_int("w",     64);
 *   p.set_int("sigma", 23);
 *   p.set_int("L",     11);     // lag in 64-bit words
 *   p.set_int("b",     0x5555555555555555ULL);
 *   auto t = create_transformation("lag", p);
 * @endcode
 *
 * @see :py:class:`regpoly.core.transformation.Transformation`  — the Python wrapper.
 * @ingroup core
 */
class LaggedTempering : public Transformation {
public:
    /**
     * @brief Construct a LaggedTempering with explicit parameters.
     *
     * @param w      Word width in bits.
     * @param sigma  Left-shift width for T1 (`0 < sigma < w`).
     * @param L      Word-index lag for T2 (`0 < L < num_words`).
     * @param b      AND-mask for T2 (`w`-bit, optimizable).
     */
    LaggedTempering(int w, int sigma, int L, uint64_t b);

    /**
     * @brief Factory hook — build from a Params dict.
     * @param params  Parameter dict with keys `w, sigma, L, b`.
     * @return        Owning pointer to the new LaggedTempering.
     * @throws std::runtime_error  If a required parameter is missing.
     */
    static std::unique_ptr<Transformation> from_params(const Params& params);

    /** @brief Parameter specs declared by this transformation. */
    static std::vector<ParamSpec> param_specs();

    /** @brief Canonical type name — returns `"lag"`. */
    std::string name() const override;

    /** @brief Human-readable parameter dump. */
    std::string display_str() const override;

    /** @brief Apply the two-step tempering in place to the output word. */
    void apply(BitVect& state) const override;

    /** @brief Deep clone — independent of `*this`. */
    std::unique_ptr<Transformation> copy() const override;

    /**
     * @brief Mutate `b` from `params` (tempering-optimiser hook).
     *
     * Structural fields (`w`, `sigma`, `L`) are not touched.
     */
    void update(const Params& params) override;

private:
    int sigma_;
    int L_;
    uint64_t b_;
};

}  // namespace regpoly::core
