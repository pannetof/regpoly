// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "generator.h"
#include "param_spec.h"
#include <memory>
#include <string>

/**
 * @file polylcg.h
 * @brief Polynomial LCG over GF(2) — the canonical "z·state mod P(z)" recurrence.
 *
 * `PolyLCGGen` advances a `k`-bit state by polynomial multiplication
 * by `z` modulo a fixed polynomial `P(z)` of degree `k` over GF(2).
 * The recurrence is the simplest F_2-linear primitive in the
 * codebase: when `P(z)` is primitive of degree `k`, the generator
 * has period `2^k - 1`. Registered as `"PolyLCGGen"` (alias
 * `"PolyLCG"`).
 *
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief Polynomial LCG over GF(2): state(t+1) = z · state(t) mod P(z).
 *
 * The single structural argument is the connection polynomial
 * `poly`, packed in a `BitVect` whose bit i is the coefficient of
 * `z^i`. The state size `k` is the polynomial's degree. Output is
 * the top `L` bits of the state.
 *
 * @code{.cpp}
 *   // Construct via the catalog or the YAML config loader — see
 *   // regpoly::library::Catalog::generator() or
 *   // yaml_config::load_seek_config().  Direct factory construction:
 *   //   auto gen = create_generator("PolyLCGGen", params, L);
 *   // where `params` carries `k` and `poly` (the latter is a
 *   // bit-packed primitive polynomial of degree k).
 * @endcode
 *
 * @see :py:class:`regpoly.core.generator.Generator`
 * @ingroup core
 */
class PolyLCGGen : public Generator {
public:
    /**
     * @brief Construct a PolyLCGGen with explicit connection polynomial.
     *
     * Most callers should go through `create_generator("PolyLCGGen", ...)`
     * rather than this constructor directly.
     *
     * @param k     State size in bits (= degree of `poly`).
     * @param poly  Connection polynomial, packed as a BitVect.
     * @param L     Output resolution in bits.
     */
    PolyLCGGen(int k, const BitVect& poly, int L);

    /**
     * @brief Build a PolyLCGGen from a Params dict (registry factory hook).
     * @param params  Parameter dict with keys k and poly.
     * @param L       Output resolution in bits.
     * @return        A constructed PolyLCGGen as a polymorphic Generator pointer.
     * @throws std::runtime_error  If a required parameter is missing or invalid.
     */
    static std::unique_ptr<Generator> from_params(const Params& params, int L);

    /**
     * @brief Parameter specs declared by this family.
     * @return  Vector of ParamSpec records consumed by the factory and the
     *          Python introspection helpers.
     */
    static std::vector<ParamSpec> param_specs();

    /** @brief Family display name — returns the canonical "PolyLCGGen" string. */
    std::string name() const override;
    /** @brief Human-readable parameter summary (includes the polynomial). */
    std::string display_str() const override;
    /** @brief Seed the state from the leading bits of `init_bv`. */
    void init(const BitVect& init_bv) override;
    /** @brief Multiply the state by `z` modulo `poly()` once. */
    void next() override;
    /** @brief Deep copy this generator (state included). */
    std::unique_ptr<Generator> copy() const override;
    /** @brief Returns `poly()` (the recurrence is z modulo poly). */
    BitVect char_poly() const override;

    /** @brief The connection polynomial used by the recurrence. */
    const BitVect& poly() const { return poly_; }

private:
    BitVect poly_;
};

}  // namespace regpoly::core
