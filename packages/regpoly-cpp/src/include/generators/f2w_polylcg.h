// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "f2w_base.h"
#include "param_spec.h"
#include <memory>
#include <string>

/**
 * @file f2w_polylcg.h
 * @brief Polynomial LCG over F_{2^w} (Panneton & L'Ecuyer 2004, PolyLCG form).
 *
 * Concrete subclass of `F2wBaseGen` implementing the polynomial-LCG
 * form of the GF(2^w) recurrence. Equivalent to `F2wLFSRGen` under
 * the characteristic-polynomial equivalence stated in Panneton &
 * L'Ecuyer 2004 §1 — see `docs/library/panneton-lecuyer-2004-f2w.yaml`
 * for the published parameter sets (each row appears under both
 * forms).
 *
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief Polynomial LCG over F_{2^w}[z] / P(z) (Panneton & L'Ecuyer 2004).
 *
 * The state is interpreted as a polynomial of degree `r-1` over
 * GF(2^w); each step multiplies by `z` modulo the connection
 * polynomial P(z). Registered as `"F2wPolyLCGGen"` (alias
 * `"GenF2wPolyLCG"`).
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   Params p;
 *   p.set_int("w", 32);
 *   p.set_int("r", 3);
 *   p.set_int("modM", 0xccb06f34);
 *   p.set_bool("normal_basis", false);
 *   p.set_int("step", 1);
 *   p.set_int_vec("nocoeff", {1, 0});
 *   p.set_int_vec("coeff", {0x30a72fa7, 0x537a531f});
 *   auto gen = create_generator("F2wPolyLCGGen", p, 32);
 * @endcode
 *
 * @see :py:class:`regpoly.core.generator.Generator`
 * @ingroup core
 */
class F2wPolyLCGGen : public F2wBaseGen {
public:
    /**
     * @brief Construct an F2wPolyLCGGen with explicit structural parameters.
     *
     * Arguments mirror `F2wLFSRGen`'s constructor exactly — the two
     * forms share the same parameter set and the difference is in
     * the `next()` semantics.
     *
     * @param w             Word width in bits (size of each F_{2^w} element).
     * @param r             Number of F_{2^w} state elements.
     * @param nbcoeff       Number of nonzero interior coefficients.
     * @param nocoeff       Positions of the interior nonzero coefficients.
     * @param coeff         Values of those coefficients in F_{2^w}.
     * @param modM          Reduction polynomial M(z) bit-packed.
     * @param normal_basis  True iff GF(2^w) uses a normal basis.
     * @param step_count    Base-recurrence steps per `next()`.
     * @param L             Output resolution in bits.
     */
    F2wPolyLCGGen(int w, int r, int nbcoeff,
                  const std::vector<int>& nocoeff,
                  const std::vector<uint64_t>& coeff,
                  uint64_t modM, bool normal_basis,
                  int step_count, int L);

    /**
     * @brief Build an F2wPolyLCGGen from a Params dict (registry factory hook).
     * @param params  Parameter dict with the keys w, r, nocoeff, coeff, modM,
     *                normal_basis, step.
     * @param L       Output resolution in bits.
     * @return        A constructed F2wPolyLCGGen as a polymorphic Generator pointer.
     * @throws std::runtime_error  If a required parameter is missing or invalid.
     */
    static std::unique_ptr<Generator> from_params(const Params& params, int L);

    /**
     * @brief Parameter specs declared by this family.
     * @return  Vector of ParamSpec records consumed by the factory and the
     *          Python introspection helpers.
     */
    static std::vector<ParamSpec> param_specs();

    /** @brief Family form-prefix used inside `display_str()`. */
    std::string type_name() const override { return "Polynomial LCG in F_{2^w}[z]/P(z)"; }
    /** @brief Seed the state from `init_bv`; the first `w * r` bits are used. */
    void init(const BitVect& init_bv) override;
    /** @brief Advance the polynomial-LCG recurrence by `step_count` base steps. */
    void next() override;
    /** @brief Deep copy this generator (state included). */
    std::unique_ptr<Generator> copy() const override;
};

}  // namespace regpoly::core
