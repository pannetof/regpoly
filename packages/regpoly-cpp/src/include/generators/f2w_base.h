// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "generator.h"
#include "gf2w_arith.h"
#include <vector>
#include <cstdint>

/**
 * @file f2w_base.h
 * @brief Abstract base for generators in F_{2^w} (Panneton & L'Ecuyer 2004).
 *
 * Shared structural parameters and arithmetic helpers for generators
 * that combine word-level GF(2^w) arithmetic with a base polynomial
 * recurrence. Concrete subclasses are `F2wLFSRGen` (LFSR form) and
 * `F2wPolyLCGGen` (polynomial-LCG form); paper Tables 1 & 2 list
 * each parameter set in both forms because they share the same
 * characteristic polynomial under the equidistribution analysis.
 *
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief Discriminator tag for the concrete F2w generator type (LFSR vs PolyLCG).
 */
enum GenF2wType { GENF2W_POLYLCG = 0, GENF2W_LFSR = 1 };

/**
 * @brief Abstract base for F_2[t] / GF(2^w) generators (Panneton & L'Ecuyer 2004).
 *
 * Carries the structural parameters common to LFSR and PolyLCG
 * forms: word width `w`, state size `r` (in F_{2^w} elements), the
 * connection polynomial encoded as `(nbcoeff, nocoeff, coeff)`, the
 * modular reduction polynomial `modM`, and the basis flag
 * `normal_basis`. Multiplication and state access in GF(2^w) are
 * shared across the two subclasses via the protected helpers
 * `multiply()`, `V()`, and `SetV()`.
 *
 * Not user-constructible directly — see `F2wLFSRGen` and
 * `F2wPolyLCGGen` for the concrete factories.
 *
 * @see :py:class:`regpoly.core.generator.Generator`
 * @ingroup core
 */
class F2wBaseGen : public Generator {
protected:
    int w_, r_;
    int nbcoeff_;
    std::vector<int> nocoeff_;
    std::vector<uint64_t> coeff_;
    uint64_t modM_;
    bool normal_basis_;
    std::vector<uint64_t> table_;
    uint64_t maskw_;
    int step_count_;

public:
    /**
     * @brief Construct the shared F_{2^w} state and arithmetic tables.
     *
     * Called from each concrete subclass's constructor. Most callers
     * should go through `create_generator(...)` on the subclass name.
     *
     * @param w             Word width in bits (size of each F_{2^w} element).
     * @param r             Number of F_{2^w} state elements.
     * @param nbcoeff       Number of nonzero interior coefficients in `nocoeff`/`coeff`.
     * @param nocoeff       Positions (degrees) of the interior nonzero coefficients.
     * @param coeff         Values of those coefficients in F_{2^w}.
     * @param modM          Reduction polynomial M(z) of degree w (bit-packed, LSB-first).
     * @param normal_basis  True iff GF(2^w) is represented in normal basis.
     * @param step_count    Number of base-recurrence steps to apply per `next()`.
     * @param L             Output resolution in bits.
     */
    F2wBaseGen(int w, int r, int nbcoeff,
               const std::vector<int>& nocoeff,
               const std::vector<uint64_t>& coeff,
               uint64_t modM, bool normal_basis,
               int step_count, int L);

    /** @brief Family base-name; concrete subclasses override `type_name()` for more detail. */
    std::string name() const override { return "Generator in F_{2^w}"; }
    /** @brief Pure virtual — subclasses implement state seeding from `init_bv`. */
    void init(const BitVect& init_bv) override = 0;
    /** @brief Pure virtual — subclasses implement the per-step advance. */
    void next() override = 0;
    /** @brief Pure virtual — subclasses implement deep copy. */
    std::unique_ptr<Generator> copy() const override = 0;

    /** @brief Human-readable parameter summary (calls `type_name()` for the form prefix). */
    std::string display_str() const override;
    /** @brief Subclass-supplied "LFSR in F_{2^w}" / "Polynomial LCG in F_{2^w}[z]/P(z)" prefix. */
    virtual std::string type_name() const = 0;

    /** @brief Word width `w` of GF(2^w). */
    int gf2w_w() const { return w_; }
    /** @brief State size `r` in F_{2^w} elements. */
    int gf2w_r() const { return r_; }
    /** @brief Number of nonzero interior coefficients in the connection polynomial. */
    int nbcoeff() const { return nbcoeff_; }
    /** @brief Positions (degrees) of the nonzero interior coefficients. */
    const std::vector<int>& nocoeff() const { return nocoeff_; }
    /** @brief Values of the nonzero interior coefficients in F_{2^w}. */
    const std::vector<uint64_t>& coeff() const { return coeff_; }
    /** @brief Modular reduction polynomial M(z) of degree w (bit-packed). */
    uint64_t modM() const { return modM_; }
    /** @brief True iff GF(2^w) is represented in normal basis. */
    bool normal_basis() const { return normal_basis_; }
    /** @brief Number of base-recurrence steps applied per `next()`. */
    int step_count() const { return step_count_; }
    /** @brief Precomputed multiplication table (basis-dependent). */
    const std::vector<uint64_t>& table() const { return table_; }

protected:
    /** @brief GF(2^w) multiplication of two w-bit values modulo `modM_`. */
    uint64_t multiply(uint64_t a, uint64_t b) const;
    /** @brief Read the F_{2^w} state element at index `idx`. */
    uint64_t V(int idx) const;
    /** @brief Write `val` into the F_{2^w} state element at index `idx`. */
    void SetV(int idx, uint64_t val);
};

}  // namespace regpoly::core
