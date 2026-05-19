// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "generator.h"
#include "param_spec.h"
#include <memory>
#include <string>

/**
 * @file tgfsr.h
 * @brief Twisted Generalized Feedback Shift Register (Matsumoto & Kurita 1992).
 *
 * `TGFSRGen` implements the original twisted GFSR recurrence: a
 * word-level F_2-linear recurrence whose state is `r` words of width
 * `w`, advanced by XORing `state[i]` with a twisted copy of
 * `state[(i+m) % r]`. The twist matrix is encoded as a bit-vector
 * `a` of length `w` (its bottom row); upper rows are the identity
 * shift. TGFSR is the direct predecessor of MT — MT replaces TGFSR's
 * full-word twist with a partial-bit (split) recurrence.
 *
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief Twisted GFSR generator (Matsumoto & Kurita 1992).
 *
 * Structural parameters: `w` (word width), `r` (state size in
 * words), `m` (twist offset), and `a` (twist matrix bottom row,
 * carried as a `BitVect` of length `w`). Output is the active state
 * word truncated to `L` bits. Registered as `"TGFSRGen"` (alias
 * `"TGFSR"`).
 *
 * @code{.cpp}
 *   // Construct via the catalog or the YAML config loader — see
 *   // regpoly::library::Catalog::generator() or
 *   // yaml_config::load_seek_config().  Direct factory construction:
 *   //   auto gen = create_generator("TGFSRGen", params, L);
 *   // where `params` carries w, r, m, a — see the family's
 *   // param_specs() for the per-field constraints.
 * @endcode
 *
 * @see :py:class:`regpoly.core.generator.Generator`
 * @ingroup core
 */
class TGFSRGen : public Generator {
public:
    /**
     * @brief Construct a TGFSRGen with explicit twist parameters.
     *
     * Most callers should go through `create_generator("TGFSRGen", ...)`
     * rather than this constructor directly.
     *
     * @param w  Word width in bits.
     * @param r  State size in words.
     * @param m  Twist offset.
     * @param a  Twist matrix bottom row (BitVect of length `w`).
     * @param L  Output resolution in bits.
     */
    TGFSRGen(int w, int r, int m, const BitVect& a, int L);

    /**
     * @brief Build a TGFSRGen from a Params dict (registry factory hook).
     * @param params  Parameter dict with keys w, r, m, a.
     * @param L       Output resolution in bits.
     * @return        A constructed TGFSRGen as a polymorphic Generator pointer.
     * @throws std::runtime_error  If a required parameter is missing or invalid.
     */
    static std::unique_ptr<Generator> from_params(const Params& params, int L);

    /**
     * @brief Parameter specs declared by this family.
     * @return  Vector of ParamSpec records consumed by the factory and the
     *          Python introspection helpers.
     */
    static std::vector<ParamSpec> param_specs();

    /** @brief Family display name — returns the canonical "TGFSRGen" string. */
    std::string name() const override;
    /** @brief Human-readable parameter summary including twist constant. */
    std::string display_str() const override;
    /** @brief Seed the state from `init_bv` (`w * r` bits consumed). */
    void init(const BitVect& init_bv) override;
    /** @brief Apply one TGFSR recurrence step. */
    void next() override;
    /** @brief Deep copy this generator (state included). */
    std::unique_ptr<Generator> copy() const override;
    /** @brief Characteristic polynomial of the F_2-linear recurrence. */
    BitVect char_poly() const override;

    /** @brief Word width in bits. */
    int w() const { return w_; }
    /** @brief State size in words. */
    int r() const { return r_; }
    /** @brief Twist offset. */
    int m() const { return m_; }
    /** @brief Twist matrix bottom row as a BitVect. */
    const BitVect& a() const { return a_; }

private:
    int w_, r_, m_;
    BitVect a_;
};

}  // namespace regpoly::core
