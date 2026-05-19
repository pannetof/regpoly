// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "generator.h"
#include "param_spec.h"
#include "params.h"
#include <cstdint>
#include <memory>
#include <string>

/**
 * @file tinymt32.h
 * @brief 127-bit-state TinyMT generator (Saito & Matsumoto 2011).
 *
 * `TinyMT32Gen` is the 127-bit-state Mersenne-Twister-style generator
 * with one tempering parameter, from Saito & Matsumoto 2011. The
 * reference implementation lives in
 * `MTToolBox/samples/TinyMTDC/tinymt32search.hpp`. The catalog file
 * `docs/library/saito-matsumoto-2011-tinymt.yaml` carries the
 * published TinyMT32 parameter set.
 *
 * State layout in `state_` (BitVect, MSB-first; status[0] high bit
 * is not part of the F_2-linear state per the upstream isZero()
 * mask 0x7fffffff, so the effective recurrence dimension is 127
 * even though the BitVect carries 128 bits):
 *
 * @verbatim
 * bits   0.. 31 = status[0]
 * bits  32.. 63 = status[1]
 * bits  64.. 95 = status[2]
 * bits  96..127 = status[3]
 * @endverbatim
 *
 * Recurrence (from `tinymt32search.hpp::next_state`):
 *
 * @verbatim
 * y = status[3]
 * x = (status[0] & 0x7fffffff) ^ status[1] ^ status[2]
 * x ^= (x << 1)
 * y ^= (y >> 1) ^ x
 * status[0] = status[1]
 * status[1] = status[2]
 * status[2] = x ^ (y << 10)
 * status[3] = y
 * if (y & 1) {
 *     status[1] ^= mat1
 *     status[2] ^= mat2
 * }
 * @endverbatim
 *
 * Tempering (from `tinymt32search.hpp::temper`):
 *
 * @verbatim
 * t1 = status[0] ^ (status[2] >> 8)
 * out = status[3] ^ t1
 * if (t1 & 1) out ^= tmat
 * return out
 * @endverbatim
 *
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief TinyMT32 generator (Saito & Matsumoto 2011).
 *
 * Structural parameters: `mat1`, `mat2` (32-bit recurrence matrix
 * constants) and `tmat` (32-bit tempering constant). Registered as
 * `"TinyMT32Gen"` (alias `"TinyMT32"`). The default catalog entry
 * uses `mat1 = 0x8f7011ee, mat2 = 0xfc78ff1f, tmat = 0x3793fdff`.
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   Params p;
 *   p.set_int("mat1", 0x8f7011ee);
 *   p.set_int("mat2", 0xfc78ff1f);
 *   p.set_int("tmat", 0x3793fdff);
 *   auto gen = create_generator("TinyMT32Gen", p, 32);
 * @endcode
 *
 * @see :py:class:`regpoly.core.generator.Generator`
 * @ingroup core
 */
class TinyMT32Gen : public Generator {
public:
    /**
     * @brief Construct a TinyMT32Gen with explicit matrix constants.
     *
     * Most callers should go through `create_generator("TinyMT32Gen", ...)`
     * rather than this constructor directly.
     *
     * @param mat1  First recurrence matrix constant (32 bits).
     * @param mat2  Second recurrence matrix constant (32 bits).
     * @param tmat  Tempering constant (32 bits).
     * @param L     Output resolution in bits.
     */
    TinyMT32Gen(uint32_t mat1, uint32_t mat2, uint32_t tmat, int L);

    /**
     * @brief Build a TinyMT32Gen from a Params dict (registry factory hook).
     * @param params  Parameter dict with keys mat1, mat2, tmat.
     * @param L       Output resolution in bits.
     * @return        A constructed TinyMT32Gen as a polymorphic Generator pointer.
     * @throws std::runtime_error  If a required parameter is missing or invalid.
     */
    static std::unique_ptr<Generator> from_params(const Params& params, int L);

    /**
     * @brief Parameter specs declared by this family.
     * @return  Vector of ParamSpec records consumed by the factory and the
     *          Python introspection helpers.
     */
    static std::vector<ParamSpec> param_specs();

    /** @brief Family display name — returns the canonical "TinyMT32Gen" string. */
    std::string name() const override;
    /** @brief Human-readable parameter summary (mat1, mat2, tmat). */
    std::string display_str() const override;
    /** @brief Seed the state from the leading 128 bits of `init_bv`. */
    void init(const BitVect& init_bv) override;
    /** @brief Apply one TinyMT32 step (recurrence + tempering). */
    void next() override;
    /** @brief Deep copy this generator (state included). */
    std::unique_ptr<Generator> copy() const override;
    /** @brief Most recently produced tempered output truncated to `L()` bits. */
    BitVect get_output() const override;

    /** @brief First recurrence matrix constant. */
    uint32_t mat1() const { return mat1_; }
    /** @brief Second recurrence matrix constant. */
    uint32_t mat2() const { return mat2_; }
    /** @brief Tempering constant. */
    uint32_t tmat() const { return tmat_; }

    /** @brief Period exponent — fixed at 127 for TinyMT32. */
    int period_exponent() const override { return 127; }

private:
    uint32_t mat1_, mat2_, tmat_;
    uint32_t last_output_;

    void load_status(uint32_t s[4]) const;
    void store_status(const uint32_t s[4]);
    uint32_t compute_temper(const uint32_t s[4]) const;
};

}  // namespace regpoly::core
