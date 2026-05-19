// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "generator.h"
#include "param_spec.h"
#include "params.h"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

/**
 * @file mtgp.h
 * @brief Mersenne Twister for Graphic Processors (Saito & Matsumoto 2013).
 *
 * `MTGPGen` is the Mersenne Twister for Graphic Processors (Saito &
 * Matsumoto 2013). Parameters (per published configuration —
 * MTGP11213, MTGP23209, MTGP44497):
 *
 *  - `mexp` (structural)   Mersenne exponent; state length `N = mexp / 32 + 1`.
 *  - `pos`  (structural)   pick-up offset in the recursion.
 *  - `sh1`, `sh2`          shift amounts in the recursion.
 *  - `mask`                upper-bit mask applied to `status[i]` — holds
 *                          the `32 - r` high bits where `r = 32*N - mexp`.
 *  - `tbl[16]`             recursion lookup table (low 4 bits of `y`).
 *  - `tmp_tbl[16]`         tempering lookup table (low 4 bits of `t`).
 *
 * Recurrence (reference: mtgp-dc `mtgp32_recursion`):
 *
 * @verbatim
 * y = (x1 & mask) ^ x2 ^ (y_in << sh1)
 * MAT = tbl[y & 0x0f]
 * status_new = y ^ (y >> sh2) ^ MAT
 * @endverbatim
 *
 * where `x1 = state[i], x2 = state[(i+1) % N], y_in = state[(i+pos) % N]`.
 *
 * Tempering (reference: `mtgp32_temper`):
 *
 * @verbatim
 * t' = t ^ (t >> 16)
 * output = status_new ^ tmp_tbl[t' & 0x0f]
 * @endverbatim
 *
 * where `t = state[(i + pos - 1) % N]` (the word just before the pick-up).
 *
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief Mersenne Twister for Graphic Processors (Saito & Matsumoto 2013).
 *
 * Registered as `"MTGPGen"` (alias `"MTGP"`). Constructed against
 * pre-computed `tbl` / `tmp_tbl` tables, so direct factory
 * construction is uncommon — most callers go through the catalog.
 *
 * @code{.cpp}
 *   // Construct via the catalog or the YAML config loader — see
 *   // regpoly::library::Catalog::generator() or
 *   // yaml_config::load_seek_config().  Direct factory construction
 *   // requires the 16-entry tbl / tmp_tbl lookup tables produced by
 *   // MTGPDC; see docs/library/saito-matsumoto-2013-mtgp.yaml for
 *   // ready-to-use parameter sets.
 *   //   auto gen = create_generator("MTGPGen", params, L);
 * @endcode
 *
 * @see :py:class:`regpoly.core.generator.Generator`
 * @ingroup core
 */
class MTGPGen : public Generator {
public:
    /**
     * @brief Construct an MTGPGen with explicit recursion + tempering tables.
     *
     * Most callers should go through `create_generator("MTGPGen", ...)`
     * rather than this constructor directly.
     *
     * @param mexp     Mersenne exponent.
     * @param pos      Pick-up offset.
     * @param sh1      Left-shift amount in the recursion.
     * @param sh2      Right-shift amount in the recursion.
     * @param mask     Upper-bit mask applied to status[i].
     * @param tbl      Recursion lookup table (16 entries).
     * @param tmp_tbl  Tempering lookup table (16 entries).
     * @param L        Output resolution in bits.
     */
    MTGPGen(int mexp, int pos, int sh1, int sh2, uint32_t mask,
         const std::vector<uint32_t>& tbl,
         const std::vector<uint32_t>& tmp_tbl,
         int L);

    /**
     * @brief Build an MTGPGen from a Params dict (registry factory hook).
     * @param params  Parameter dict with keys mexp, pos, sh1, sh2, mask,
     *                tbl, tmp_tbl.
     * @param L       Output resolution in bits.
     * @return        A constructed MTGPGen as a polymorphic Generator pointer.
     * @throws std::runtime_error  If a required parameter is missing or invalid.
     */
    static std::unique_ptr<Generator> from_params(
        const Params& params, int L);

    /**
     * @brief Parameter specs declared by this family.
     * @return  Vector of ParamSpec records consumed by the factory and the
     *          Python introspection helpers.
     */
    static std::vector<ParamSpec> param_specs();

    /** @brief Family display name — returns the canonical "MTGPGen" string. */
    std::string name() const override;
    /** @brief Human-readable parameter summary. */
    std::string display_str() const override;
    /** @brief Seed the state from the leading `32 * N` bits of `init_bv`. */
    void init(const BitVect& init_bv) override;
    /** @brief Apply one MTGP step (recurrence + tempering). */
    void next() override;
    /** @brief Deep copy this generator (state included). */
    std::unique_ptr<Generator> copy() const override;
    /** @brief Most recently produced tempered output truncated to `L()` bits. */
    BitVect get_output() const override;

    /** @brief Mersenne exponent. */
    int mexp() const { return mexp_; }
    /** @brief State length in 32-bit words. */
    int N() const { return N_; }

protected:
    /** @brief Default test method — `simd_notprimitive` (char poly is reducible). */
    std::optional<std::string> compute_default_test_method(const std::string& test_type) const override;

private:
    int mexp_;
    int N_;                    // = mexp/32 + 1
    int pos_, sh1_, sh2_;
    uint32_t mask_;
    std::vector<uint32_t> tbl_;      // 16 entries
    std::vector<uint32_t> tmp_tbl_;  // 16 entries
    int idx_;                  // ring head 0 ≤ idx_ < N_
    uint32_t last_output_;     // cached result of the most recent step

    uint32_t get_word(int i) const {
        return (uint32_t)state_.get_word(i, 32);
    }
    void set_word(int i, uint32_t v) {
        state_.set_word(i, 32, (uint64_t)v);
    }
};

}  // namespace regpoly::core
