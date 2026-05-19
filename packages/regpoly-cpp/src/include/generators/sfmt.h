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
 * @file sfmt.h
 * @brief SIMD-oriented Fast Mersenne Twister (Saito & Matsumoto 2008).
 *
 * `SFMTGen` is the SIMD-oriented Fast Mersenne Twister of Saito &
 * Matsumoto 2008. The state is an array of `N` 128-bit super-words,
 * where `N = MEXP / 128 + 1`. Each step applies the "do_recursion"
 * update to one array slot, advancing the ring index; one full
 * sweep of `N` steps corresponds to one `generate_all` block in
 * the reference implementation. Catalog file:
 * `docs/library/saito-matsumoto-2008.yaml`.
 *
 * This class stores the state as a `BitVect` of `128 * N` bits and
 * exposes each 128-bit word as four 32-bit lanes `u[0..3]`
 * (matching the reference implementation's little-endian within-word
 * layout).
 *
 * Parameters are the ten tabulated SFMTGen constants: the shift /
 * rotation amounts `POS1, SL1, SL2, SR1, SR2` and the 4-lane bitmask
 * `MSK1..MSK4`. The parity constants from the reference code are not
 * used here — the regpoly pipeline initialises the state to a
 * non-zero bit pattern via the search loop.
 *
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief SIMD-oriented Fast Mersenne Twister (Saito & Matsumoto 2008).
 *
 * Registered as `"SFMTGen"` (alias `"SFMT"`). The flagship
 * configuration `SFMT19937` uses `mexp = 19937, pos1 = 122, sl1 = 18,
 * sl2 = 1, sr1 = 11, sr2 = 1` with the catalog-supplied four-lane
 * masks.
 *
 * The characteristic polynomial of the `128 * N`-bit recurrence
 * factors as a degree-`MEXP` primitive factor times a low-degree
 * non-primitive remainder; equidistribution analysis must therefore
 * use the non-primitive method.
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   Params p;
 *   p.set_int("mexp", 19937);
 *   p.set_int("pos1", 122);
 *   p.set_int("sl1", 18);
 *   p.set_int("sl2", 1);
 *   p.set_int("sr1", 11);
 *   p.set_int("sr2", 1);
 *   p.set_int("msk1", 0xdfffffef);
 *   p.set_int("msk2", 0xddfecb7f);
 *   p.set_int("msk3", 0xbffaffff);
 *   p.set_int("msk4", 0xbffffff6);
 *   auto gen = create_generator("SFMTGen", p, 32);
 * @endcode
 *
 * @see :py:class:`regpoly.core.generator.Generator`
 * @ingroup core
 */
class SFMTGen : public Generator {
public:
    /**
     * @brief Construct an SFMTGen with explicit structural parameters.
     *
     * Most callers should go through `create_generator("SFMTGen", ...)`
     * rather than this constructor directly.
     *
     * @param mexp  Mersenne exponent (period exponent).
     * @param pos1  Pick-up offset in the recurrence.
     * @param sl1   First left shift amount.
     * @param sl2   Second left shift amount.
     * @param sr1   First right shift amount.
     * @param sr2   Second right shift amount.
     * @param msk1  First 32-bit lane mask.
     * @param msk2  Second 32-bit lane mask.
     * @param msk3  Third 32-bit lane mask.
     * @param msk4  Fourth 32-bit lane mask.
     * @param L     Output resolution in bits.
     */
    SFMTGen(int mexp, int pos1, int sl1, int sl2, int sr1, int sr2,
         uint32_t msk1, uint32_t msk2, uint32_t msk3, uint32_t msk4,
         int L);

    /**
     * @brief Build an SFMTGen from a Params dict (registry factory hook).
     * @param params  Parameter dict with keys mexp, pos1, sl1, sl2, sr1, sr2,
     *                msk1..msk4.
     * @param L       Output resolution in bits.
     * @return        A constructed SFMTGen as a polymorphic Generator pointer.
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

    /** @brief Family display name — returns the canonical "SFMTGen" string. */
    std::string name() const override;
    /** @brief Human-readable parameter summary. */
    std::string display_str() const override;
    /** @brief Seed the state from the leading `128 * N` bits of `init_bv`. */
    void init(const BitVect& init_bv) override;
    /** @brief Apply one SFMT recurrence step (one slot per call). */
    void next() override;
    /** @brief Deep copy this generator (state included). */
    std::unique_ptr<Generator> copy() const override;
    /** @brief Current `L`-bit output word from the active lane. */
    BitVect get_output() const override;

    /** @brief Mersenne exponent. */
    int mexp() const { return mexp_; }
    /** @brief State length in 128-bit super-words. */
    int N() const { return N_; }
    /** @brief Period exponent — equal to `mexp`. */
    int period_exponent() const override { return mexp_; }

    /**
     * @brief Number of output phases per state-cycle (= 128 / L).
     *
     * Per Saito-Matsumoto 2008, §3.2 Proposition 2: SFMTGen's `k(v)`
     * is computed on the augmented automaton `S' = S × {0, 1, …, P-1}`
     * where `P = 128 / L` is the number of L-bit lanes per 128-bit
     * word (`P = 4` for `L = 32`). Each phase `i`'s sequence reads
     * only lane `i` of every output 128-bit word, and
     * `k(v) = min_i k_i(v)`. The notprimitive slow path consults
     * this via the abstract `Generator` API and partitions the
     * output stream into `P` phases.
     */
    int output_phases() const override { return 128 / L_; }

    /**
     * @brief SIMD lane count for the SIMD-aware notprimitive method.
     *
     * SFMTGen packs `128 / L` lanes per 128-bit super-word (4 for
     * `L = 32`, 2 for `L = 64`).
     */
    int simd_lane_count() const override { return 128 / L_; }

    /**
     * @brief Install `s` as the raw state without triggering `generate_all`.
     *
     * Unlike `init()`, this does not advance the recurrence. Mirrors
     * the new state into `work_` so subsequent `next()` calls read
     * lanes of the freshly-set state. Required by the SIMD-PIS
     * reduction which needs to XOR generator states between
     * iterations without advancing.
     */
    void set_raw_state(const BitVect& s) override;

    /**
     * @brief Advance the SFMTGen by one 128-bit super-word (SIMD-aware override).
     *
     * Maintains the per-word advance index `pword_idx_` so successive
     * calls walk the ring without consuming lane outputs.
     */
    void simd_advance_one_word() override;
    /**
     * @brief Read the current 128-bit super-word for SIMD-PIS reduction.
     * @param start_mode  Lane-phase selector matching MTToolBox's convention.
     */
    BitVect simd_read_super_word(int start_mode) const override;
    /** @brief XOR `other`'s state into this one's (lane-aware super-word add). */
    void simd_add_state(const Generator& other) override;
    /** @brief Reset the per-word advance index to 0. */
    void simd_reset_word_index() override { pword_idx_ = 0; }

protected:
    /** @brief Default test method — `simd_notprimitive` for SFMT. */
    std::optional<std::string> compute_default_test_method(const std::string& test_type) const override;

private:
    int mexp_;
    int N_;
    int pos1_, sl1_, sl2_, sr1_, sr2_;
    uint32_t msk1_, msk2_, msk3_, msk4_;
    uint32_t parity0_, parity1_, parity2_, parity3_;

    // `state_` (inherited) holds exactly k = 128·N bits — the true
    // F2-linear state of the SFMTGen recurrence.  N = MEXP/128 + 1, so
    // for SFMT607 this is 640 = 5·128 = 10·64 bits, not 607.
    BitVect work_;        // 128 · N bits, lane-addressable view of state_
    int buf_idx_;         // flat 4·N 32-bit lane pointer
    uint32_t last_output_;
    int pword_idx_ = 0;   // SIMD-PIS per-word advance index, matching
                          // MTToolBox sfmtsearch.hpp's `index` field.

    void generate_all();
    void period_certify();            // apply SFMTGen parity fix to work_
    void rebuild_work_from_state();   // state_ → work_ + period cert
    void sync_state_from_work();      // work_[0..MEXP-1] → state_

    // Lane access into the 128·N-bit work buffer.  SFMTGen treats u[0]
    // as the low 32 bits of a 128-bit word.  BitVect is MSB-first, so
    // the first 32-bit slot in each 128-bit chunk is the HIGH 32 bits;
    // map lane 0 → slot 3 to align u[0] with the low bits the
    // lshift128 / rshift128 helpers expect.
    uint32_t lane(int word_idx, int l) const {
        return (uint32_t)work_.get_word(4 * word_idx + (3 - l), 32);
    }
    void set_lane(int word_idx, int l, uint32_t val) {
        work_.set_word(4 * word_idx + (3 - l), 32, (uint64_t)val);
    }
};

}  // namespace regpoly::core
