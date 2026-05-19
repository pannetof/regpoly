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
 * @file dsfmt.h
 * @brief Double-precision SIMD-oriented Fast Mersenne Twister (Saito 2009).
 *
 * `DSFMTGen` is the double-precision variant of SFMT (Saito 2009).
 * State: `N + 1` super-words of 128 bits each, where
 * `N = (mexp - 128) / 104 + 1`. The first `N` super-words hold
 * `status[0..N-1]`; `status[N]` is the "lung" register, modified
 * in-place by every `do_recursion` call (the lung carries state
 * across recursion steps in a way SFMTGen's stateless update does
 * not). Catalog file: `docs/library/saito-2009-dsfmt.yaml`.
 *
 *  - Lanes per super-word: 2 × `uint64_t` (lane 0 = `u[0]`, lane 1 = `u[1]`).
 *  - Output bits per lane: 52 (the low 52 bits — the high 12 are
 *    reserved for the IEEE-754 exponent in the actual dSFMT, but
 *    the F_2-linear analysis uses the low-52 view directly).
 *  - Shift parameters: only `sl1` (`sr` is hardcoded = 12 in the
 *    upstream sample).
 *  - Masks: 2 × `uint64_t` (`msk1` for `u[0]`, `msk2` for `u[1]`).
 *
 * SIMD-PIS overrides: `simd_lane_count = 2`, `output_phases = 2`,
 * plus `simd_advance_one_word` / `simd_read_super_word` /
 * `simd_add_state` / `simd_reset_word_index` following the SFMTGen
 * pattern but with the dSFMT recurrence and the lung carry.
 *
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief Double-precision SIMD Fast Mersenne Twister (Saito 2009).
 *
 * Registered as `"DSFMTGen"` (alias `"dSFMTGen"`). The flagship
 * configuration `dSFMT19937` uses `mexp = 19937, pos1 = 117,
 * sl1 = 19, msk1 = 0x000ffafffffffb3f, msk2 = 0x000ffdfffc90fffd`.
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   Params p;
 *   p.set_int("mexp", 19937);
 *   p.set_int("pos1", 117);
 *   p.set_int("sl1", 19);
 *   p.set_int("msk1", 0x000ffafffffffb3f);
 *   p.set_int("msk2", 0x000ffdfffc90fffd);
 *   auto gen = create_generator("DSFMTGen", p, 52);
 * @endcode
 *
 * @see :py:class:`regpoly.core.generator.Generator`
 * @ingroup core
 */
class DSFMTGen : public Generator {
public:
    /**
     * @brief Construct a DSFMTGen with explicit structural parameters.
     *
     * Most callers should go through `create_generator("DSFMTGen", ...)`
     * rather than this constructor directly.
     *
     * @param mexp  Mersenne exponent.
     * @param pos1  Pick-up offset in the recurrence.
     * @param sl1   Left-shift amount.
     * @param msk1  Mask for lane 0 (u[0]).
     * @param msk2  Mask for lane 1 (u[1]).
     * @param L     Output resolution in bits (typically 52).
     */
    DSFMTGen(int mexp, int pos1, int sl1,
             uint64_t msk1, uint64_t msk2, int L);

    /**
     * @brief Build a DSFMTGen from a Params dict (registry factory hook).
     * @param params  Parameter dict with keys mexp, pos1, sl1, msk1, msk2.
     * @param L       Output resolution in bits.
     * @return        A constructed DSFMTGen as a polymorphic Generator pointer.
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

    /** @brief Family display name — returns the canonical "DSFMTGen" string. */
    std::string name() const override;
    /** @brief Human-readable parameter summary. */
    std::string display_str() const override;
    /** @brief Seed the state from the leading `128 * (N + 1)` bits of `init_bv`. */
    void init(const BitVect& init_bv) override;
    /** @brief Apply one dSFMT recurrence step (one lane per call). */
    void next() override;
    /** @brief Deep copy this generator (state included). */
    std::unique_ptr<Generator> copy() const override;
    /** @brief Current low-52-bit output from the active lane, packed in a BitVect. */
    BitVect get_output() const override;

    /** @brief Mersenne exponent. */
    int mexp() const { return mexp_; }
    /** @brief State length in 128-bit super-words (status array, excluding lung). */
    int N() const { return N_; }
    /** @brief Period exponent — equal to `mexp`. */
    int period_exponent() const override { return mexp_; }

    /**
     * @brief Output phases per state-cycle (= 2 lanes per super-word).
     *
     * dSFMT packs 2 lanes per 128-bit super-word; one full state
     * cycle emits 2 lanes per word over N words. `output_phases = 2`
     * partitions the stream by lane within a state-cycle (analogous
     * to SFMTGen's `128 / L = 4` phases at `L = 32`).
     */
    int output_phases() const override { return 2; }

    /** @brief SIMD lane count: 2 × 64-bit lanes per 128-bit super-word. */
    int simd_lane_count() const override { return 2; }

    /**
     * @brief Install `s` as the raw state without triggering `generate_all`.
     *
     * Matches SFMTGen semantics — required by SIMD-PIS to XOR
     * basis-vector states without advancing the recurrence.
     */
    void set_raw_state(const BitVect& s) override;

    /** @brief Advance the dSFMT by one super-word (SIMD-aware override). */
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

    /**
     * @brief Number of words to advance per full state-step.
     *
     * dSFMT's F_2-linear state has size `128 * (N + 1)` bits
     * (`status[0..N-1]` plus the lung), but one full `generate_all`
     * batch is N `do_recursion` calls — the lung is updated in-place
     * by every `do_recursion` and is not a separate "word" in the
     * per-word advance loop.
     */
    int simd_full_step_words() const override { return N_; }

protected:
    /** @brief Default test method — `simd_notprimitive` for dSFMT. */
    std::optional<std::string> compute_default_test_method(const std::string& test_type) const override;

private:
    int mexp_;
    int N_;
    int pos1_;
    int sl1_;
    uint64_t msk1_;
    uint64_t msk2_;

    // state_ (inherited) holds 128 * (N + 1) bits — status[0..N-1] in the
    // first N super-words plus the lung in super-word index N.  work_ is
    // the lane-addressable parallel buffer (kept in sync with state_).
    BitVect work_;
    int buf_idx_;          // flat 2·N 64-bit lane pointer for next()/get_output
    uint64_t last_output_;  // low 52 bits of the most recently emitted lane
    int pword_idx_ = 0;    // SIMD-PIS per-word advance index (matches
                           // MTToolBox dSFMTsearch.hpp's `index` field)

    void generate_all();          // one full N-step batch, mirrors dsfmt_gen_rand_all
    void rebuild_work_from_state(); // state_ → work_ + all-zero guard
    void sync_state_from_work();    // work_ → state_ (identity copy)

    // Single-word lung-carry recursion: status[idx] = do_recursion(
    //   status[idx], status[(idx+pos1) % N], lung).  Reads + writes lung.
    void do_recursion(int idx);

    // Lane access into the 128·(N+1)-bit work buffer.  dSFMT C uses
    // u[0] = LOW 64 bits of a 128-bit word.  In the MSB-first BitVect,
    // the LOW 64 bits of a 128-bit chunk are the SECOND 64-bit slot
    // within that chunk, so map lane 0 → slot 1, lane 1 → slot 0.
    uint64_t lane_u64(int word_idx, int l) const {
        return work_.get_word(2 * word_idx + (1 - l), 64);
    }
    void set_lane_u64(int word_idx, int l, uint64_t val) {
        work_.set_word(2 * word_idx + (1 - l), 64, val);
    }
};

}  // namespace regpoly::core
