#pragma once
#include "generator.h"
#include "param_spec.h"
#include "params.h"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

// dSFMT (Saito 2009) — double-precision SIMD-oriented Fast Mersenne Twister.
//
// State: N+1 super-words of 128 bits each, where N = (mexp-128)/104 + 1.
// The first N super-words hold status[0..N-1]; status[N] is the "lung"
// register, modified in-place by every do_recursion call (lung carries
// state across recursion steps in a way SFMTGen's stateless update does not).
//
// Lanes per super-word: 2 × uint64_t (lane 0 = u[0], lane 1 = u[1]).
// Output bits per lane: 52 (the low 52 bits — high 12 are reserved for
// IEEE-754 exponent in the actual dSFMT, but the F2-linear analysis uses
// the low-52 view directly).
// Shift parameters: only sl1 (sr is hardcoded = 12 in the upstream sample).
// Masks: 2 × uint64_t (msk1 for u[0], msk2 for u[1]).
//
// SIMD-PIS overrides: lane_count = 2, output_phases = 2, plus
// simd_advance_one_word / simd_read_super_word / simd_add_state /
// simd_reset_word_index following the SFMTGen pattern but with the dSFMT
// recurrence and the lung carry.
class DSFMTGen : public Generator {
public:
    DSFMTGen(int mexp, int pos1, int sl1,
             uint64_t msk1, uint64_t msk2, int L);

    static std::unique_ptr<Generator> from_params(
        const Params& params, int L);
    static std::vector<ParamSpec> param_specs();

    std::string name() const override;
    std::string display_str() const override;
    void init(const BitVect& init_bv) override;
    void next() override;
    std::unique_ptr<Generator> copy() const override;
    BitVect get_output() const override;

    int mexp() const { return mexp_; }
    int N() const { return N_; }
    int period_exponent() const override { return mexp_; }

    // dSFMT packs 2 lanes per 128-bit super-word; one full state cycle
    // emits 2 lanes per word over N words.  output_phases = 2 partitions
    // the stream by lane within a state-cycle (analogous to SFMTGen's
    // 128/L = 4 phases at L=32).
    int output_phases() const override { return 2; }

    // SIMD lane count: 2 × 64-bit lanes per 128-bit super-word.
    int simd_lane_count() const override { return 2; }

    // Set state without triggering generate_all (matches SFMTGen semantics
    // — required by SIMD-PIS to XOR basis-vector states without advancing
    // the recurrence).
    void set_raw_state(const BitVect& s) override;

    // ── SIMD-aware overrides for METHOD_SIMD_NOTPRIMITIVE ──────────────
    void simd_advance_one_word() override;
    BitVect simd_read_super_word(int start_mode) const override;
    void simd_add_state(const Generator& other) override;
    void simd_reset_word_index() override { pword_idx_ = 0; }

    // dSFMT's F2-linear state has size 128*(N+1) bits (status[0..N-1]
    // plus lung), but one full generate_all batch is N do_recursion
    // calls — the lung is updated in-place by every do_recursion and
    // is not a separate "word" in the per-word advance loop.
    int simd_full_step_words() const override { return N_; }

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
