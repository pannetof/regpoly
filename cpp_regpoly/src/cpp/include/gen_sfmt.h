#pragma once
#include "generateur.h"
#include "param_spec.h"
#include "params.h"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

// SIMD-oriented Fast Mersenne Twister (Saito & Matsumoto 2008).
//
// The state is an array of N 128-bit words, where N = MEXP/128 + 1.
// Each step applies the "do_recursion" update to one array slot,
// advancing the ring index; one full sweep of N steps corresponds to
// one "generate_all" block in the reference implementation.
//
// This class stores the state as a BitVect of 128*N bits and exposes
// each 128-bit word as four 32-bit lanes `u[0..3]` (matching the
// reference implementation's little-endian within-word layout).
//
// Parameters are the ten tabulated SFMT constants: the shift /
// rotation amounts POS1, SL1, SL2, SR1, SR2 and the 4-lane bitmask
// MSK1..MSK4.  The parity constants from the reference code are not
// used here — the regpoly pipeline initialises the state to a
// non-zero bit pattern via Combinaison.
class SFMT : public Generateur {
public:
    SFMT(int mexp, int pos1, int sl1, int sl2, int sr1, int sr2,
         uint32_t msk1, uint32_t msk2, uint32_t msk3, uint32_t msk4,
         int L);

    static std::unique_ptr<Generateur> from_params(
        const Params& params, int L);
    static std::vector<ParamSpec> param_specs();

    std::string name() const override;
    std::string display_str() const override;
    void init(const BitVect& init_bv) override;
    void next() override;
    std::unique_ptr<Generateur> copy() const override;
    BitVect get_output() const override;

    int mexp() const { return mexp_; }
    int N() const { return N_; }
    int period_exponent() const override { return mexp_; }

    // Per Saito–Matsumoto 2008, §3.2 Proposition 2: SFMT's k(v) is
    // computed on the augmented automaton S' = S × {0, 1, …, P-1}
    // where P = 128/L is the number of L-bit lanes per 128-bit word
    // (P = 4 for L=32).  Each phase i's sequence reads only lane i of
    // every output 128-bit word, and k(v) = min_i k_i(v).  Notprimitive's
    // slow-path consults this via the abstract Generateur API and
    // partitions the output stream into P phases.
    int output_phases() const override { return 128 / L_; }

    // SIMD lane count for the SIMD-aware notprimitive method.  SFMT
    // packs 128/L lanes per 128-bit super-word (4 for L=32, 2 for L=64).
    int simd_lane_count() const override { return 128 / L_; }

    // Set state without triggering generate_all (unlike init()).  Mirrors
    // the new state into work_ so subsequent next() calls read lanes of
    // the freshly-set state.  Required by the SIMD-PIS reduction which
    // needs to XOR generator states between iterations without advancing
    // the recurrence.
    void set_raw_state(const BitVect& s) override;

    // ── SIMD-aware overrides for METHOD_SIMD_NOTPRIMITIVE (Phase 3) ──
    void simd_advance_one_word() override;
    BitVect simd_read_super_word(int start_mode) const override;
    void simd_add_state(const Generateur& other) override;
    void simd_reset_word_index() override { pword_idx_ = 0; }

private:
    int mexp_;
    int N_;
    int pos1_, sl1_, sl2_, sr1_, sr2_;
    uint32_t msk1_, msk2_, msk3_, msk4_;
    uint32_t parity0_, parity1_, parity2_, parity3_;

    // `state_` (inherited) holds exactly k = 128·N bits — the true
    // F2-linear state of the SFMT recurrence.  N = MEXP/128 + 1, so
    // for SFMT607 this is 640 = 5·128 = 10·64 bits, not 607.
    // (The earlier convention used k = MEXP, which is the period
    // exponent but not the dimension of the linear state space.  The
    // characteristic polynomial of the 128·N-bit recurrence factors
    // as a degree-MEXP primitive factor times a low-degree non-
    // primitive remainder; equidistribution analysis must therefore
    // use the non-primitive method.)  work_ is a parallel buffer of
    // the same size used by the lane-oriented recurrence code.
    BitVect work_;        // 128 · N bits, lane-addressable view of state_
    int buf_idx_;         // flat 4·N 32-bit lane pointer
    uint32_t last_output_;
    int pword_idx_ = 0;   // SIMD-PIS per-word advance index, matching
                          // MTToolBox sfmtsearch.hpp's `index` field.
                          // simd_advance_one_word does (idx+1)%N then
                          // updates state[idx] (matches MTToolBox::next_state).
                          // Independent of buf_idx_/last_output_.

    void generate_all();
    void period_certify();            // apply SFMT parity fix to work_
    void rebuild_work_from_state();   // state_ → work_ + period cert
    void sync_state_from_work();      // work_[0..MEXP-1] → state_

    // Lane access into the 128·N-bit work buffer.  SFMT treats u[0]
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
