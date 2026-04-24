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

private:
    int mexp_;
    int N_;
    int pos1_, sl1_, sl2_, sr1_, sr2_;
    uint32_t msk1_, msk2_, msk3_, msk4_;
    uint32_t parity0_, parity1_, parity2_, parity3_;

    // `state_` (inherited) holds exactly k = MEXP bits — the true
    // state as far as regpoly is concerned.  The SFMT recurrence also
    // needs 31 extra scratch bits to carry out the 128-bit word-level
    // operations; those live here in a private buffer and are
    // reconstructed from `state_` whenever `state_` is externally
    // modified (e.g. by init() or transition_matrix() perturbations).
    BitVect work_;        // 128 · N bits, private scratch
    int buf_idx_;         // flat 4·N 32-bit lane pointer
    uint32_t last_output_;

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
