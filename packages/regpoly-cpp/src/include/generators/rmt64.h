#pragma once
#include "generator.h"
#include "param_spec.h"
#include "params.h"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

// RMT64Gen — 64-bit Reducible Mersenne Twister (Saito 2011 et seq).
//
// Reference: MTToolBox samples/RMT/rmt64.hpp.  The "Reducible" name
// signals that the characteristic polynomial of the F2-linear state
// is *not* primitive in general — period is 2^mexp - 1 within an
// invariant subspace of the (mexp/64+1)*64 = state_size F2-linear
// state.  Equidistribution uses METHOD_NOTPRIMITIVE.
//
// State layout in `state_` (BitVect, MSB-first):
//     bits   0..  63  = status[0]
//     bits  64.. 127  = status[1]
//      ...
//     bits  64*(N-1)..64*N-1 = status[N-1]
//   where N = mexp/64 + 1.
//
// We do not roll the index into the BitVect; instead, like
// MTToolBox, we keep an `index_` integer that points to the
// most-recently-updated word.
//
// Recurrence (from rmt64.hpp::generate, sh1=29, sh2=17, sh3=37, sh4=43):
//     index = (index + 1) % size
//     x = state[index]
//     y = state[(index + pos) % size]
//     y ^= (y << 17)
//     x = y ^ (x >> 1) ^ ((x & 1) ? mata : 0)
//     state[index] = x
//     // tempering:
//     x ^= (x >> 29) & 0x5555_5555_5555_5555
//     x ^= (x << 17) & maskb
//     x ^= (x << 37) & maskc
//     x ^= (x >> 43)
//     return x      // 64-bit output
class RMT64Gen : public Generator {
public:
    RMT64Gen(int mexp, int pos, uint64_t mata,
          uint64_t maskb, uint64_t maskc, int L);

    static std::unique_ptr<Generator> from_params(const Params& params, int L);
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

private:
    int mexp_;
    int N_;
    int pos_;
    uint64_t mata_, maskb_, maskc_;
    int index_;          // most-recently-updated word
    uint64_t last_output_;

    uint64_t get_word_u64(int i) const {
        return state_.get_word(i, 64);
    }
    void set_word_u64(int i, uint64_t v) {
        state_.set_word(i, 64, v);
    }
};
