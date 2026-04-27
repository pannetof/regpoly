#pragma once
#include "generator.h"
#include "param_spec.h"
#include "params.h"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

// Mersenne Twister for Graphic Processors (Saito & Matsumoto 2013).
//
// Parameters (per published configuration — MTGP11213, MTGP23209,
// MTGP44497):
//   mexp (structural)   Mersenne exponent; state length N = mexp/32 + 1.
//   pos  (structural)   pick-up offset in the recursion.
//   sh1, sh2            shift amounts in the recursion.
//   mask                upper-bit mask applied to status[i] — holds
//                       the 32 - r high bits where r = 32*N - mexp.
//   tbl[16]             recursion lookup table (low 4 bits of y).
//   tmp_tbl[16]         tempering lookup table (low 4 bits of t).
//
// Recurrence (reference: mtgp-dc `mtgp32_recursion`):
//   y = (x1 & mask) ^ x2 ^ (y_in << sh1)
//   MAT = tbl[y & 0x0f]
//   status_new = y ^ (y >> sh2) ^ MAT
// where x1 = state[i], x2 = state[(i+1) % N], y_in = state[(i+pos) % N].
//
// Tempering (reference: mtgp32_temper):
//   t' = t ^ (t >> 16)
//   output = status_new ^ tmp_tbl[t' & 0x0f]
// where t = state[(i + pos - 1) % N] (the word just before the pick-up).
class MTGPGen : public Generator {
public:
    MTGPGen(int mexp, int pos, int sh1, int sh2, uint32_t mask,
         const std::vector<uint32_t>& tbl,
         const std::vector<uint32_t>& tmp_tbl,
         int L);

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
