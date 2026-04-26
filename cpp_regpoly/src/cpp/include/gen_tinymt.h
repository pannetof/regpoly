#pragma once
#include "generateur.h"
#include "param_spec.h"
#include "params.h"
#include <cstdint>
#include <memory>
#include <string>

// TinyMT32 (Saito-Matsumoto 2011) — 127-bit-state Mersenne-twister-style
// generator with one tempering parameter.
//
// Reference: MTToolBox samples/TinyMTDC/tinymt32search.hpp.
//
// State layout in `state_` (BitVect, MSB-first; status[0] high bit is
// not part of the F2-linear state per the upstream isZero() mask
// 0x7fffffff, so the effective recurrence dimension is 127 even though
// the BitVect carries 128 bits):
//     bits   0.. 31 = status[0]
//     bits  32.. 63 = status[1]
//     bits  64.. 95 = status[2]
//     bits  96..127 = status[3]
//
// Recurrence (from tinymt32search.hpp::next_state):
//     y = status[3]
//     x = (status[0] & 0x7fffffff) ^ status[1] ^ status[2]
//     x ^= (x << 1)
//     y ^= (y >> 1) ^ x
//     status[0] = status[1]
//     status[1] = status[2]
//     status[2] = x ^ (y << 10)
//     status[3] = y
//     if (y & 1) {
//         status[1] ^= mat1
//         status[2] ^= mat2
//     }
//
// Tempering (from tinymt32search.hpp::temper):
//     t1 = status[0] ^ (status[2] >> 8)
//     out = status[3] ^ t1
//     if (t1 & 1) out ^= tmat
//     return out
class TinyMT32 : public Generateur {
public:
    TinyMT32(uint32_t mat1, uint32_t mat2, uint32_t tmat, int L);

    static std::unique_ptr<Generateur> from_params(const Params& params, int L);
    static std::vector<ParamSpec> param_specs();

    std::string name() const override;
    std::string display_str() const override;
    void init(const BitVect& init_bv) override;
    void next() override;
    std::unique_ptr<Generateur> copy() const override;
    BitVect get_output() const override;

    uint32_t mat1() const { return mat1_; }
    uint32_t mat2() const { return mat2_; }
    uint32_t tmat() const { return tmat_; }

    int period_exponent() const override { return 127; }

private:
    uint32_t mat1_, mat2_, tmat_;
    uint32_t last_output_;

    void load_status(uint32_t s[4]) const;
    void store_status(const uint32_t s[4]);
    uint32_t compute_temper(const uint32_t s[4]) const;
};
