#include "tgfsr.h"
#include <cstdio>
#include <NTL/GF2X.h>

TGFSRGen::TGFSRGen(int w, int r, int m, const BitVect& a, int L)
    : Generator(w * r, L), w_(w), r_(r), m_(m), a_(a.copy()) {}

std::string TGFSRGen::name() const { return "TGFSRGen"; }

std::string TGFSRGen::display_str() const {
    uint32_t a_word;
    if (a_.nbits() >= 32)
        a_word = (uint32_t)a_.get_word(0, 32);
    else
        a_word = (uint32_t)a_.get_word(0, a_.nbits());

    char buf[128];
    snprintf(buf, sizeof(buf),
             " w= %3d  r= %3d  m= %3d  a= %08x",
             w_, r_, m_, a_word);
    return std::string(buf);
}

void TGFSRGen::init(const BitVect& init_bv) {
    state_ = BitVect(k_);
    state_.copy_part_from(init_bv, k_);
}

void TGFSRGen::next() {
    BitVect temp2 = state_.copy();
    temp2.lshift(w_ * (r_ - m_ - 1));
    temp2.and_mask(w_);

    BitVect temp = state_.copy();
    temp.lshift(w_ * (r_ - 1));
    temp.and_mask(w_);

    if (temp.get_bit(w_ - 1) == 0) {
        temp.rshift(1);
        temp.xor_with(temp2);
    } else {
        temp.rshift(1);
        temp.xor_with(temp2);
        temp.xor_with(a_);
    }

    state_.rshift(w_);
    temp.and_mask(w_);
    state_.xor_with(temp);
}

std::unique_ptr<Generator> TGFSRGen::copy() const {
    auto p = std::make_unique<TGFSRGen>(w_, r_, m_, a_, L_);
    p->state_ = state_.copy();
    return p;
}

BitVect TGFSRGen::char_poly() const {
    // CharTGFSR: P(t) = (t^r + t^m)^w + sum_{j: bit j of a} (t^r + t^m)^j
    // where a is the w-bit twist mask (bit 0 = MSB in our BitVect).
    int K = k_;  // w * r

    // Build tntm = t^r + t^m in NTL
    NTL::GF2X tntm;
    NTL::SetCoeff(tntm, r_);
    NTL::SetCoeff(tntm, m_);

    // res = (t^r + t^m)^w
    NTL::GF2X res;
    NTL::power(res, tntm, w_);

    // Add (t^r + t^m)^j for each set bit j of a
    for (int j = 0; j < w_; j++) {
        if (a_.get_bit(j)) {
            NTL::GF2X term;
            NTL::power(term, tntm, j);
            res += term;
        }
    }

    // Convert to BitVect: bit j = coefficient of z^j, no leading term
    BitVect bv(K);
    for (int j = 0; j < K; j++)
        bv.set_bit(j, IsOne(coeff(res, j)) ? 1 : 0);
    return bv;
}

// ── Factory methods ────────────────────────────────────────────────────

std::unique_ptr<Generator> TGFSRGen::from_params(const Params& params, int L) {
    int w = (int)params.get_int("w");
    int r = (int)params.get_int("r");
    int m = (int)params.get_int("m");
    uint64_t a_val = (uint64_t)params.get_int("a");
    int k = w * r;
    BitVect a_bv(k);
    if (k > 32) {
        for (int i = 0; i < 32; i++)
            if ((a_val >> (31 - i)) & 1)
                a_bv.set_bit(i, 1);
    } else {
        for (int i = 0; i < k; i++)
            if ((a_val >> (k - 1 - i)) & 1)
                a_bv.set_bit(i, 1);
    }
    return std::make_unique<TGFSRGen>(w, r, m, a_bv, std::min(w, L));
}

std::vector<ParamSpec> TGFSRGen::param_specs() {
    return {
        {"w", "int", true,  false, 0, "",        "", false},
        {"r", "int", true,  false, 0, "",        "", false},
        {"m", "int", false, false, 0, "range",   "1,r-1", false},
        {"a", "int", false, false, 0, "bitmask", "w", false},
    };
}
