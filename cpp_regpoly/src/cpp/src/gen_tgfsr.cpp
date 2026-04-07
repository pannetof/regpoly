#include "gen_tgfsr.h"
#include <cstdio>

TGFSRGen::TGFSRGen(int w, int r, int m, const BitVect& a, int L)
    : Generateur(w * r, L), w_(w), r_(r), m_(m), a_(a.copy()) {}

std::string TGFSRGen::name() const { return "TGFSR"; }

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

std::unique_ptr<Generateur> TGFSRGen::copy() const {
    auto p = std::make_unique<TGFSRGen>(w_, r_, m_, a_, L_);
    p->state_ = state_.copy();
    return p;
}
