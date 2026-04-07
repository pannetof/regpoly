#include "gen_mt.h"
#include <algorithm>
#include <cstdio>

MersenneTwister::MersenneTwister(int w, int r, int m, int p, uint64_t a, int L)
    : Generateur(w * r - p, L),
      w_(w), r_(r), m_(m), p_(p), a_(a),
      i_(0), uu_(0), ll_(0), maskw_(0),
      state_bits_(w * r)
{
    state_ = BitVect(state_bits_);
}

std::string MersenneTwister::name() const { return "Mersenne Twister"; }

std::string MersenneTwister::display_str() const {
    char buf[128];
    snprintf(buf, sizeof(buf),
             " w= %3d  r= %3d  m= %3d  p= %3d a= %10x  ",
             w_, r_, m_, p_, (unsigned)a_);
    return std::string(buf);
}

void MersenneTwister::init(const BitVect& init_bv) {
    state_ = BitVect(state_bits_);
    state_.copy_part_from(init_bv, k_);

    i_ = 0;
    maskw_ = (w_ == 64) ? ~0ULL : ((1ULL << w_) - 1);
    uu_ = maskw_ << p_;

    if (p_ != 0)
        ll_ = (1ULL << p_) - 1;
    else
        ll_ = 0;
}

void MersenneTwister::next() {
    uint64_t Y = (V(i_ % r_) & uu_) | (V((i_ + 1) % r_) & ll_);

    uint64_t V_i;
    if (Y & 1ULL)
        V_i = V((i_ + m_) % r_) ^ (Y >> 1) ^ a_;
    else
        V_i = V((i_ + m_) % r_) ^ (Y >> 1);

    SetV(i_, V_i);
    i_ = (i_ + 1) % r_;
}

BitVect MersenneTwister::rotated_state() const {
    int total = state_bits_;
    int rotation = (i_ * w_) % total;

    if (rotation == 0)
        return state_.copy();

    BitVect part1 = state_.copy();
    part1.rshift(rotation);
    BitVect part2 = state_.copy();
    part2.lshift(total - rotation);
    part1.xor_with(part2);
    return part1;
}

BitVect MersenneTwister::get_output() const {
    BitVect rot = rotated_state();
    BitVect out(L_);
    out.copy_part_from(rot, L_);
    return out;
}

void MersenneTwister::get_transition_state(uint64_t* out_words, int out_nwords) const {
    BitVect rot = rotated_state();
    BitVect tmp(k_);
    tmp.copy_part_from(rot, k_);
    int n = std::min(out_nwords, tmp.nwords());
    for (int i = 0; i < n; i++)
        out_words[i] = tmp.data()[i];
    for (int i = n; i < out_nwords; i++)
        out_words[i] = 0;
}

std::unique_ptr<Generateur> MersenneTwister::copy() const {
    auto p = std::make_unique<MersenneTwister>(w_, r_, m_, p_, a_, L_);
    p->state_ = state_.copy();
    p->i_ = i_;
    p->uu_ = uu_;
    p->ll_ = ll_;
    p->maskw_ = maskw_;
    p->state_bits_ = state_bits_;
    return p;
}

uint64_t MersenneTwister::V(int idx) const {
    int ii = r_ - idx - 1;
    return state_.get_word(ii, w_);
}

void MersenneTwister::SetV(int idx, uint64_t val) {
    int ii = r_ - idx - 1;
    state_.set_word(ii, w_, val & maskw_);
}
