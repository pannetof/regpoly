#include "gen_tausworthe.h"
#include <algorithm>
#include <sstream>
#include <string>
#include <cstdio>

Tausworthe::Tausworthe(int k, const std::vector<int>& Q, int s, bool quicktaus, int L)
    : Generateur(k, L), Q_(Q), NbCoeff_((int)Q.size()),
      s_(s), quicktaus_(quicktaus), gen_kms_(k - s)
{
    state_ = BitVect(L);
}

std::string Tausworthe::name() const { return "Tausworthe Generator"; }

std::string Tausworthe::display_str() const {
    // Build polynomial string: iterate Q_ from index 0 to NbCoeff_-1
    // If Q_[i]==1 print " x +", else print " x^{Q[i]} +"
    // Then append " 1 "
    std::string poly_str;
    for (int i = 0; i < NbCoeff_; i++) {
        if (Q_[i] == 1)
            poly_str += " x +";
        else
            poly_str += " x^" + std::to_string(Q_[i]) + " +";
    }
    poly_str += " 1 ";

    // Right-pad to 40 chars
    while ((int)poly_str.size() < 40)
        poly_str += " ";

    // Append "   s={s_}"
    poly_str += "   s=" + std::to_string(s_);

    return poly_str;
}

void Tausworthe::init(const BitVect& init_bv) {
    state_ = BitVect(L_);
    state_.copy_part_from(init_bv, k_);

    for (int j = k_; j < L_; j++) {
        int bit = 0;
        for (int i = 0; i < NbCoeff_ - 1; i++)
            bit ^= state_.get_bit(j - (k_ - Q_[i]));
        if (bit & 1)
            state_.set_bit(j, 1);
    }
}

void Tausworthe::next() {
    if (quicktaus_)
        next_quick();
    else
        next_general();
}

std::unique_ptr<Generateur> Tausworthe::copy() const {
    auto p = std::make_unique<Tausworthe>(k_, Q_, s_, quicktaus_, L_);
    p->state_ = state_.copy();
    p->gen_kms_ = gen_kms_;
    return p;
}

BitVect Tausworthe::char_poly() const {
    // The characteristic polynomial is z^k + z^Q[0] + z^Q[1] + ... + 1
    // Q_ contains the non-leading, non-constant exponents, plus k itself.
    // Q_ is sorted. Q_.back() == k_.
    // The constant term (z^0) is always present.
    // Return k bits where bit j = coefficient of z^j (no leading term).
    BitVect bv(k_);
    bv.set_bit(0, 1);  // constant term z^0
    for (int i = 0; i < NbCoeff_ - 1; i++)
        bv.set_bit(Q_[i], 1);
    return bv;
}

void Tausworthe::get_transition_state(uint64_t* out_words, int out_nwords) const {
    BitVect tmp(k_);
    tmp.copy_part_from(state_, k_);
    int n = std::min(out_nwords, tmp.nwords());
    for (int i = 0; i < n; i++)
        out_words[i] = tmp.data()[i];
    for (int i = n; i < out_nwords; i++)
        out_words[i] = 0;
}

void Tausworthe::next_quick() {
    BitVect gen_B(L_);
    gen_B.zero();

    for (int j = 1; j < NbCoeff_ - 1; j++) {
        BitVect shifted = state_.copy();
        shifted.lshift(Q_[j]);
        gen_B.xor_with(shifted);
    }
    gen_B.xor_with(state_);
    gen_B.rshift(gen_kms_);

    state_.and_mask(k_);
    state_.lshift(s_);
    state_.xor_with(gen_B);
}

void Tausworthe::next_general() {
    int ss = std::min(s_, L_ - k_);
    int m = s_;

    while (m > 0) {
        state_.lshift(ss);
        for (int j = L_ - ss; j < L_; j++) {
            int bit = 0;
            for (int i = 0; i < NbCoeff_ - 1; i++)
                bit ^= state_.get_bit(j - (k_ - Q_[i]));
            if (bit & 1)
                state_.set_bit(j, 1);
        }
        m -= ss;
        if (m < ss)
            ss = m;
    }
}
