#include "gen_f2w_lfsr.h"

GenF2wLFSR::GenF2wLFSR(int w, int r, int nbcoeff,
                         const std::vector<int>& nocoeff,
                         const std::vector<uint64_t>& coeff,
                         uint64_t modM, bool normal_basis,
                         int step_count, int L)
    : GenF2wBase(w, r, nbcoeff, nocoeff, coeff, modM, normal_basis, step_count, L)
{
    int k_bits = w * r;
    if (k_bits < 64) {
        L_ = ((L - k_bits) / w) * w + k_bits;
    }
    state_bits_ = std::max(L_, k_bits);
}

void GenF2wLFSR::init(const BitVect& init_bv) {
    state_ = BitVect(state_bits_);
    state_.copy_part_from(init_bv, k_);

    int p = L_ / w_ - r_;
    for (int m = 0; m < p; m++) {
        uint64_t res = 0;
        for (int j = 0; j < nbcoeff_; j++)
            res ^= multiply(V(nocoeff_[j] + m), coeff_[j]);
        SetV(r_ + m, res);
    }
}

void GenF2wLFSR::next() {
    for (int i = 0; i < step_count_; i++) {
        uint64_t res = 0;
        for (int j = 0; j < nbcoeff_; j++)
            res ^= multiply(V(nocoeff_[j]), coeff_[j]);
        state_.lshift(w_);
        SetV(r_ - 1, res);
    }

    int p = L_ / w_ - r_;
    for (int m = 0; m < p; m++) {
        uint64_t res = 0;
        for (int j = 0; j < nbcoeff_; j++)
            res ^= multiply(V(nocoeff_[j] + m), coeff_[j]);
        SetV(r_ + m, res);
    }
}

void GenF2wLFSR::get_transition_state(uint64_t* out_words, int out_nwords) const {
    BitVect tmp(k_);
    tmp.copy_part_from(state_, k_);
    int n = std::min(out_nwords, tmp.nwords());
    for (int i = 0; i < n; i++)
        out_words[i] = tmp.data()[i];
    for (int i = n; i < out_nwords; i++)
        out_words[i] = 0;
}

std::unique_ptr<Generateur> GenF2wLFSR::copy() const {
    auto p = std::make_unique<GenF2wLFSR>(
        w_, r_, nbcoeff_, nocoeff_, coeff_,
        modM_, normal_basis_, step_count_, L_);
    p->state_ = state_.copy();
    p->table_ = table_;
    p->state_bits_ = state_bits_;
    return p;
}
