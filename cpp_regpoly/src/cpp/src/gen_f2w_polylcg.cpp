#include "gen_f2w_polylcg.h"

GenF2wPolyLCG::GenF2wPolyLCG(int w, int r, int nbcoeff,
                               const std::vector<int>& nocoeff,
                               const std::vector<uint64_t>& coeff,
                               uint64_t modM, bool normal_basis,
                               int step_count, int L)
    : GenF2wBase(w, r, nbcoeff, nocoeff, coeff, modM, normal_basis, step_count, L) {}

void GenF2wPolyLCG::init(const BitVect& init_bv) {
    state_ = BitVect(k_);
    state_.copy_part_from(init_bv, k_);
}

void GenF2wPolyLCG::next() {
    for (int i = 0; i < step_count_; i++) {
        uint64_t VP = V(r_ - 1);
        state_.rshift(w_);
        if (VP) {
            for (int j = 0; j < nbcoeff_; j++) {
                uint64_t res = multiply(VP, coeff_[j]);
                res ^= V(nocoeff_[j]);
                SetV(nocoeff_[j], res);
            }
        }
    }
}

std::unique_ptr<Generateur> GenF2wPolyLCG::copy() const {
    auto p = std::make_unique<GenF2wPolyLCG>(
        w_, r_, nbcoeff_, nocoeff_, coeff_,
        modM_, normal_basis_, step_count_, L_);
    p->state_ = state_.copy();
    p->table_ = table_;
    return p;
}
