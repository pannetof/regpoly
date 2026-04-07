#include "trans_permutation.h"
#include <sstream>

PermutationTrans::PermutationTrans(int w, int p, int q)
    : p_(p), q_(q)
{
    w_ = w;
}

std::string PermutationTrans::name() const {
    return "Permutation";
}

std::string PermutationTrans::display_str() const {
    std::ostringstream oss;
    oss << "Permutation(" << p_ << "," << q_ << ")";
    return oss.str();
}

void PermutationTrans::apply(BitVect& state) const {
    int n = state.nbits();
    BitVect permuted(n);
    for (int i = 0; i < w_; i++) {
        int u = (i * p_ + q_) % w_;
        if (state.get_bit(u))
            permuted.set_bit(i, 1);
    }
    state.and_invmask(w_);
    state.xor_with(permuted);
}

std::unique_ptr<Transformation> PermutationTrans::copy() const {
    return std::make_unique<PermutationTrans>(w_, p_, q_);
}

void PermutationTrans::update(const Params& params) {
    w_ = (int)params.get_int("w", w_);
    p_ = (int)params.get_int("p", p_);
    q_ = (int)params.get_int("q", q_);
}
