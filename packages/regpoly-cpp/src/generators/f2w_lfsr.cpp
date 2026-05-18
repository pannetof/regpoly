// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "f2w_lfsr.h"
#include "f2w_polylcg.h"

using namespace regpoly::core;


namespace regpoly::core {

F2wLFSRGen::F2wLFSRGen(int w, int r, int nbcoeff,
                         const std::vector<int>& nocoeff,
                         const std::vector<uint64_t>& coeff,
                         uint64_t modM, bool normal_basis,
                         int step_count, int L)
    : F2wBaseGen(w, r, nbcoeff, nocoeff, coeff, modM, normal_basis, step_count, L)
{
    int k_bits = w * r;
    if (k_bits < 64) {
        L_ = ((L - k_bits) / w) * w + k_bits;
    }
    state_bits_ = std::max(L_, k_bits);
}

void F2wLFSRGen::init(const BitVect& init_bv) {
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

void F2wLFSRGen::next() {
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

void F2wLFSRGen::get_transition_state(uint64_t* out_words, int out_nwords) const {
    BitVect tmp(k_);
    tmp.copy_part_from(state_, k_);
    int n = std::min(out_nwords, tmp.nwords());
    for (int i = 0; i < n; i++)
        out_words[i] = tmp.data()[i];
    for (int i = n; i < out_nwords; i++)
        out_words[i] = 0;
}

std::unique_ptr<Generator> F2wLFSRGen::copy() const {
    auto p = std::make_unique<F2wLFSRGen>(
        w_, r_, nbcoeff_, nocoeff_, coeff_,
        modM_, normal_basis_, step_count_, L_);
    p->state_ = state_.copy();
    p->table_ = table_;
    p->state_bits_ = state_bits_;
    return p;
}

// ── Factory methods ────────────────────────────────────────────────────

std::unique_ptr<Generator> F2wLFSRGen::from_params(const Params& params, int L) {
    int w = (int)params.get_int("w");
    int r = (int)params.get_int("r");
    uint64_t modM = (uint64_t)params.get_int("modM");
    bool normal_basis = params.get_bool("normal_basis", false);
    int step_count = (int)params.get_int("step", 1);

    // Mirror F2wPolyLCGGen::from_params: accept either an explicit
    // nocoeff vector (catalog / direct C++ path) or derive it from
    // {nb_terms, t, q} (web paper-notation search form path).
    std::vector<int> nocoeff_vals;
    if (params.has("nocoeff") && !params.get_int_vec("nocoeff").empty()) {
        nocoeff_vals = params.get_int_vec("nocoeff");
    } else {
        int nb_terms = (int)params.get_int("nb_terms", 3);
        if (nb_terms != 2 && nb_terms != 3)
            throw std::invalid_argument(
                "F2w from_params: nb_terms must be 2 or 3 (got "
                + std::to_string(nb_terms) + ")");
        int t = (int)params.get_int("t");
        if (t < 1 || t > r - 1)
            throw std::invalid_argument(
                "F2w from_params: t must be in [1, r-1]");
        if (nb_terms == 3) {
            int q = (int)params.get_int("q");
            if (q < 1 || q >= t)
                throw std::invalid_argument(
                    "F2w from_params: q must be in [1, t-1]");
            nocoeff_vals = {t, q, 0};
        } else {
            nocoeff_vals = {t, 0};
        }
    }
    auto coeff_vals = params.get_uint_vec("coeff");
    int nbcoeff = (int)nocoeff_vals.size();
    return std::make_unique<F2wLFSRGen>(
        w, r, nbcoeff, nocoeff_vals, coeff_vals,
        modM, normal_basis, step_count, L);
}

std::vector<ParamSpec> F2wLFSRGen::param_specs() {
    return F2wPolyLCGGen::param_specs();  // same parameter layout
}

}  // namespace regpoly::core
