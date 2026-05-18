// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "f2w_polylcg.h"

using namespace regpoly::core;


namespace regpoly::core {

F2wPolyLCGGen::F2wPolyLCGGen(int w, int r, int nbcoeff,
                               const std::vector<int>& nocoeff,
                               const std::vector<uint64_t>& coeff,
                               uint64_t modM, bool normal_basis,
                               int step_count, int L)
    : F2wBaseGen(w, r, nbcoeff, nocoeff, coeff, modM, normal_basis, step_count, L) {}

void F2wPolyLCGGen::init(const BitVect& init_bv) {
    state_ = BitVect(k_);
    state_.copy_part_from(init_bv, k_);
}

void F2wPolyLCGGen::next() {
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

std::unique_ptr<Generator> F2wPolyLCGGen::copy() const {
    auto p = std::make_unique<F2wPolyLCGGen>(
        w_, r_, nbcoeff_, nocoeff_, coeff_,
        modM_, normal_basis_, step_count_, L_);
    p->state_ = state_.copy();
    p->table_ = table_;
    return p;
}

// ── Factory methods ────────────────────────────────────────────────────

std::unique_ptr<Generator> F2wPolyLCGGen::from_params(const Params& params, int L) {
    int w = (int)params.get_int("w");
    int r = (int)params.get_int("r");
    uint64_t modM = (uint64_t)params.get_int("modM");
    bool normal_basis = params.get_bool("normal_basis", false);
    int step_count = (int)params.get_int("step", 1);

    std::vector<int> nocoeff_vals;
    if (params.has("nocoeff") && !params.get_int_vec("nocoeff").empty()) {
        // Catalog / direct-C++ path — explicit nocoeff list.
        nocoeff_vals = params.get_int_vec("nocoeff");
    } else {
        // Paper-notation path — derive nocoeff from {nb_terms, t, q}.
        // P(z) = z^r + b_{r-t} z^t + (b_{r-q} z^q if nb_terms=3) + b_r
        // ⇒ nocoeff = [t, q, 0] (3-term) or [t, 0] (2-term).
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
    return std::make_unique<F2wPolyLCGGen>(
        w, r, nbcoeff, nocoeff_vals, coeff_vals,
        modM, normal_basis, step_count, L);
}

std::vector<ParamSpec> F2wPolyLCGGen::param_specs() {
    // Two parallel entry paths share these specs:
    //   - Catalog / direct C++: caller supplies {nocoeff, coeff, modM};
    //     `from_params` short-circuits past nb_terms/t/q.
    //   - Web search form (paper notation): caller fixes {w, r, nb_terms},
    //     optionally fixes {t, q}; the search loop samples {t, q, modM,
    //     coeff} per iteration via the rand_types below. `from_params`
    //     then derives `nocoeff = [t, q, 0]` or `[t, 0]` from
    //     {nb_terms, t, q}.
    //
    //   t, q, modM have has_default=false because the search loop
    //   (`primitive_search.cpp`) only invokes the sampler for non-
    //   defaulted specs.  nb_terms has has_default=true with default 3
    //   so the "search 3-term polynomials" workflow doesn't require it
    //   in structural_params — but the catalog path skips the loop
    //   entirely and reads nb_terms only when nocoeff is absent.
    return {
        {"w",            "int",      true,  false, 0,  "",                "",          false},
        {"r",            "int",      true,  false, 0,  "",                "",          false},
        {"nb_terms",     "int",      true,  true,  3,  "",                "",          false},
        {"modM",         "int",      false, false, 0,  "irreducible_gf2", "w",         false},
        {"normal_basis", "bool",     true,  true,  0,  "",                "",          false},
        {"step",         "int",      true,  true,  1,  "",                "",          false},
        {"t",            "int",      false, false, 0,  "range",           "1,r-1",     false},
        {"q",            "int",      false, false, 0,  "range",           "1,t-1",     false},
        {"nocoeff",      "int_vec",  false, true,  0,  "",                "",          false},
        {"coeff",        "uint_vec", false, false, 0,  "bitmask_vec",     "w,nb_terms", false},
    };
}

}  // namespace regpoly::core
