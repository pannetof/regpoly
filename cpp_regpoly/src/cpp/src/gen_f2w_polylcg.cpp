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

// ── Factory methods ────────────────────────────────────────────────────

std::unique_ptr<Generateur> GenF2wPolyLCG::from_params(const Params& params, int L) {
    int w = (int)params.get_int("w");
    int r = (int)params.get_int("r");
    uint64_t modM = (uint64_t)params.get_int("modM");
    bool normal_basis = params.get_bool("normal_basis", false);
    int step_count = (int)params.get_int("step", 1);
    auto nocoeff_vals = params.get_int_vec("nocoeff");
    auto coeff_vals = params.get_uint_vec("coeff");
    int nbcoeff = (int)nocoeff_vals.size();
    return std::make_unique<GenF2wPolyLCG>(
        w, r, nbcoeff, nocoeff_vals, coeff_vals,
        modM, normal_basis, step_count, L);
}

std::vector<ParamSpec> GenF2wPolyLCG::param_specs() {
    return {
        {"w",            "int",      true,  false, 0, "",             "", false},
        {"r",            "int",      true,  false, 0, "",             "", false},
        {"modM",         "int",      true,  false, 0, "",             "", false},
        {"normal_basis", "bool",     true,  true,  0, "",             "", false},
        {"step",         "int",      true,  true,  1, "",             "", false},
        {"nocoeff",      "int_vec",  true,  false, 0, "",             "", false},
        {"coeff",        "uint_vec", false, false, 0, "bitmask_vec",  "w,nocoeff", false},
    };
}
