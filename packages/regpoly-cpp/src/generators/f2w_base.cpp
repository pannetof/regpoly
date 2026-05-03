#include "f2w_base.h"
#include <sstream>
#include <iomanip>
#include <cstdio>

F2wBaseGen::F2wBaseGen(int w, int r, int nbcoeff,
                       const std::vector<int>& nocoeff,
                       const std::vector<uint64_t>& coeff,
                       uint64_t modM, bool normal_basis,
                       int step_count, int L)
    : Generator(w * r, std::min(L, w * r)),
      w_(w), r_(r), nbcoeff_(nbcoeff),
      nocoeff_(nocoeff), coeff_(coeff),
      modM_(modM), normal_basis_(normal_basis),
      maskw_((w == 64) ? ~0ULL : ((1ULL << w) - 1)),
      step_count_(step_count)
{
    if (normal_basis)
        table_ = gf2w_make_table(w, modM);
}

std::string F2wBaseGen::display_str() const {
    std::ostringstream oss;

    // Line 1: type_name()
    oss << type_name() << "\n";

    // Line 2: "z^{r_} + " then coefficients then " -- coeffs in F_2/P(z) where P(z)={modM:08x}"
    oss << "z^" << r_ << " + ";
    for (int j = 0; j < nbcoeff_; j++) {
        char coeff_hex[16];
        snprintf(coeff_hex, sizeof(coeff_hex), "%08x", (unsigned)coeff_[j]);
        if (j < nbcoeff_ - 1) {
            oss << "(" << coeff_hex << ")z^" << nocoeff_[j] << " +";
        } else {
            // Last coefficient
            if (nocoeff_[j] == 0) {
                oss << "(" << coeff_hex << ")";
            } else {
                oss << "(" << coeff_hex << ")z^" << nocoeff_[j];
            }
        }
    }
    char modM_hex[16];
    snprintf(modM_hex, sizeof(modM_hex), "%08x", (unsigned)modM_);
    oss << " -- coeffs in F_2/P(z) where P(z)=" << modM_hex << "\n";

    // Line 3: "Step = {step_count_}"
    oss << "Step = " << step_count_;

    return oss.str();
}

uint64_t F2wBaseGen::multiply(uint64_t a, uint64_t b) const {
    if (normal_basis_)
        return gf2w_multiply_normal(a, b, w_, table_.data());
    else
        return gf2w_multiply_poly(a, b, w_, modM_);
}

uint64_t F2wBaseGen::V(int idx) const {
    return state_.get_word(idx, w_);
}

void F2wBaseGen::SetV(int idx, uint64_t val) {
    state_.set_word(idx, w_, val & maskw_);
}
