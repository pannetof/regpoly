#include "gen_ac1d.h"
#include <sstream>

AC1DGen::AC1DGen(int n, const std::vector<std::vector<int>>& matrix, int L)
    : Generateur(n, L), n_(n), matrix_(matrix)
{
    state_ = BitVect(k_);
}

std::string AC1DGen::name() const { return "AC-1D"; }

std::string AC1DGen::display_str() const {
    std::ostringstream oss;
    oss << "{";
    for (int i = 0; i < n_; i++) {
        if (i > 0) oss << ",";
        oss << "{";
        for (int j = 0; j < n_; j++) {
            if (j > 0) oss << ",";
            oss << matrix_[i][j];
        }
        oss << "}";
    }
    oss << "}";
    return oss.str();
}

void AC1DGen::init(const BitVect& init_bv) {
    state_ = BitVect(k_);
    state_.copy_part_from(init_bv, k_);
}

void AC1DGen::next() {
    BitVect new_state(n_);
    for (int i = 0; i < n_; i++) {
        int bit = 0;
        for (int j = 0; j < n_; j++) {
            if (matrix_[i][j] == 1) {
                bit ^= state_.get_bit(j);
            }
        }
        new_state.set_bit(i, bit);
    }
    state_ = new_state;
}

std::unique_ptr<Generateur> AC1DGen::copy() const {
    auto g = std::make_unique<AC1DGen>(n_, matrix_, L_);
    g->state_ = state_.copy();
    return g;
}

// ── Factory methods ────────────────────────────────────────────────────

std::unique_ptr<Generateur> AC1DGen::from_params(const Params& params, int L) {
    int n = (int)params.get_int("n");
    auto flat = params.get_int_vec("matrix");
    std::vector<std::vector<int>> matrix(n, std::vector<int>(n, 0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (i * n + j < (int)flat.size())
                matrix[i][j] = flat[i * n + j];
    return std::make_unique<AC1DGen>(n, matrix, L);
}

std::vector<ParamSpec> AC1DGen::param_specs() {
    return {
        {"n",      "int",     true,  false, 0, "",     "", false},
        {"matrix", "int_vec", false, false, 0, "none", "", false},
    };
}
