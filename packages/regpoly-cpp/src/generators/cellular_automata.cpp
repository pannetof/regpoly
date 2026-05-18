// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "cellular_automata.h"
#include <algorithm>
#include <sstream>
#include <stdexcept>

CellularAutomataGen::CellularAutomataGen(
    int k, const std::vector<int>& rule150_positions, int s, int L)
    : Generator(k, L),
      s_(s),
      rule150_positions_(rule150_positions),
      rules_(k)
{
    if (k <= 0)
        throw std::invalid_argument("CellularAutomataGen: k must be > 0");
    if (s_ < 1)
        throw std::invalid_argument("CellularAutomataGen: s must be >= 1");
    for (int pos : rule150_positions_) {
        if (pos < 0 || pos >= k)
            throw std::invalid_argument(
                "CellularAutomataGen: rule150_positions out of range");
        rules_.set_bit(pos, 1);
    }
}

std::string CellularAutomataGen::name() const {
    return "CellularAutomataGen";
}

std::string CellularAutomataGen::display_str() const {
    std::ostringstream oss;
    oss << "CA(k=" << k_ << ",s=" << s_ << ",r150=[";
    for (size_t i = 0; i < rule150_positions_.size(); ++i) {
        if (i > 0) oss << ",";
        oss << rule150_positions_[i];
    }
    oss << "])";
    return oss.str();
}

void CellularAutomataGen::init(const BitVect& init_bv) {
    state_ = BitVect(k_);
    state_.copy_part_from(init_bv, k_);
}

void CellularAutomataGen::next_one_step_() {
    BitVect L = state_.copy(); L.lshift(1);   // brings cell i+1 → cell i
    BitVect R = state_.copy(); R.rshift(1);   // brings cell i-1 → cell i
    BitVect M = state_.copy(); M.and_with(rules_);  // middle bit at r150 cells
    state_ = std::move(L);
    state_.xor_with(R);
    state_.xor_with(M);
}

void CellularAutomataGen::next() {
    for (int i = 0; i < s_; ++i)
        next_one_step_();
}

std::unique_ptr<Generator> CellularAutomataGen::copy() const {
    auto g = std::make_unique<CellularAutomataGen>(
        k_, rule150_positions_, s_, L_);
    g->state_ = state_.copy();
    return g;
}

// ── Factory ────────────────────────────────────────────────────────────

std::unique_ptr<Generator> CellularAutomataGen::from_params(
    const Params& params, int L)
{
    int k = (int)params.get_int("k");
    int s = (int)params.get_int("s", 1);
    auto positions = params.get_int_vec("rule150_positions");
    return std::make_unique<CellularAutomataGen>(k, positions, s, L);
}

std::vector<ParamSpec> CellularAutomataGen::param_specs() {
    return {
        {"k",                 "int",     true,  false, 0, "",     "", false},
        {"rule150_positions", "int_vec", true,  false, 0, "none", "", false},
        {"s",                 "int",     true,  true,  1, "",     "", false},
    };
}

std::optional<std::string>
CellularAutomataGen::compute_default_test_method(const std::string& test_type) const {
    if (test_type == "equidistribution")
        return std::string("matricial");
    return Generator::compute_default_test_method(test_type);
}
