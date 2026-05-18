// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "xoshiro.h"

#include "gen_enumerator.h"

#include <NTL/ZZ.h>

#include <sstream>
#include <stdexcept>

namespace {

int validate_ctor_(int w, int r, int A, int B) {
    if (w < 1 || w > 64)
        throw std::invalid_argument(
            "XoshiroGen: w must be in [1, 64]");
    if (r != 4 && r != 8)
        throw std::invalid_argument(
            "XoshiroGen: r (word count) must be 4 or 8");
    if (A < 1 || A > w - 1 || B < 1 || B > w - 1)
        throw std::invalid_argument(
            "XoshiroGen: A, B must each lie in [1, w-1]");
    return w;
}

uint64_t compute_wmask_(int w) {
    return (w >= 64) ? ~uint64_t{0} : ((uint64_t{1} << w) - 1);
}

}  // anonymous namespace

XoshiroGen::XoshiroGen(int w, int r, int A, int B, int L)
    : Generator(w * r, L),
      w_(validate_ctor_(w, r, A, B)),
      r_(r),
      A_(A), B_(B),
      wmask_(compute_wmask_(w))
{
    state_ = BitVect(k_);
}

uint64_t XoshiroGen::rotl_(uint64_t x, int rot) const {
    if (rot == 0) return x & wmask_;
    return ((x << rot) | (x >> (w_ - rot))) & wmask_;
}

std::string XoshiroGen::name() const { return "Xoshiro Generator"; }

std::string XoshiroGen::display_str() const {
    std::ostringstream oss;
    oss << "xoshiro" << (w_ * r_)
        << " (w=" << w_ << ", r=" << r_
        << ", A=" << A_ << ", B=" << B_ << ")";
    return oss.str();
}

void XoshiroGen::init(const BitVect& init_bv) {
    state_ = BitVect(k_);
    state_.copy_part_from(init_bv, k_);
}

void XoshiroGen::next_4_() {
    uint64_t s0 = state_.get_word(0, w_);
    uint64_t s1 = state_.get_word(1, w_);
    uint64_t s2 = state_.get_word(2, w_);
    uint64_t s3 = state_.get_word(3, w_);

    const uint64_t t = (s1 << A_) & wmask_;
    s2 ^= s0;
    s3 ^= s1;
    s1 ^= s2;
    s0 ^= s3;
    s2 ^= t;
    s3 = rotl_(s3, B_);

    state_.set_word(0, w_, s0 & wmask_);
    state_.set_word(1, w_, s1 & wmask_);
    state_.set_word(2, w_, s2 & wmask_);
    state_.set_word(3, w_, s3 & wmask_);
}

void XoshiroGen::next_8_() {
    uint64_t s0 = state_.get_word(0, w_);
    uint64_t s1 = state_.get_word(1, w_);
    uint64_t s2 = state_.get_word(2, w_);
    uint64_t s3 = state_.get_word(3, w_);
    uint64_t s4 = state_.get_word(4, w_);
    uint64_t s5 = state_.get_word(5, w_);
    uint64_t s6 = state_.get_word(6, w_);
    uint64_t s7 = state_.get_word(7, w_);

    const uint64_t t = (s1 << A_) & wmask_;
    s2 ^= s0;
    s5 ^= s1;
    s1 ^= s2;
    s7 ^= s3;
    s3 ^= s4;
    s4 ^= s5;
    s0 ^= s6;
    s6 ^= s7;
    s6 ^= t;
    s7 = rotl_(s7, B_);

    state_.set_word(0, w_, s0 & wmask_);
    state_.set_word(1, w_, s1 & wmask_);
    state_.set_word(2, w_, s2 & wmask_);
    state_.set_word(3, w_, s3 & wmask_);
    state_.set_word(4, w_, s4 & wmask_);
    state_.set_word(5, w_, s5 & wmask_);
    state_.set_word(6, w_, s6 & wmask_);
    state_.set_word(7, w_, s7 & wmask_);
}

void XoshiroGen::next() {
    if (r_ == 4) next_4_();
    else         next_8_();
}

std::unique_ptr<Generator> XoshiroGen::copy() const {
    auto p = std::make_unique<XoshiroGen>(w_, r_, A_, B_, L_);
    p->state_ = state_.copy();
    return p;
}

std::unique_ptr<Generator> XoshiroGen::from_params(
    const Params& params, int L)
{
    if (!params.has("w") || !params.has("r") ||
        !params.has("A") || !params.has("B"))
        throw std::invalid_argument(
            "XoshiroGen: requires w, r, A, B");
    int w = (int)params.get_int("w");
    int r = (int)params.get_int("r");
    int A = (int)params.get_int("A");
    int B = (int)params.get_int("B");
    return std::make_unique<XoshiroGen>(w, r, A, B, L);
}

std::vector<ParamSpec> XoshiroGen::param_specs() {
    return {
        {"w", "int", true,  false, 0, "",      "",      false},
        {"r", "int", true,  false, 0, "",      "",      false},
        {"A", "int", false, false, 0, "range", "1,w-1", false},
        {"B", "int", false, false, 0, "range", "1,w-1", false},
    };
}

// ═══════════════════════════════════════════════════════════════════════════
// Exhaustive enumerator: A × B over [1, w-1]^2
// ═══════════════════════════════════════════════════════════════════════════

namespace {

class XoshiroEnumerator : public GenEnumerator {
public:
    XoshiroEnumerator(int w, int r) : w_(w), r_(r) {
        NTL::ZZ axis_size(w_ - 1);
        total_ = axis_size * axis_size;
    }

    std::string size_dec() const override { return zz_to_dec(total_); }

    std::vector<Axis> axes() const override {
        std::vector<Axis> out;
        const std::string size = zz_to_dec(NTL::ZZ((long)(w_ - 1)));
        std::ostringstream desc;
        desc << "shift / rotation in [1, " << (w_ - 1) << "]";
        for (const char* name : {"A", "B"}) {
            Axis ax;
            ax.name = name;
            ax.size_dec = size;
            ax.describe = desc.str();
            out.push_back(std::move(ax));
        }
        return out;
    }

    Params at(const std::string& idx_dec) const override {
        NTL::ZZ idx = parse_zz(idx_dec);
        if (idx < 0 || idx >= total_)
            throw std::out_of_range(
                "XoshiroGen enumerator: idx out of range");
        NTL::ZZ axis_size(w_ - 1);
        std::vector<NTL::ZZ> sizes = {axis_size, axis_size};
        auto digits = mixed_radix_decode(sizes, idx);
        Params out;
        out.set_int("w", w_);
        out.set_int("r", r_);
        out.set_int("A", NTL::to_long(digits[0]) + 1);
        out.set_int("B", NTL::to_long(digits[1]) + 1);
        return out;
    }

private:
    int w_;
    int r_;
    NTL::ZZ total_;
};

}  // anonymous namespace

std::unique_ptr<GenEnumerator> XoshiroGen::make_enumerator(
    const Params& resolved, int /*L*/)
{
    if (!resolved.has("w"))
        throw std::invalid_argument("needs_w");
    if (!resolved.has("r"))
        throw std::invalid_argument("needs_r");
    int w = (int)resolved.get_int("w");
    int r = (int)resolved.get_int("r");
    if (w < 2 || w > 64)
        throw std::invalid_argument("needs_w_in_2_to_64");
    if (r != 4 && r != 8)
        throw std::invalid_argument("needs_r_in_4_or_8");
    return std::make_unique<XoshiroEnumerator>(w, r);
}
