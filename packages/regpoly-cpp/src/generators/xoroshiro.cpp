// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "xoroshiro.h"

#include "gen_enumerator.h"

#include <NTL/ZZ.h>

#include <sstream>
#include <stdexcept>

namespace {

int validate_ctor_(int w, int r, int A, int B, int C) {
    if (w < 1 || w > 64)
        throw std::invalid_argument(
            "XoroshiroGen: w must be in [1, 64]");
    if (r < 2)
        throw std::invalid_argument(
            "XoroshiroGen: r (word count) must be >= 2");
    if (A < 1 || A > w - 1 || B < 1 || B > w - 1 || C < 1 || C > w - 1)
        throw std::invalid_argument(
            "XoroshiroGen: A, B, C must each lie in [1, w-1]");
    return w;
}

uint64_t compute_wmask_(int w) {
    return (w >= 64) ? ~uint64_t{0} : ((uint64_t{1} << w) - 1);
}

}  // anonymous namespace

XoroshiroGen::XoroshiroGen(int w, int r, int A, int B, int C, int L)
    : Generator(w * r, L),
      w_(validate_ctor_(w, r, A, B, C)),
      r_(r),
      A_(A), B_(B), C_(C),
      wmask_(compute_wmask_(w))
{
    state_ = BitVect(k_);
}

uint64_t XoroshiroGen::rotl_(uint64_t x, int rot) const {
    if (rot == 0) return x & wmask_;
    return ((x << rot) | (x >> (w_ - rot))) & wmask_;
}

std::string XoroshiroGen::name() const { return "Xoroshiro Generator"; }

std::string XoroshiroGen::display_str() const {
    std::ostringstream oss;
    oss << "xoroshiro" << (w_ * r_)
        << " (w=" << w_ << ", r=" << r_
        << ", A=" << A_ << ", B=" << B_ << ", C=" << C_ << ")";
    return oss.str();
}

void XoroshiroGen::init(const BitVect& init_bv) {
    state_ = BitVect(k_);
    state_.copy_part_from(init_bv, k_);
}

void XoroshiroGen::next() {
    // State layout: position 0 holds the paper's s_0, position r-1 holds
    // s_{r-1}.  After the step, positions 0..r-3 receive the old 1..r-2
    // (logical left-shift by one word), and the final two positions are
    // computed from s_0 and s_{r-1}.
    const uint64_t s0 = state_.get_word(0, w_);
    const uint64_t s_last = state_.get_word(r_ - 1, w_);
    const uint64_t t = (s0 ^ s_last) & wmask_;

    // Shift state left by one word (block 0 falls off; block r-1 becomes
    // the old block r-2 after this).  In BitVect's bit-0-is-MSB
    // convention "left" means toward bit 0, which is what lshift does.
    state_.lshift(w_);

    state_.set_word(r_ - 2, w_,
                    (rotl_(s0, A_) ^ t ^ (t << B_)) & wmask_);
    state_.set_word(r_ - 1, w_, rotl_(t, C_));
}

std::unique_ptr<Generator> XoroshiroGen::copy() const {
    auto p = std::make_unique<XoroshiroGen>(w_, r_, A_, B_, C_, L_);
    p->state_ = state_.copy();
    return p;
}

std::unique_ptr<Generator> XoroshiroGen::from_params(
    const Params& params, int L)
{
    if (!params.has("w") || !params.has("r") ||
        !params.has("A") || !params.has("B") || !params.has("C"))
        throw std::invalid_argument(
            "XoroshiroGen: requires w, r, A, B, C");
    int w = (int)params.get_int("w");
    int r = (int)params.get_int("r");
    int A = (int)params.get_int("A");
    int B = (int)params.get_int("B");
    int C = (int)params.get_int("C");
    return std::make_unique<XoroshiroGen>(w, r, A, B, C, L);
}

std::vector<ParamSpec> XoroshiroGen::param_specs() {
    // w and r are structural (they determine the state size r*w).  A/B/C
    // are the engine triple, each drawn from [1, w-1] in random mode.
    return {
        {"w", "int", true,  false, 0, "",      "",     false},
        {"r", "int", true,  false, 0, "",      "",     false},
        {"A", "int", false, false, 0, "range", "1,w-1", false},
        {"B", "int", false, false, 0, "range", "1,w-1", false},
        {"C", "int", false, false, 0, "range", "1,w-1", false},
    };
}

// ═══════════════════════════════════════════════════════════════════════════
// Exhaustive enumerator: A × B × C over [1, w-1]^3
// ═══════════════════════════════════════════════════════════════════════════

namespace {

class XoroshiroEnumerator : public GenEnumerator {
public:
    XoroshiroEnumerator(int w, int r) : w_(w), r_(r) {
        NTL::ZZ axis_size(w_ - 1);
        total_ = axis_size * axis_size * axis_size;
    }

    std::string size_dec() const override { return zz_to_dec(total_); }

    std::vector<Axis> axes() const override {
        std::vector<Axis> out;
        const std::string size = zz_to_dec(NTL::ZZ((long)(w_ - 1)));
        std::ostringstream desc;
        desc << "shift / rotation in [1, " << (w_ - 1) << "]";
        for (const char* name : {"A", "B", "C"}) {
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
                "XoroshiroGen enumerator: idx out of range");
        NTL::ZZ axis_size(w_ - 1);
        std::vector<NTL::ZZ> sizes = {axis_size, axis_size, axis_size};
        auto digits = mixed_radix_decode(sizes, idx);
        Params out;
        out.set_int("w", w_);
        out.set_int("r", r_);
        out.set_int("A", NTL::to_long(digits[0]) + 1);
        out.set_int("B", NTL::to_long(digits[1]) + 1);
        out.set_int("C", NTL::to_long(digits[2]) + 1);
        return out;
    }

private:
    int w_;
    int r_;
    NTL::ZZ total_;
};

}  // anonymous namespace

std::unique_ptr<GenEnumerator> XoroshiroGen::make_enumerator(
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
    if (r < 2)
        throw std::invalid_argument("needs_r_at_least_2");
    return std::make_unique<XoroshiroEnumerator>(w, r);
}
