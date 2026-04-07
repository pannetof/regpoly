#include "gen_marsaxorshift.h"
#include <cstdio>
#include <sstream>

MarsaXorshiftGen::MarsaXorshiftGen(
    int type, int w, int r, int m,
    const Type1Params& t1,
    const Type2xParams& t2x,
    const std::vector<Tap>& taps,
    const Type4Params& t4,
    const std::vector<MiEntry>& mi,
    int L)
    : Generateur(w * r, L),
      type_(type), w_(w), r_(r), m_(m),
      t1_(t1), t2x_(t2x), taps_(taps), t4_(t4), mi_(mi)
{
    state_ = BitVect(k_);
}

std::string MarsaXorshiftGen::name() const { return "Marsaglia Xor-shift"; }

std::string MarsaXorshiftGen::display_str() const {
    std::ostringstream oss;
    oss << "type = " << type_ << ", w = " << w_ << ", r = " << r_;

    if (type_ == 1) {
        oss << ", a = " << t1_.a << ", b = " << t1_.b << ", c = " << t1_.c;
    } else if (type_ >= 21 && type_ <= 25) {
        oss << ", m = " << m_;
        oss << ", p = [";
        for (size_t i = 0; i < t2x_.p.size(); i++) {
            if (i > 0) oss << ", ";
            oss << t2x_.p[i];
        }
        oss << "], q = [";
        for (size_t i = 0; i < t2x_.q.size(); i++) {
            if (i > 0) oss << ", ";
            oss << t2x_.q[i];
        }
        oss << "]";
    } else if (type_ == 3) {
        oss << ", taps = [";
        for (size_t i = 0; i < taps_.size(); i++) {
            if (i > 0) oss << ", ";
            oss << "(" << taps_[i].position << ", " << taps_[i].shift << ")";
        }
        oss << "]";
    } else if (type_ == 4) {
        oss << ", m = " << m_;
        oss << ", p = [" << t4_.p[0] << ", " << t4_.p[1] << "]";
        oss << ", q = [" << t4_.q[0] << ", " << t4_.q[1] << "]";
    } else if (type_ == 100) {
        oss << ", mi = [";
        for (size_t i = 0; i < mi_.size(); i++) {
            if (i > 0) oss << ", ";
            oss << "{pos=" << mi_[i].position << ", shifts=[";
            for (size_t j = 0; j < mi_[i].shifts.size(); j++) {
                if (j > 0) oss << ", ";
                oss << mi_[i].shifts[j];
            }
            oss << "]}";
        }
        oss << "]";
    }
    return oss.str();
}

void MarsaXorshiftGen::init(const BitVect& init_bv) {
    state_ = BitVect(k_);
    state_.copy_part_from(init_bv, k_);
}

uint64_t MarsaXorshiftGen::ShiftR(uint64_t v, int s) {
    if (s > 0) return v >> s;
    return v << (-s);
}

uint64_t MarsaXorshiftGen::V(int idx) const {
    return (uint64_t)state_.get_word(idx, 32);
}

void MarsaXorshiftGen::SetV(int idx, uint64_t val) {
    state_.set_word(idx, 32, val & 0xFFFFFFFF);
}

void MarsaXorshiftGen::next() {
    if (type_ == 1) {
        // Single word xorshift
        uint64_t v = V(0);
        v ^= ShiftR(v, t1_.a);
        v ^= ShiftR(v, t1_.b);
        v ^= ShiftR(v, t1_.c);
        SetV(0, v);

    } else if (type_ >= 21 && type_ <= 25) {
        // Two-component, r words — uses 64-bit arithmetic to match C's ulong behavior.
        // Left shifts create upper bits that right shifts pull back down.
        uint64_t y = V(m_ - 1);
        uint64_t x = V(r_ - 1);
        state_.rshift(w_);
        for (int j = 0; j < 3; j++) {
            if (t2x_.p[j] != 0) x ^= ShiftR(x, t2x_.p[j]);
        }
        for (int j = 0; j < 3; j++) {
            if (t2x_.q[j] != 0) y ^= ShiftR(y, t2x_.q[j]);
        }
        SetV(0, x ^ y);

    } else if (type_ == 3) {
        // Multi-tap — uses 64-bit arithmetic to match C's ulong behavior.
        // C reads *(ulong*)&vect[i] which on little-endian gives
        // vect[i] | (vect[i+1] << 32).
        uint64_t t = 0;
        for (const auto& tap : taps_) {
            int idx = tap.position - 1;
            uint64_t v64 = V(idx);
            if (idx + 1 < r_)
                v64 |= V(idx + 1) << 32;
            uint64_t v32 = V(idx);  // 32-bit operand in C's XOR
            t ^= v32 ^ ShiftR(v64, tap.shift);
        }
        state_.rshift(w_);
        SetV(0, t);

    } else if (type_ == 4) {
        // Two-component, 2 shifts each — uses 64-bit arithmetic to match C's ulong behavior.
        uint64_t y = V(m_ - 1);
        uint64_t x = V(r_ - 1);
        state_.rshift(w_);
        x ^= ShiftR(x, t4_.p[0]);
        x ^= ShiftR(x, t4_.p[1]);
        y ^= ShiftR(y, t4_.q[0]);
        y ^= ShiftR(y, t4_.q[1]);
        SetV(0, x ^ y);

    } else if (type_ == 100) {
        // General — uses 64-bit arithmetic to match C's ulong behavior.
        uint64_t t = 0;
        for (const auto& entry : mi_) {
            uint64_t temp = V(entry.position - 1);
            for (int shift : entry.shifts) {
                temp ^= ShiftR(temp, shift);
            }
            t ^= temp;
        }
        state_.rshift(w_);
        SetV(0, t);
    }
}

std::unique_ptr<Generateur> MarsaXorshiftGen::copy() const {
    auto g = std::make_unique<MarsaXorshiftGen>(
        type_, w_, r_, m_, t1_, t2x_, taps_, t4_, mi_, L_);
    g->state_ = state_.copy();
    return g;
}
