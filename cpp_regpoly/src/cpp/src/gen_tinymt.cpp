#include "gen_tinymt.h"
#include <cstdio>

TinyMT32::TinyMT32(uint32_t mat1, uint32_t mat2, uint32_t tmat, int L)
    : Generateur(128, L),
      mat1_(mat1), mat2_(mat2), tmat_(tmat),
      last_output_(0) {}

std::string TinyMT32::name() const { return "TinyMT32"; }

std::string TinyMT32::display_str() const {
    char buf[80];
    std::snprintf(buf, sizeof(buf),
                  " mat1=%08x mat2=%08x tmat=%08x", mat1_, mat2_, tmat_);
    return std::string(buf);
}

void TinyMT32::load_status(uint32_t s[4]) const {
    s[0] = (uint32_t)state_.get_word(0, 32);
    s[1] = (uint32_t)state_.get_word(1, 32);
    s[2] = (uint32_t)state_.get_word(2, 32);
    s[3] = (uint32_t)state_.get_word(3, 32);
}

void TinyMT32::store_status(const uint32_t s[4]) {
    state_.set_word(0, 32, (uint64_t)s[0]);
    state_.set_word(1, 32, (uint64_t)s[1]);
    state_.set_word(2, 32, (uint64_t)s[2]);
    state_.set_word(3, 32, (uint64_t)s[3]);
}

uint32_t TinyMT32::compute_temper(const uint32_t s[4]) const {
    uint32_t t1 = s[0] ^ (s[2] >> 8);
    uint32_t out = s[3] ^ t1;
    if (t1 & 1) out ^= tmat_;
    return out;
}

void TinyMT32::init(const BitVect& init_bv) {
    state_ = BitVect(k_);
    state_.copy_part_from(init_bv, k_);
    uint32_t s[4];
    load_status(s);
    // The all-zero state is a fixed point of the recurrence; the
    // BM-based notprimitive method needs a non-degenerate stream to
    // recover the characteristic polynomial.  Mirror MTToolBox's
    // tinymt32::seed(1): if the supplied init is degenerate, set
    // status[3] = 1.
    if ((s[0] & 0x7fffffffu) == 0 && s[1] == 0 && s[2] == 0 && s[3] == 0) {
        s[3] = 1;
        store_status(s);
    }
    last_output_ = compute_temper(s);
}

void TinyMT32::next() {
    uint32_t s[4];
    load_status(s);
    uint32_t y = s[3];
    uint32_t x = (s[0] & 0x7fffffffu) ^ s[1] ^ s[2];
    x ^= (x << 1);
    y ^= (y >> 1) ^ x;
    uint32_t ns[4];
    ns[0] = s[1];
    ns[1] = s[2];
    ns[2] = x ^ (y << 10);
    ns[3] = y;
    if (y & 1u) {
        ns[1] ^= mat1_;
        ns[2] ^= mat2_;
    }
    store_status(ns);
    last_output_ = compute_temper(ns);
}

std::unique_ptr<Generateur> TinyMT32::copy() const {
    auto p = std::make_unique<TinyMT32>(mat1_, mat2_, tmat_, L_);
    p->state_ = state_.copy();
    p->last_output_ = last_output_;
    return p;
}

BitVect TinyMT32::get_output() const {
    BitVect out(L_);
    int n = std::min(L_, 32);
    for (int i = 0; i < n; i++) {
        if ((last_output_ >> (31 - i)) & 1u)
            out.set_bit(i, 1);
    }
    return out;
}

std::unique_ptr<Generateur> TinyMT32::from_params(
    const Params& params, int L) {
    uint32_t mat1 = (uint32_t)params.get_int("mat1");
    uint32_t mat2 = (uint32_t)params.get_int("mat2");
    uint32_t tmat = (uint32_t)params.get_int("tmat");
    return std::make_unique<TinyMT32>(mat1, mat2, tmat, std::min(L, 32));
}

std::vector<ParamSpec> TinyMT32::param_specs() {
    return {
        {"mat1", "int", true, false, 0, "", "", false},
        {"mat2", "int", true, false, 0, "", "", false},
        {"tmat", "int", true, false, 0, "", "", false},
    };
}
