#include "xorshift128.h"
#include <cstdio>

XorShift128Gen::XorShift128Gen(int a, int b, int c, int L, int pattern)
    : Generator(128, L), a_(a), b_(b), c_(c), pattern_(pattern) {}

std::string XorShift128Gen::name() const { return "XorShift128Gen"; }

std::string XorShift128Gen::display_str() const {
    char buf[64];
    std::snprintf(buf, sizeof(buf),
                  " a=%d  b=%d  c=%d  pat=%d", a_, b_, c_, pattern_);
    return std::string(buf);
}

void XorShift128Gen::load_state(uint32_t& x, uint32_t& y,
                             uint32_t& z, uint32_t& w) const {
    w = (uint32_t)state_.get_word(0, 32);
    z = (uint32_t)state_.get_word(1, 32);
    y = (uint32_t)state_.get_word(2, 32);
    x = (uint32_t)state_.get_word(3, 32);
}

void XorShift128Gen::store_state(uint32_t x, uint32_t y,
                              uint32_t z, uint32_t w) {
    state_.set_word(0, 32, (uint64_t)w);
    state_.set_word(1, 32, (uint64_t)z);
    state_.set_word(2, 32, (uint64_t)y);
    state_.set_word(3, 32, (uint64_t)x);
}

void XorShift128Gen::init(const BitVect& init_bv) {
    state_ = BitVect(k_);
    state_.copy_part_from(init_bv, k_);
}

void XorShift128Gen::next() {
    uint32_t x, y, z, w;
    load_state(x, y, z, w);
    uint32_t t;
    if (pattern_ == 0) {
        t = x ^ (x << a_);
        x = y;
        y = z;
        z = w;
        w = (w ^ (w >> c_)) ^ (t ^ (t >> b_));
    } else {
        // pattern == 1: matches MTToolBox xorshift-5.cpp.
        t = x ^ (x << b_);
        t ^= t >> c_;
        x = y;
        y = z;
        z = w;
        w = (w ^ (w << a_)) ^ t;
    }
    store_state(x, y, z, w);
}

std::unique_ptr<Generator> XorShift128Gen::copy() const {
    auto p = std::make_unique<XorShift128Gen>(a_, b_, c_, L_, pattern_);
    p->state_ = state_.copy();
    return p;
}

std::unique_ptr<Generator> XorShift128Gen::from_params(
    const Params& params, int L) {
    int a = (int)params.get_int("a");
    int b = (int)params.get_int("b");
    int c = (int)params.get_int("c");
    int pattern = (int)params.get_int("pattern");
    return std::make_unique<XorShift128Gen>(a, b, c, std::min(L, 32), pattern);
}

std::vector<ParamSpec> XorShift128Gen::param_specs() {
    return {
        {"a",       "int", true, false, 0, "range", "1,31", false},
        {"b",       "int", true, false, 0, "range", "1,31", false},
        {"c",       "int", true, false, 0, "range", "1,31", false},
        {"pattern", "int", false, true, 0, "",      "",     false},
    };
}
