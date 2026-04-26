#include "gen_rmt.h"
#include <cstdio>

RMT64::RMT64(int mexp, int pos, uint64_t mata,
             uint64_t maskb, uint64_t maskc, int L)
    : Generateur(((mexp / 64) + 1) * 64, L),
      mexp_(mexp),
      N_((mexp / 64) + 1),
      pos_(pos), mata_(mata), maskb_(maskb), maskc_(maskc),
      index_(N_ - 1),
      last_output_(0) {}

std::string RMT64::name() const { return "RMT64"; }

std::string RMT64::display_str() const {
    char buf[160];
    std::snprintf(buf, sizeof(buf),
                  " MEXP=%d N=%d pos=%d mata=%016lx maskb=%016lx maskc=%016lx",
                  mexp_, N_, pos_,
                  (unsigned long)mata_,
                  (unsigned long)maskb_,
                  (unsigned long)maskc_);
    return std::string(buf);
}

void RMT64::init(const BitVect& init_bv) {
    state_ = BitVect(k_);
    state_.copy_part_from(init_bv, k_);
    index_ = N_ - 1;
    // Avoid the all-zero fixed point so the BM-based notprimitive
    // method has a non-degenerate stream to recover the characteristic
    // polynomial from.  Mirrors MTToolBox's seed(1).
    bool all_zero = true;
    for (int i = 0; i < N_; i++) {
        if (get_word_u64(i) != 0) { all_zero = false; break; }
    }
    if (all_zero) {
        set_word_u64(0, 1);
    }
    last_output_ = 0;
}

void RMT64::next() {
    static const uint64_t maska = 0x5555555555555555ULL;
    index_ = (index_ + 1) % N_;
    uint64_t x = get_word_u64(index_);
    uint64_t y = get_word_u64((index_ + pos_) % N_);
    y ^= (y << 17);
    uint64_t mat = (x & 1ULL) ? mata_ : 0ULL;
    x = y ^ (x >> 1) ^ mat;
    set_word_u64(index_, x);
    // Tempering.
    x ^= (x >> 29) & maska;
    x ^= (x << 17) & maskb_;
    x ^= (x << 37) & maskc_;
    x ^= (x >> 43);
    last_output_ = x;
}

std::unique_ptr<Generateur> RMT64::copy() const {
    auto p = std::make_unique<RMT64>(mexp_, pos_, mata_, maskb_, maskc_, L_);
    p->state_ = state_.copy();
    p->index_ = index_;
    p->last_output_ = last_output_;
    return p;
}

BitVect RMT64::get_output() const {
    BitVect out(L_);
    int n = std::min(L_, 64);
    for (int i = 0; i < n; i++) {
        if ((last_output_ >> (63 - i)) & 1ULL)
            out.set_bit(i, 1);
    }
    return out;
}

std::unique_ptr<Generateur> RMT64::from_params(const Params& params, int L) {
    int mexp = (int)params.get_int("mexp");
    int pos  = (int)params.get_int("pos");
    uint64_t mata  = (uint64_t)params.get_int("mata");
    uint64_t maskb = (uint64_t)params.get_int("maskb");
    uint64_t maskc = (uint64_t)params.get_int("maskc");
    return std::make_unique<RMT64>(mexp, pos, mata, maskb, maskc,
                                   std::min(L, 64));
}

std::vector<ParamSpec> RMT64::param_specs() {
    return {
        {"mexp",  "int", true, false, 0, "", "", false},
        {"pos",   "int", true, false, 0, "", "", false},
        {"mata",  "int", true, false, 0, "", "", false},
        {"maskb", "int", true, false, 0, "", "", false},
        {"maskc", "int", true, false, 0, "", "", false},
    };
}
