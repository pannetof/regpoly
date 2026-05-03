#include "dsfmt.h"
#include <algorithm>
#include <cstdio>
#include <stdexcept>

// dSFMT recurrence (from dSMFT/dSFMT-common.h:101-113 and
// MTToolBox/samples/dSFMTdc/dSFMTsearch.hpp:254-272):
//
//   t0 = a.u[0];  t1 = a.u[1];  L0 = lung.u[0];  L1 = lung.u[1];
//   lung.u[0] = (t0 << SL1) ^ (L1 >> 32) ^ (L1 << 32) ^ b.u[0];
//   lung.u[1] = (t1 << SL1) ^ (L0 >> 32) ^ (L0 << 32) ^ b.u[1];
//   r.u[0]    = (lung.u[0] >> SR)  ^ (lung.u[0] & MSK1) ^ t0;
//   r.u[1]    = (lung.u[1] >> SR)  ^ (lung.u[1] & MSK2) ^ t1;
//
// where SR = 12 is hardcoded.

static constexpr int kSR = 12;

DSFMTGen::DSFMTGen(int mexp, int pos1, int sl1,
                   uint64_t msk1, uint64_t msk2, int L)
    : Generator(128 * ((mexp - 128) / 104 + 2), L),
      mexp_(mexp),
      N_((mexp - 128) / 104 + 1),
      pos1_(pos1), sl1_(sl1),
      msk1_(msk1), msk2_(msk2),
      buf_idx_(2 * ((mexp - 128) / 104 + 1)),
      last_output_(0)
{
    state_ = BitVect(128 * (N_ + 1));
    work_  = BitVect(128 * (N_ + 1));
}

std::string DSFMTGen::name() const { return "dSFMT"; }

std::string DSFMTGen::display_str() const {
    char buf[160];
    snprintf(buf, sizeof(buf),
             " MEXP=%d N=%d POS1=%d SL1=%d "
             "MSK1=%016llx MSK2=%016llx",
             mexp_, N_, pos1_, sl1_,
             (unsigned long long)msk1_,
             (unsigned long long)msk2_);
    return std::string(buf);
}

// ── work buffer ↔ state_ sync ────────────────────────────────────────────

void DSFMTGen::rebuild_work_from_state() {
    work_ = state_.copy();
    bool any = false;
    int slots = 2 * (N_ + 1);
    for (int s = 0; s < slots && !any; s++) {
        if (work_.get_word(s, 64) != 0) any = true;
    }
    if (!any) {
        // Mirror dsfmt period_certification's all-zero fallback by
        // injecting a 1-bit into the lung's u[0].  This breaks the
        // trivial orbit; the matricial test never feeds an all-zero
        // perturbation so this path only triggers from init().
        set_lane_u64(N_, 0, 1ULL);
    }
}

void DSFMTGen::set_raw_state(const BitVect& s) {
    // Direct copy — no all-zero guard (would break F2-linearity for
    // the SIMD-PIS standard basis vectors).
    BitVect ns(128 * (N_ + 1));
    ns.copy_part_from(s, std::min(s.nbits(), 128 * (N_ + 1)));
    state_ = std::move(ns);
    work_ = state_.copy();
}

void DSFMTGen::sync_state_from_work() {
    state_ = work_.copy();
}

// ── recurrence ───────────────────────────────────────────────────────────

void DSFMTGen::do_recursion(int idx) {
    uint64_t a0 = lane_u64(idx, 0);
    uint64_t a1 = lane_u64(idx, 1);
    int b_idx = (idx + pos1_) % N_;
    uint64_t b0 = lane_u64(b_idx, 0);
    uint64_t b1 = lane_u64(b_idx, 1);
    uint64_t L0 = lane_u64(N_, 0);
    uint64_t L1 = lane_u64(N_, 1);

    uint64_t newL0 = (a0 << sl1_) ^ (L1 >> 32) ^ (L1 << 32) ^ b0;
    uint64_t newL1 = (a1 << sl1_) ^ (L0 >> 32) ^ (L0 << 32) ^ b1;
    uint64_t r0 = (newL0 >> kSR) ^ (newL0 & msk1_) ^ a0;
    uint64_t r1 = (newL1 >> kSR) ^ (newL1 & msk2_) ^ a1;

    set_lane_u64(idx, 0, r0);
    set_lane_u64(idx, 1, r1);
    set_lane_u64(N_, 0, newL0);
    set_lane_u64(N_, 1, newL1);
}

void DSFMTGen::generate_all() {
    for (int i = 0; i < N_; i++) {
        do_recursion(i);
    }
}

void DSFMTGen::init(const BitVect& init_bv) {
    state_ = BitVect(128 * (N_ + 1));
    state_.copy_part_from(init_bv,
                          std::min(init_bv.nbits(), 128 * (N_ + 1)));
    rebuild_work_from_state();
    generate_all();
    sync_state_from_work();
    buf_idx_ = 0;
    last_output_ = lane_u64(0, 0) & ((1ULL << 52) - 1);
    buf_idx_ = 1;
}

void DSFMTGen::next() {
    if (buf_idx_ >= 2 * N_) {
        generate_all();
        buf_idx_ = 0;
        sync_state_from_work();
    }
    last_output_ = lane_u64(buf_idx_ / 2, buf_idx_ % 2) & ((1ULL << 52) - 1);
    buf_idx_++;
}

std::unique_ptr<Generator> DSFMTGen::copy() const {
    auto p = std::make_unique<DSFMTGen>(
        mexp_, pos1_, sl1_, msk1_, msk2_, L_);
    p->state_ = state_.copy();
    p->work_ = work_.copy();
    p->buf_idx_ = buf_idx_;
    p->last_output_ = last_output_;
    p->pword_idx_ = pword_idx_;
    return p;
}

BitVect DSFMTGen::get_output() const {
    // Output L_ bits of the most recent lane.  For L=52 (the natural
    // dSFMT output width), the entire low-52-bit value is returned with
    // bit 51 (lane MSB) at BitVect position 0.  For L=64, the high 12
    // bits are zero.
    BitVect out(L_);
    int n = std::min(L_, 52);
    for (int b = 0; b < n; b++) {
        if (last_output_ & (1ULL << (51 - b))) out.set_bit(b, 1);
    }
    return out;
}

// ── SIMD-aware overrides for METHOD_SIMD_NOTPRIMITIVE ────────────────────

void DSFMTGen::simd_advance_one_word() {
    // Mirror MTToolBox dSFMTsearch.hpp::next_state (lines 277-283):
    //   index = (index + 1) % size;
    //   do_recursion(state[index], state[index],
    //                state[(index + pos1) % size], lung);
    pword_idx_ = (pword_idx_ + 1) % N_;
    do_recursion(pword_idx_);
    sync_state_from_work();
}

BitVect DSFMTGen::simd_read_super_word(int sm) const {
    // Mirror MTToolBox dSFMT::generate (dSFMTsearch.hpp:289-332):
    //   case 0: r.u[0] = state[idx].u[0];  r.u[1] = state[idx].u[1]
    //   case 1: r.u[0] = state[prev].u[1]; r.u[1] = state[idx].u[0]
    // Then pack the v MSBs of each lane (bits 51..51-v+1 in the uint64)
    // into the BitVect at MSB-first positions [i*L_, (i+1)*L_).
    //
    // L_ is the per-lane bit-width passed in by the test (52 for dSFMT).
    // For lane i, BitVect position i*L_ + b receives bit (51-b) of u[i].
    BitVect out(128);
    int idx = pword_idx_;
    int prev = (idx + N_ - 1) % N_;
    uint64_t lanes[2];
    if (sm == 0) {
        lanes[0] = lane_u64(idx, 0);
        lanes[1] = lane_u64(idx, 1);
    } else {
        // sm == 1
        lanes[0] = lane_u64(prev, 1);
        lanes[1] = lane_u64(idx, 0);
    }
    for (int i = 0; i < 2; i++) {
        uint64_t v = lanes[i];
        int base = i * L_;
        int n = std::min(L_, 52);
        for (int b = 0; b < n; b++) {
            if (v & (1ULL << (51 - b))) out.set_bit(base + b, 1);
        }
    }
    return out;
}

void DSFMTGen::simd_add_state(const Generator& other) {
    // INDEX-ALIGNED XOR per MTToolBox dSFMT::add (dSFMTsearch.hpp:472-481):
    //   for i in [0, size): state[(i + index) % size] ^= that->state[(i + that->index) % size]
    //   lung ^= that->lung
    //   previous ^= that->previous   (we don't track previous on the
    //                                 generator — SimdLinVec carries it)
    const DSFMTGen* that = dynamic_cast<const DSFMTGen*>(&other);
    if (!that || that->N_ != N_) {
        BitVect s = state_.copy();
        s.xor_with(other.state());
        state_ = std::move(s);
        work_ = state_.copy();
        return;
    }
    int my_idx = pword_idx_;
    int other_idx = that->pword_idx_;
    for (int i = 0; i < N_; i++) {
        int dst = (i + my_idx) % N_;
        int src = (i + other_idx) % N_;
        for (int l = 0; l < 2; l++) {
            uint64_t v = lane_u64(dst, l) ^ that->lane_u64(src, l);
            set_lane_u64(dst, l, v);
        }
    }
    // Lung: index-independent XOR.
    set_lane_u64(N_, 0, lane_u64(N_, 0) ^ that->lane_u64(N_, 0));
    set_lane_u64(N_, 1, lane_u64(N_, 1) ^ that->lane_u64(N_, 1));
    state_ = work_.copy();
}

// ── factory ─────────────────────────────────────────────────────────────

std::unique_ptr<Generator> DSFMTGen::from_params(
    const Params& params, int L)
{
    const int mexp = (int)params.get_int("mexp");
    const int pos1 = (int)params.get_int("pos1");
    const int sl1  = (int)params.get_int("sl1");
    const uint64_t msk1 = (uint64_t)params.get_int("msk1");
    const uint64_t msk2 = (uint64_t)params.get_int("msk2");
    return std::make_unique<DSFMTGen>(mexp, pos1, sl1, msk1, msk2, L);
}

std::vector<ParamSpec> DSFMTGen::param_specs() {
    return {
        {"mexp", "int", true, false, 0, "", "", false},
        {"pos1", "int", true, false, 0, "", "", false},
        {"sl1",  "int", true, false, 0, "", "", false},
        {"msk1", "int", true, false, 0, "", "", false},
        {"msk2", "int", true, false, 0, "", "", false},
    };
}
