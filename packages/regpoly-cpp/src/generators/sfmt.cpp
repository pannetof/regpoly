#include "sfmt.h"
#include <algorithm>
#include <cstdio>
#include <stdexcept>

// ── 128-bit byte-shifts (lshift128 / rshift128 from the SFMTGen spec) ──
// Both operate on a 128-bit quantity represented as four 32-bit lanes
// u[0..3] with u[0] the low 32 bits.

static void lshift128(uint32_t out[4], const uint32_t in[4], int bytes) {
    const int bits = bytes * 8;
    uint64_t tl = ((uint64_t)in[1] << 32) | in[0];
    uint64_t th = ((uint64_t)in[3] << 32) | in[2];
    uint64_t oh, ol;
    if (bits == 0)       { oh = th;                                    ol = tl; }
    else if (bits < 64)  { oh = (th << bits) | (tl >> (64 - bits));    ol = tl << bits; }
    else if (bits == 64) { oh = tl;                                    ol = 0; }
    else if (bits < 128) { oh = tl << (bits - 64);                     ol = 0; }
    else                 { oh = 0;                                     ol = 0; }
    out[0] = (uint32_t)ol;
    out[1] = (uint32_t)(ol >> 32);
    out[2] = (uint32_t)oh;
    out[3] = (uint32_t)(oh >> 32);
}

static void rshift128(uint32_t out[4], const uint32_t in[4], int bytes) {
    const int bits = bytes * 8;
    uint64_t tl = ((uint64_t)in[1] << 32) | in[0];
    uint64_t th = ((uint64_t)in[3] << 32) | in[2];
    uint64_t oh, ol;
    if (bits == 0)       { oh = th;                                    ol = tl; }
    else if (bits < 64)  { ol = (tl >> bits) | (th << (64 - bits));    oh = th >> bits; }
    else if (bits == 64) { ol = th;                                    oh = 0; }
    else if (bits < 128) { ol = th >> (bits - 64);                     oh = 0; }
    else                 { ol = 0;                                     oh = 0; }
    out[0] = (uint32_t)ol;
    out[1] = (uint32_t)(ol >> 32);
    out[2] = (uint32_t)oh;
    out[3] = (uint32_t)(oh >> 32);
}

// ── Per-MEXP parity constants (from SFMTGen-paramsNNN.h) ──────────────
// Resolved once at construction time into {parity0..parity3}.

static void parity_for(int mexp, uint32_t p[4]) {
    p[0] = p[1] = p[2] = p[3] = 0;
    switch (mexp) {
        case 607:    p[0]=0x00000001; p[3]=0x5986f054; break;
        case 1279:   p[0]=0x00000001; p[3]=0x20000000; break;
        case 2281:   p[0]=0x00000001; p[3]=0x41dfa600; break;
        case 4253:   p[0]=0xa8000001; p[1]=0xaf5390a3;
                     p[2]=0xb740b3f8; p[3]=0x6c11486d; break;
        case 11213:  p[0]=0x00000001; p[2]=0xe8148000;
                     p[3]=0xd0c7afa3; break;
        case 19937:  p[0]=0x00000001; p[3]=0x13c9e684; break;
        case 44497:  p[0]=0x00000001; p[2]=0xa3ac4000;
                     p[3]=0xecc1327a; break;
        case 86243:  p[0]=0x00000001; p[3]=0xe9528d85; break;
        case 132049: p[0]=0x00000001; p[2]=0xcb520000;
                     p[3]=0xc7e91c7d; break;
        case 216091: p[0]=0xf8000001; p[1]=0x89e80709;
                     p[2]=0x3bd2b64b; p[3]=0x0c64b1e4; break;
        default:     break;  // unknown MEXP — no parity table
    }
}

SFMTGen::SFMTGen(int mexp, int pos1, int sl1, int sl2, int sr1, int sr2,
           uint32_t msk1, uint32_t msk2, uint32_t msk3, uint32_t msk4,
           int L)
    : Generator(128 * (mexp / 128 + 1), L),
      mexp_(mexp), N_(mexp / 128 + 1),
      pos1_(pos1), sl1_(sl1), sl2_(sl2), sr1_(sr1), sr2_(sr2),
      msk1_(msk1), msk2_(msk2), msk3_(msk3), msk4_(msk4),
      buf_idx_(4 * (mexp / 128 + 1)),
      last_output_(0)
{
    // state_ (inherited) is k = 128·N bits — the true F2-linear state of
    // the SFMTGen recurrence.  The earlier convention used k = MEXP, which
    // is the period exponent but not the dimension of the linear state
    // space; that was wrong for non-full-period analyses.  work_ is kept
    // as a parallel buffer for the lane-oriented recurrence code below
    // and is the same size as state_; the two are kept in sync via
    // rebuild_work_from_state() / sync_state_from_work().
    state_ = BitVect(128 * N_);
    work_  = BitVect(128 * N_);

    uint32_t p[4];
    parity_for(mexp, p);
    parity0_ = p[0]; parity1_ = p[1]; parity2_ = p[2]; parity3_ = p[3];
}

std::string SFMTGen::name() const { return "SFMTGen"; }

std::string SFMTGen::display_str() const {
    char buf[192];
    snprintf(buf, sizeof(buf),
             " MEXP=%d N=%d POS1=%d SL1=%d SL2=%d SR1=%d SR2=%d "
             "MSK=%08x %08x %08x %08x",
             mexp_, N_, pos1_, sl1_, sl2_, sr1_, sr2_,
             msk1_, msk2_, msk3_, msk4_);
    return std::string(buf);
}

// ── work buffer ↔ state_ sync ─────────────────────────────────────────

void SFMTGen::rebuild_work_from_state() {
    // state_ and work_ now have the same size (128·N bits) and the
    // same meaning — work_ is just the lane-addressable view used by
    // the recurrence code.  The sync is a straight copy, which is
    // exactly linear in state_, so the matricial test (which perturbs
    // state_ by canonical basis vectors) sees the true F2-linear
    // transition.  The all-zero guard remains for the only realistic
    // case it triggers — init() called with an all-zero seed; the
    // matricial test never feeds an all-zero perturbation.
    work_ = state_.copy();
    bool any = false;
    for (int i = 0; i < 4 * N_ && !any; i++) {
        if (work_.get_word(i, 32) != 0) any = true;
    }
    if (!any) set_lane(0, 0, 1);
}

void SFMTGen::set_raw_state(const BitVect& s) {
    // Set state_ to s (truncated/padded to 128*N bits), mirror state_
    // into work_ via DIRECT identity copy.  buf_idx and last_output_
    // are intentionally left untouched: PIS expects subsequent next()
    // calls to read the SAME lane positions as before the state XOR,
    // just with the new (XOR'd) state values.
    //
    // Critical: do NOT call rebuild_work_from_state — its all-zero
    // guard inserts a synthetic 1-bit when state collapses to zero,
    // which BREAKS F2-linearity for the SIMD-PIS standard basis vectors
    // (which legitimately need zero gen state with nonzero `next` bits).
    BitVect new_state(128 * N_);
    new_state.copy_part_from(s, std::min(s.nbits(), 128 * N_));
    state_ = std::move(new_state);
    work_ = state_.copy();
}

void SFMTGen::sync_state_from_work() {
    // Identity copy: state_ and work_ are the same size and mean the
    // same thing.
    state_ = work_.copy();
}

void SFMTGen::period_certify() {
    if (!(parity0_ | parity1_ | parity2_ | parity3_)) return;
    const uint32_t parity[4] = {parity0_, parity1_, parity2_, parity3_};
    uint32_t inner = 0;
    for (int i = 0; i < 4; i++) inner ^= lane(0, i) & parity[i];
    for (int s = 16; s > 0; s >>= 1) inner ^= inner >> s;
    if ((inner & 1) == 1) return;   // already in main orbit
    for (int i = 0; i < 4; i++) {
        for (int b = 0; b < 32; b++) {
            uint32_t m = 1U << b;
            if (m & parity[i]) {
                set_lane(0, i, lane(0, i) ^ m);
                return;
            }
        }
    }
}

void SFMTGen::init(const BitVect& init_bv) {
    // Accept any bit-length init; take the first 128·N bits as the
    // true state.  Shorter inits get zero-padded; longer inits get
    // truncated.  rebuild_work_from_state() handles the all-zero
    // guard so the SFMTGen recurrence never starts in the trivial orbit.
    state_ = BitVect(128 * N_);
    state_.copy_part_from(init_bv,
                          std::min(init_bv.nbits(), 128 * N_));
    rebuild_work_from_state();
    // Run one generate-all immediately so get_output() returns a valid
    // first output word without requiring a next() call.  Matches how
    // the reference's first gen_rand32() call triggers gen_rand_all
    // and returns state[0].u[0]; essential for the matricial test,
    // which reads get_output() at step mc=0 before any next().
    generate_all();
    sync_state_from_work();
    buf_idx_ = 0;
    last_output_ = lane(0, 0);
    buf_idx_ = 1;
}

void SFMTGen::generate_all() {
    const uint32_t msk[4] = { msk1_, msk2_, msk3_, msk4_ };
    for (int i = 0; i < N_; i++) {
        uint32_t a[4], b[4], c[4], d[4];
        for (int l = 0; l < 4; l++) {
            a[l] = lane(i,                  l);
            b[l] = lane((i + pos1_) % N_,   l);
            c[l] = lane((i + N_ - 2) % N_,  l);
            d[l] = lane((i + N_ - 1) % N_,  l);
        }
        uint32_t x[4], y[4];
        lshift128(x, a, sl2_);
        rshift128(y, c, sr2_);
        for (int l = 0; l < 4; l++) {
            uint32_t r =
                a[l]
              ^ x[l]
              ^ ((b[l] >> sr1_) & msk[l])
              ^ y[l]
              ^ (d[l] << sl1_);
            set_lane(i, l, r);
        }
    }
}

void SFMTGen::next() {
    if (buf_idx_ >= 4 * N_) {
        generate_all();
        buf_idx_ = 0;
        // Mirror the MEXP low bits of the freshly regenerated buffer
        // into state_ so any external observer (transition_matrix,
        // lattice test) sees the up-to-date true state.
        sync_state_from_work();
    }
    last_output_ = lane(buf_idx_ / 4, buf_idx_ % 4);
    buf_idx_++;
}

std::unique_ptr<Generator> SFMTGen::copy() const {
    auto p = std::make_unique<SFMTGen>(
        mexp_, pos1_, sl1_, sl2_, sr1_, sr2_,
        msk1_, msk2_, msk3_, msk4_, L_);
    p->state_ = state_.copy();
    p->work_ = work_.copy();
    p->buf_idx_ = buf_idx_;
    p->last_output_ = last_output_;
    p->pword_idx_ = pword_idx_;
    return p;
}

// ── SIMD-aware overrides for METHOD_SIMD_NOTPRIMITIVE (Phase 3) ───────

void SFMTGen::simd_advance_one_word() {
    // Increment pword_idx_ modulo N, then update work_[pword_idx_] via
    // do_recursion using the same lane reads as one iteration of
    // generate_all's body (i = pword_idx_).  Mirrors MTToolBox SFMTGen
    // next_state (sfmtsearch.hpp lines 309–317).
    pword_idx_ = (pword_idx_ + 1) % N_;
    int i = pword_idx_;
    const uint32_t msk[4] = { msk1_, msk2_, msk3_, msk4_ };
    uint32_t a[4], b[4], c[4], d[4];
    for (int l = 0; l < 4; l++) {
        a[l] = lane(i,                  l);
        b[l] = lane((i + pos1_) % N_,   l);
        c[l] = lane((i + N_ - 2) % N_,  l);
        d[l] = lane((i + N_ - 1) % N_,  l);
    }
    uint32_t x[4], y[4];
    auto lshift128 = [](uint32_t out[4], const uint32_t in[4], int bytes) {
        const int bits = bytes * 8;
        uint64_t tl = ((uint64_t)in[1] << 32) | in[0];
        uint64_t th = ((uint64_t)in[3] << 32) | in[2];
        uint64_t oh, ol;
        if (bits == 0)       { oh = th; ol = tl; }
        else if (bits < 64)  { oh = (th << bits) | (tl >> (64 - bits)); ol = tl << bits; }
        else if (bits == 64) { oh = tl; ol = 0; }
        else if (bits < 128) { oh = tl << (bits - 64); ol = 0; }
        else                 { oh = 0; ol = 0; }
        out[0] = (uint32_t)ol; out[1] = (uint32_t)(ol >> 32);
        out[2] = (uint32_t)oh; out[3] = (uint32_t)(oh >> 32);
    };
    auto rshift128 = [](uint32_t out[4], const uint32_t in[4], int bytes) {
        const int bits = bytes * 8;
        uint64_t tl = ((uint64_t)in[1] << 32) | in[0];
        uint64_t th = ((uint64_t)in[3] << 32) | in[2];
        uint64_t oh, ol;
        if (bits == 0)       { oh = th; ol = tl; }
        else if (bits < 64)  { ol = (tl >> bits) | (th << (64 - bits)); oh = th >> bits; }
        else if (bits == 64) { ol = th; oh = 0; }
        else if (bits < 128) { ol = th >> (bits - 64); oh = 0; }
        else                 { ol = 0; oh = 0; }
        out[0] = (uint32_t)ol; out[1] = (uint32_t)(ol >> 32);
        out[2] = (uint32_t)oh; out[3] = (uint32_t)(oh >> 32);
    };
    lshift128(x, a, sl2_);
    rshift128(y, c, sr2_);
    for (int l = 0; l < 4; l++) {
        uint32_t r =
            a[l]
          ^ x[l]
          ^ ((b[l] >> sr1_) & msk[l])
          ^ y[l]
          ^ (d[l] << sl1_);
        set_lane(i, l, r);
    }
    // state_ mirrors work_ (set_raw_state convention); keep them in sync.
    sync_state_from_work();
}

BitVect SFMTGen::simd_read_super_word(int sm) const {
    // Mirror MTToolBox sfmt::generate (sfmtsearch.hpp:334-359):
    //   sm=0:  r.u[0..3] = state[idx].u[0..3]                          (no rotation)
    //   sm=1:  r.u[0..2] = state[prev].u[1..3], r.u[3]   = state[idx].u[0]
    //   sm=2:  r.u[0..1] = state[prev].u[2..3], r.u[2..3]= state[idx].u[0..1]
    //   sm=3:  r.u[0]    = state[prev].u[3],    r.u[1..3]= state[idx].u[0..2]
    // Pack lane u[i] MSB-first into BitVect bits [i·L .. (i+1)·L - 1].
    BitVect out(128);
    int idx  = pword_idx_;
    int prev = (idx + N_ - 1) % N_;
    int lc = 128 / L_;
    for (int i = 0; i < lc; i++) {
        int src_word, src_lane;
        if (sm == 0) {
            src_word = idx;
            src_lane = i;
        } else if (i + sm < lc) {
            src_word = prev;
            src_lane = i + sm;
        } else {
            src_word = idx;
            src_lane = i + sm - lc;
        }
        uint32_t v = (uint32_t)lane(src_word, src_lane);
        for (int b = 0; b < L_; b++) {
            if (v & (1u << (L_ - 1 - b)))
                out.set_bit(i * L_ + b, 1);
        }
    }
    return out;
}

void SFMTGen::simd_add_state(const Generator& other) {
    // INDEX-ALIGNED XOR per MTToolBox sfmt::add (sfmtsearch.hpp:549-562):
    //   for i in [0, N): state[(i + my_idx) % N] ^= other.state[(i + other_idx) % N]
    // The "logical position 0" of each generator is at array index = its
    // pword_idx_ (= MTToolBox's `index`).  Aligning by index ensures the
    // SFMTGen recurrence remains coherent after add(), so subsequent
    // simd_advance_one_word calls advance the summed state correctly.
    //
    // pword_idx_ itself is per-vector metadata, not F2-linear state, and
    // is left unchanged by add (matches MTToolBox).
    const SFMTGen* that = dynamic_cast<const SFMTGen*>(&other);
    if (!that || that->N_ != N_) {
        // Defensive fallback: if `other` isn't an SFMTGen (shouldn't happen
        // in PIS), fall back to flat XOR.
        BitVect s = state_.copy();
        s.xor_with(other.state());
        state_ = std::move(s);
        work_ = state_.copy();
        return;
    }
    int my_idx    = pword_idx_;
    int other_idx = that->pword_idx_;
    for (int i = 0; i < N_; i++) {
        int dst_w = (i + my_idx)    % N_;
        int src_w = (i + other_idx) % N_;
        for (int l = 0; l < 4; l++) {
            uint32_t v = (uint32_t)lane(dst_w, l)
                       ^ (uint32_t)that->lane(src_w, l);
            set_lane(dst_w, l, v);
        }
    }
    // Mirror state_ into work_ directly (identity copy, no zero guard
    // — PIS legitimately produces zero states; rebuild_work_from_state's
    // synthetic 1-bit guard would break F2-linearity).
    state_ = work_.copy();
}

BitVect SFMTGen::get_output() const {
    BitVect out(32);
    out.set_word(0, 32, (uint64_t)last_output_);
    return out;
}

std::unique_ptr<Generator> SFMTGen::from_params(
    const Params& params, int L)
{
    const int mexp = (int)params.get_int("mexp");
    const int pos1 = (int)params.get_int("pos1");
    const int sl1  = (int)params.get_int("sl1");
    const int sl2  = (int)params.get_int("sl2");
    const int sr1  = (int)params.get_int("sr1");
    const int sr2  = (int)params.get_int("sr2");
    auto msk = params.get_uint_vec("msk");
    if (msk.size() != 4)
        throw std::invalid_argument(
            "SFMTGen: 'msk' must be a 4-element uint_vec");
    return std::make_unique<SFMTGen>(
        mexp, pos1, sl1, sl2, sr1, sr2,
        (uint32_t)msk[0], (uint32_t)msk[1],
        (uint32_t)msk[2], (uint32_t)msk[3], L);
}

std::vector<ParamSpec> SFMTGen::param_specs() {
    return {
        {"mexp", "int",      true,  false, 0, "",     "", false},
        {"pos1", "int",      true,  false, 0, "",     "", false},
        {"sl1",  "int",      true,  false, 0, "",     "", false},
        {"sl2",  "int",      true,  false, 0, "",     "", false},
        {"sr1",  "int",      true,  false, 0, "",     "", false},
        {"sr2",  "int",      true,  false, 0, "",     "", false},
        {"msk",  "uint_vec", true,  false, 0, "none", "", false},
    };
}
