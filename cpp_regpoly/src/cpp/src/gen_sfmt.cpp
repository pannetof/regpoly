#include "gen_sfmt.h"
#include <algorithm>
#include <cstdio>
#include <stdexcept>

// ── 128-bit byte-shifts (lshift128 / rshift128 from the SFMT spec) ──
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

// ── Per-MEXP parity constants (from SFMT-paramsNNN.h) ──────────────
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

SFMT::SFMT(int mexp, int pos1, int sl1, int sl2, int sr1, int sr2,
           uint32_t msk1, uint32_t msk2, uint32_t msk3, uint32_t msk4,
           int L)
    : Generateur(mexp, L),
      mexp_(mexp), N_(mexp / 128 + 1),
      pos1_(pos1), sl1_(sl1), sl2_(sl2), sr1_(sr1), sr2_(sr2),
      msk1_(msk1), msk2_(msk2), msk3_(msk3), msk4_(msk4),
      buf_idx_(4 * (mexp / 128 + 1)),
      last_output_(0)
{
    // state_ (inherited) is k = MEXP bits — the true state.
    state_ = BitVect(mexp);
    // work_ is 128·N bits of private scratch for the recurrence.
    work_ = BitVect(128 * N_);

    uint32_t p[4];
    parity_for(mexp, p);
    parity0_ = p[0]; parity1_ = p[1]; parity2_ = p[2]; parity3_ = p[3];
}

std::string SFMT::name() const { return "SFMT"; }

std::string SFMT::display_str() const {
    char buf[192];
    snprintf(buf, sizeof(buf),
             " MEXP=%d N=%d POS1=%d SL1=%d SL2=%d SR1=%d SR2=%d "
             "MSK=%08x %08x %08x %08x",
             mexp_, N_, pos1_, sl1_, sl2_, sr1_, sr2_,
             msk1_, msk2_, msk3_, msk4_);
    return std::string(buf);
}

// ── work buffer ↔ state_ sync ─────────────────────────────────────────

void SFMT::rebuild_work_from_state() {
    // Populate work_ from state_ by a straightforward linear
    // embedding: state_ bits land in the low MEXP bits of work_, the
    // remaining 128·N − MEXP bits are set to 0.  This must stay
    // strictly linear (no conditional XOR) because regpoly's
    // matricial test perturbs state_ bit-by-bit and requires a
    // linear transition function.  Period certification is therefore
    // NOT applied here; the caller is responsible for supplying a
    // state that already satisfies the parity constraint, which is
    // the standard convention for F2-linear generators in this
    // framework.
    work_ = BitVect(128 * N_);
    work_.copy_part_from(state_, mexp_);

    // Avoid the all-zero trap — a no-op rather than a conditional
    // XOR, so still linear on the vast majority of non-zero states.
    bool any = false;
    for (int i = 0; i < 4 * N_ && !any; i++) {
        if (work_.get_word(i, 32) != 0) any = true;
    }
    if (!any) set_lane(0, 0, 1);
}

void SFMT::sync_state_from_work() {
    // state_ takes the first MEXP bits of work_.
    BitVect new_state(mexp_);
    new_state.copy_part_from(work_, mexp_);
    state_ = new_state;
}

void SFMT::period_certify() {
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

void SFMT::init(const BitVect& init_bv) {
    // Accept any bit-length init; take the first MEXP bits as the
    // true state and reconstruct the private work buffer from there.
    state_ = BitVect(mexp_);
    state_.copy_part_from(init_bv,
                          std::min(init_bv.nbits(), mexp_));
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

void SFMT::generate_all() {
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

void SFMT::next() {
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

std::unique_ptr<Generateur> SFMT::copy() const {
    auto p = std::make_unique<SFMT>(
        mexp_, pos1_, sl1_, sl2_, sr1_, sr2_,
        msk1_, msk2_, msk3_, msk4_, L_);
    p->state_ = state_.copy();
    p->work_ = work_.copy();
    p->buf_idx_ = buf_idx_;
    p->last_output_ = last_output_;
    return p;
}

BitVect SFMT::get_output() const {
    BitVect out(32);
    out.set_word(0, 32, (uint64_t)last_output_);
    return out;
}

std::unique_ptr<Generateur> SFMT::from_params(
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
            "SFMT: 'msk' must be a 4-element uint_vec");
    return std::make_unique<SFMT>(
        mexp, pos1, sl1, sl2, sr1, sr2,
        (uint32_t)msk[0], (uint32_t)msk[1],
        (uint32_t)msk[2], (uint32_t)msk[3], L);
}

std::vector<ParamSpec> SFMT::param_specs() {
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
