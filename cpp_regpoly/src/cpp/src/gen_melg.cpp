#include "gen_melg.h"
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cstring>

MELG::MELG(int w, int N, int M, int r, int sigma1, int sigma2,
           uint64_t a, int L)
    : Generateur(N * w - r, L),
      w_(w), N_(N), M_(M), r_(r),
      sigma1_(sigma1), sigma2_(sigma2), a_(a),
      i_(0), arr_size_(N - 1),
      state_bits_(N * w)
{
    if (w >= 64)
        word_mask_ = ~0ULL;
    else
        word_mask_ = (1ULL << w) - 1;

    if (r == 0) {
        upper_mask_ = word_mask_;
        lower_mask_ = 0;
    } else {
        lower_mask_ = (1ULL << r) - 1;
        upper_mask_ = word_mask_ & ~lower_mask_;
    }

    // State stores N words of w bits: (N-1) array words + 1 word for v.
    // Total = N*w bits.  Of these, k_ = N*w - r are independent.
    state_ = BitVect(state_bits_);
}

std::string MELG::name() const { return "MELG"; }

std::string MELG::display_str() const {
    std::ostringstream oss;
    oss << k_ << "\n"
        << " w=" << w_ << "  N=" << N_ << "  M=" << M_
        << "  r=" << r_
        << "  sigma1=" << sigma1_ << "  sigma2=" << sigma2_
        << "  a=" << std::hex << std::setfill('0') << std::setw(w_ / 4) << a_
        << std::dec;
    return oss.str();
}

// ── Word access ─────────────────────────────────────────────────────────
// Direct storage: w[0] at word 0, w[1] at word 1, ..., w[N-2] at word N-2,
// v at word N-1.  Total N words of w bits.

uint64_t MELG::V(int idx) const {
    return state_.get_word(idx, w_);
}

void MELG::SetV(int idx, uint64_t val) {
    state_.set_word(idx, w_, val & word_mask_);
}

// ── Init ────────────────────────────────────────────────────────────────

void MELG::init(const BitVect& init_bv) {
    state_ = BitVect(state_bits_);
    i_ = 0;

    // Unpack from canonical layout:
    //   bits [0..w-r-1]            → top (w-r) bits of w[0]
    //   bits [w-r..w-r+(N-2)*w-1]  → w[1], w[2], ..., w[N-2]
    //   bits [w-r+(N-2)*w..k-1]    → v
    int pos = 0;

    // w[0]: only top (w-r) bits
    {
        uint64_t word = 0;
        for (int b = 0; b < w_ - r_ && pos < k_; b++, pos++) {
            if (init_bv.get_bit(pos))
                word |= (1ULL << (w_ - 1 - b));
        }
        SetV(0, word);
    }

    // w[1] through w[N-2]: full w-bit words
    for (int j = 1; j < arr_size_; j++) {
        uint64_t word = 0;
        for (int b = 0; b < w_ && pos < k_; b++, pos++) {
            if (init_bv.get_bit(pos))
                word |= (1ULL << (w_ - 1 - b));
        }
        SetV(j, word);
    }

    // v: full w bits
    {
        uint64_t word = 0;
        for (int b = 0; b < w_ && pos < k_; b++, pos++) {
            if (init_bv.get_bit(pos))
                word |= (1ULL << (w_ - 1 - b));
        }
        SetV(arr_size_, word);
    }
}

// ── Rotated state ───────────────────────────────────────────────────────
// Rotate the array part by i_*w_ so the canonical representation is
// independent of the pointer position.  v (the last word) stays in place.

BitVect MELG::rotated_state() const {
    // Rotate so that the output word w[prev] (the word just written)
    // is at position 0.  After next(), i_ has been incremented, so
    // prev = (i_ - 1 + arr_size_) % arr_size_.
    int prev = (i_ + arr_size_ - 1) % arr_size_;
    if (prev == 0)
        return state_.copy();

    int arr_bits = arr_size_ * w_;
    int rotation = (prev * w_) % arr_bits;

    if (rotation == 0)
        return state_.copy();

    // Rotate just the array part; leave v in place.
    BitVect result(state_bits_);

    // Copy v (last w_ bits) unchanged
    for (int b = 0; b < w_; b++) {
        int pos = arr_bits + b;
        if (state_.get_bit(pos))
            result.set_bit(pos, 1);
    }

    // Rotate the array: shift right by 'rotation' bits, wrapping around.
    // New bit i = old bit (i + rotation) % arr_bits
    for (int i = 0; i < arr_bits; i++) {
        int src = (i + rotation) % arr_bits;
        if (state_.get_bit(src))
            result.set_bit(i, 1);
    }

    return result;
}

// ── Recurrence ──────────────────────────────────────────────────────────

void MELG::next() {
    int idx_i  = i_;
    int idx_i1 = (i_ + 1) % arr_size_;
    int idx_m  = (i_ + M_) % arr_size_;

    uint64_t wi  = V(idx_i);
    uint64_t wi1 = V(idx_i1);
    uint64_t wm  = V(idx_m);
    uint64_t v   = V(arr_size_);  // v is stored at index arr_size_

    // x = (w[i] & upper_mask) ^ (w[i+1] & lower_mask)
    uint64_t x = (wi & upper_mask_) ^ (wi1 & lower_mask_);

    // xA = (x >> 1) ^ (LSB(x) ? a : 0)   [twist matrix, same as MT]
    uint64_t xA = (x >> 1);
    if (x & 1ULL)
        xA ^= a_;
    xA &= word_mask_;

    // v_new = xA ^ w[i+M] ^ vB,  where vB = v ^ (v << sigma1)
    uint64_t vB = (v ^ (v << sigma1_)) & word_mask_;
    uint64_t v_new = (xA ^ wm ^ vB) & word_mask_;

    // w[i] = x ^ vC,  where vC = v_new ^ (v_new >> sigma2)
    uint64_t vC = (v_new ^ (v_new >> sigma2_)) & word_mask_;
    uint64_t new_wi = (x ^ vC) & word_mask_;

    SetV(idx_i, new_wi);
    SetV(arr_size_, v_new);  // update v

    i_ = idx_i1;
}

// ── Output ──────────────────────────────────────────────────────────────

BitVect MELG::get_output() const {
    // Return the full rotated state so that:
    //   - the first w bits are the output word w[prev]
    //   - subsequent words are accessible for lag-based tempering
    // The caller (gauss.cpp) will take the first L bits after applying
    // transformations.
    return rotated_state();
}

void MELG::get_transition_state(uint64_t* out_words, int out_nwords) const {
    // Return the state in canonical (logical) order so that
    // init(get_transition_state()) reproduces the same configuration.
    //
    // The canonical layout is the paper's state vector:
    //   (w^{w-r}_{i}, w_{i+1}, ..., w_{i+N-2}, v)
    // This means we read the array starting from logical position i_
    // (wrapping around), then v at the end.

    BitVect tmp(k_);
    int pos = 0;

    // Array words in logical order: i_, (i_+1)%s, ..., (i_+s-1)%s
    for (int j = 0; j < arr_size_; j++) {
        int logical = (i_ + j) % arr_size_;
        uint64_t word = V(logical);
        int bits_to_write = (j == 0) ? (w_ - r_) : w_;
        int bit_start = (j == 0) ? 0 : 0;  // always from MSB

        for (int b = 0; b < bits_to_write && pos < k_; b++, pos++) {
            if ((word >> (w_ - 1 - b)) & 1)
                tmp.set_bit(pos, 1);
        }
    }

    // v
    uint64_t v = V(arr_size_);
    for (int b = 0; b < w_ && pos < k_; b++, pos++) {
        if ((v >> (w_ - 1 - b)) & 1)
            tmp.set_bit(pos, 1);
    }

    int n = std::min(out_nwords, tmp.nwords());
    for (int i = 0; i < n; i++)
        out_words[i] = tmp.data()[i];
    for (int i = n; i < out_nwords; i++)
        out_words[i] = 0;
}

// ── Characteristic polynomial ────────────────────────────────────────────
// The output MSB doesn't observe the full state (v is hidden).
// Use BM on the first bit of get_transition_state instead, which is
// a linear functional of the canonical state and sees all dimensions.

BitVect MELG::char_poly() const {
    int K = k_;
    auto gen = copy();

    BitVect init_bv(K);
    init_bv.set_bit(0, 1);
    gen->init(init_bv);

    int nw = (K + 63) / 64;
    std::vector<uint64_t> state_words(nw, 0);

    std::vector<int> seq(2 * K, 0);
    for (int i = 0; i < 2 * K; i++) {
        gen->next();
        gen->get_transition_state(state_words.data(), nw);
        seq[i] = (state_words[0] >> 63) & 1;
    }

    // Berlekamp-Massey
    std::vector<int> C(K + 1, 0);
    std::vector<int> B(K + 1, 0);
    C[0] = B[0] = 1;
    int Len = 0;
    int x = 1;

    for (int N = 0; N < 2 * K; N++) {
        int d = seq[N];
        for (int j = 1; j <= Len; j++) {
            if (C[j]) d ^= seq[N - j];
        }
        d &= 1;

        if (d == 0) {
            x += 1;
        } else if (2 * Len > N) {
            std::vector<int> temp(K + 1, 0);
            for (int i = 0; i + x <= K; i++) temp[i + x] = B[i];
            for (int i = 0; i <= K; i++) C[i] ^= temp[i];
            x += 1;
        } else {
            std::vector<int> T = C;
            std::vector<int> temp(K + 1, 0);
            for (int i = 0; i + x <= K; i++) temp[i + x] = B[i];
            for (int i = 0; i <= K; i++) C[i] ^= temp[i];
            Len = N + 1 - Len;
            B = T;
            x = 1;
        }
    }

    BitVect bv(K);
    for (int j = 0; j < K; j++)
        bv.set_bit(j, C[K - j]);
    return bv;
}

// ── Copy ────────────────────────────────────────────────────────────────

std::unique_ptr<Generateur> MELG::copy() const {
    auto g = std::make_unique<MELG>(w_, N_, M_, r_, sigma1_, sigma2_, a_, L_);
    g->state_ = state_.copy();
    g->i_ = i_;
    return g;
}
