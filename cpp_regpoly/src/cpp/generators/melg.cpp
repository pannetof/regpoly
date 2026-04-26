#include "melg.h"
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cstring>
#include <NTL/GF2X.h>


MELGGen::MELGGen(int w, int N, int M, int r, int sigma1, int sigma2,
           uint64_t a, int L)
    : Generator(N * w - r, L),
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

std::string MELGGen::name() const { return "MELGGen"; }

std::string MELGGen::display_str() const {
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
// For w=64, each word aligns exactly with one uint64_t in state_.data().

uint64_t MELGGen::V(int idx) const {
    if (w_ == 64) return state_.data()[idx];
    return state_.get_word(idx, w_);
}

void MELGGen::SetV(int idx, uint64_t val) {
    if (w_ == 64) { state_.data()[idx] = val; return; }
    state_.set_word(idx, w_, val & word_mask_);
}

// ── Bit-stream helpers ──────────────────────────────────────────────────
// Extract / insert n bits at an arbitrary bit offset in a uint64_t array
// (MSB-first layout, same as BitVect).

namespace {

// Extract n bits (1 ≤ n ≤ 64) starting at bit position 'pos' from data[].
inline uint64_t extract_bits(const uint64_t* data, int pos, int n) {
    int wi = pos / 64;
    int bi = pos % 64;
    int avail = 64 - bi;
    uint64_t mask = (n == 64) ? ~0ULL : ((1ULL << n) - 1);
    if (n <= avail) {
        return (data[wi] >> (avail - n)) & mask;
    } else {
        uint64_t hi = data[wi] & ((1ULL << avail) - 1);
        int rest = n - avail;
        uint64_t lo = data[wi + 1] >> (64 - rest);
        return (hi << rest) | lo;
    }
}

// Insert n bits (1 ≤ n ≤ 64) at bit position 'pos' into data[].
inline void insert_bits(uint64_t* data, int pos, int n, uint64_t val) {
    int wi = pos / 64;
    int bi = pos % 64;
    int avail = 64 - bi;
    uint64_t mask = (n == 64) ? ~0ULL : ((1ULL << n) - 1);
    val &= mask;
    if (n <= avail) {
        int shift = avail - n;
        uint64_t cmask = mask << shift;
        data[wi] = (data[wi] & ~cmask) | (val << shift);
    } else {
        int rest = n - avail;
        uint64_t hi_mask = (1ULL << avail) - 1;
        data[wi] = (data[wi] & ~hi_mask) | (val >> rest);
        uint64_t lo_mask = ((1ULL << rest) - 1) << (64 - rest);
        data[wi + 1] = (data[wi + 1] & ~lo_mask) | ((val & ((1ULL << rest) - 1)) << (64 - rest));
    }
}

}  // namespace

// ── Init ────────────────────────────────────────────────────────────────

void MELGGen::init(const BitVect& init_bv) {
    state_ = BitVect(state_bits_);
    i_ = 0;

    // Unpack from canonical layout (word-level):
    //   bits [0..w-r-1]            → top (w-r) bits of w[0]
    //   bits [w-r..w-r+(N-2)*w-1]  → w[1], w[2], ..., w[N-2]
    //   bits [w-r+(N-2)*w..k-1]    → v
    const uint64_t* src = init_bv.data();
    int pos = 0;

    // w[0]: only top (w-r) bits, placed in MSB positions
    int wr = w_ - r_;
    SetV(0, (wr > 0 ? extract_bits(src, pos, wr) : 0) << r_);
    pos += wr;

    // w[1] through w[N-2]: full w-bit words
    for (int j = 1; j < arr_size_; j++) {
        SetV(j, extract_bits(src, pos, w_));
        pos += w_;
    }

    // v: full w bits
    SetV(arr_size_, extract_bits(src, pos, w_));
}

// ── Rotated state ───────────────────────────────────────────────────────
// Rotate the array part so the output word w[prev] sits at position 0.
// v (the last word) stays in place.

BitVect MELGGen::rotated_state() const {
    int prev = (i_ + arr_size_ - 1) % arr_size_;
    if (prev == 0)
        return state_.copy();

    int s = arr_size_;
    BitVect result(state_bits_);

    if (w_ == 64) {
        // Fast path: words are uint64_t-aligned.
        // Rotate the array: put words [prev..s-1, 0..prev-1] into [0..s-1].
        const uint64_t* src = state_.data();
        uint64_t* dst = result.data();
        int tail = s - prev;
        std::memcpy(dst, src + prev, tail * sizeof(uint64_t));
        std::memcpy(dst + tail, src, prev * sizeof(uint64_t));
        dst[s] = src[s];  // copy v
    } else {
        // General case: bit-level rotation.
        int arr_bits = s * w_;
        int rotation = (prev * w_) % arr_bits;

        // Copy v unchanged
        for (int b = 0; b < w_; b++) {
            int pos = arr_bits + b;
            if (state_.get_bit(pos))
                result.set_bit(pos, 1);
        }
        // Rotate array bits
        for (int i = 0; i < arr_bits; i++) {
            int src = (i + rotation) % arr_bits;
            if (state_.get_bit(src))
                result.set_bit(i, 1);
        }
    }
    return result;
}

// ── Recurrence ──────────────────────────────────────────────────────────

void MELGGen::next() {
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

BitVect MELGGen::get_output() const {
    // Return the full rotated state so that:
    //   - the first w bits are the output word w[prev]
    //   - subsequent words are accessible for lag-based tempering
    // The caller (gauss.cpp) will take the first L bits after applying
    // transformations.
    return rotated_state();
}

void MELGGen::get_transition_state(uint64_t* out_words, int out_nwords) const {
    // Return the state in canonical (logical) order so that
    // init(get_transition_state()) reproduces the same configuration.
    //
    // The canonical layout is the paper's state vector:
    //   (w^{w-r}_{i}, w_{i+1}, ..., w_{i+N-2}, v)
    // This means we read the array starting from logical position i_
    // (wrapping around), then v at the end.

    BitVect tmp(k_);
    uint64_t* dst = tmp.data();
    int pos = 0;

    // Array words in logical order: i_, (i_+1)%s, ..., (i_+s-1)%s
    for (int j = 0; j < arr_size_; j++) {
        int logical = (i_ + j) % arr_size_;
        uint64_t word = V(logical);
        int bits_to_write = (j == 0) ? (w_ - r_) : w_;
        // Extract top 'bits_to_write' bits of the word
        uint64_t val = word >> (w_ - bits_to_write);
        insert_bits(dst, pos, bits_to_write, val);
        pos += bits_to_write;
    }

    // v
    insert_bits(dst, pos, w_, V(arr_size_));

    int n = std::min(out_nwords, tmp.nwords());
    for (int i = 0; i < n; i++)
        out_words[i] = tmp.data()[i];
    for (int i = n; i < out_nwords; i++)
        out_words[i] = 0;
}

// ── Algebraic characteristic polynomial ─────────────────────────────────
// Uses the 2w×2w block characteristic matrix R(t):
//
//   R(t) = [ R11  R12 ]    R11 = t^s I + (I+B2 A)(U+tL) + t^M B2
//          [ R21  R22 ]    R12 = B2 B1
//                          R21 = A(U+tL) + t^M I
//                          R22 = (t+1)I + S_L(sigma1)
//
// where s = N-1, A = twist matrix, B1 = I+S_L(sigma1), B2 = I+S_R(sigma2),
// U = upper mask, L = lower mask.
//
// The characteristic polynomial is phi(t) = det(R(t)) / t^r.

namespace {

// Binary w×w matrix: row i is a uint64_t with bit j at position (w-1-j).
using BinRow = uint64_t;

inline int bin_bit(BinRow row, int j, int w) {
    return (row >> (w - 1 - j)) & 1;
}

inline BinRow bin_set(int j, int w) {
    return 1ULL << (w - 1 - j);
}

// Multiply two binary w×w matrices (w <= 64).
std::vector<BinRow> bin_mat_mul(const std::vector<BinRow>& A,
                                const std::vector<BinRow>& B, int w) {
    // Transpose B for efficient dot products
    std::vector<BinRow> BT(w, 0);
    for (int i = 0; i < w; i++)
        for (int j = 0; j < w; j++)
            if (bin_bit(B[j], i, w))
                BT[i] |= bin_set(j, w);

    std::vector<BinRow> C(w, 0);
    for (int i = 0; i < w; i++)
        for (int j = 0; j < w; j++)
            if (__builtin_popcountll(A[i] & BT[j]) & 1)
                C[i] |= bin_set(j, w);
    return C;
}

// Bareiss determinant of an n×n matrix of GF2X polynomials.
NTL::GF2X bareiss_det_gf2x(std::vector<std::vector<NTL::GF2X>>& M, int n) {
    NTL::GF2X prev;
    NTL::set(prev);  // prev = 1

    for (int k = 0; k < n; k++) {
        // Find pivot
        int pivot = -1;
        for (int i = k; i < n; i++) {
            if (!IsZero(M[i][k])) { pivot = i; break; }
        }
        if (pivot == -1)
            return NTL::GF2X();  // det = 0

        if (pivot != k)
            std::swap(M[k], M[pivot]);  // sign flip irrelevant in GF(2)

        for (int i = k + 1; i < n; i++) {
            for (int j = k + 1; j < n; j++) {
                // Bareiss: M[i][j] = (M[k][k]*M[i][j] + M[i][k]*M[k][j]) / prev
                NTL::GF2X num;
                NTL::mul(num, M[k][k], M[i][j]);
                NTL::GF2X tmp;
                NTL::mul(tmp, M[i][k], M[k][j]);
                num += tmp;
                // Exact division (guaranteed by Bareiss property)
                NTL::div(M[i][j], num, prev);
            }
            NTL::clear(M[i][k]);
        }
        prev = M[k][k];
    }
    return M[n - 1][n - 1];
}

}  // anonymous namespace

BitVect MELGGen::char_poly() const {
    int w = w_;
    int s = arr_size_;       // N-1
    int K = k_;              // N*w - r
    int r = r_;
    int M = M_;

    // For small k, Berlekamp-Massey is faster (O(k²) with tiny constants).
    // The algebraic Bareiss approach wins for large k (crossover ≈ 20000).
    if (K < 20000)
        return Generator::char_poly();

    // ── Build w×w binary matrices ──────────────────────────────────────
    // Bit convention: bit 0 = MSB, bit w-1 = LSB.

    // Identity
    std::vector<BinRow> I_mat(w);
    for (int i = 0; i < w; i++)
        I_mat[i] = bin_set(i, w);

    // Upper mask U: 1 for bits 0..w-r-1
    std::vector<BinRow> U_mat(w, 0);
    for (int i = 0; i < w - r; i++)
        U_mat[i] = bin_set(i, w);

    // Lower mask L: 1 for bits w-r..w-1
    std::vector<BinRow> L_mat(w, 0);
    for (int i = w - r; i < w; i++)
        L_mat[i] = bin_set(i, w);

    // Twist matrix A: A[i][j] = delta(j,i-1) + a_bit[i] * delta(j,w-1)
    std::vector<BinRow> A_mat(w, 0);
    for (int i = 1; i < w; i++)
        A_mat[i] |= bin_set(i - 1, w);    // right shift
    for (int i = 0; i < w; i++) {
        int a_bit = (a_ >> (w - 1 - i)) & 1;
        if (a_bit)
            A_mat[i] |= bin_set(w - 1, w);  // conditional XOR column
    }

    // B1 = I + S_L(sigma1)
    std::vector<BinRow> B1_mat(I_mat);
    for (int i = 0; i < w; i++) {
        int j = i + sigma1_;
        if (j < w) B1_mat[i] ^= bin_set(j, w);
    }

    // B2 = I + S_R(sigma2)
    std::vector<BinRow> B2_mat(I_mat);
    for (int i = 0; i < w; i++) {
        int j = i - sigma2_;
        if (j >= 0) B2_mat[i] ^= bin_set(j, w);
    }

    // S_L(sigma1) shift-only matrix
    std::vector<BinRow> SL_mat(w, 0);
    for (int i = 0; i < w; i++) {
        int j = i + sigma1_;
        if (j < w) SL_mat[i] = bin_set(j, w);
    }

    // ── Compute binary matrix products ────────────────────────────────
    auto B2A         = bin_mat_mul(B2_mat, A_mat, w);
    std::vector<BinRow> I_plus_B2A(w);
    for (int i = 0; i < w; i++)
        I_plus_B2A[i] = I_mat[i] ^ B2A[i];
    auto C1U  = bin_mat_mul(I_plus_B2A, U_mat, w);  // (I+B2A)·U
    auto C1L  = bin_mat_mul(I_plus_B2A, L_mat, w);  // (I+B2A)·L
    auto AU   = bin_mat_mul(A_mat, U_mat, w);
    auto AL   = bin_mat_mul(A_mat, L_mat, w);
    auto B2B1 = bin_mat_mul(B2_mat, B1_mat, w);

    // ── Schur complement reduction to w×w ─────────────────────────────
    // Instead of the full 2w×2w Bareiss, use:
    //   det(R) = det(R22) · det(S)
    // where S = R11 - R12 · R22^{-1} · R21 (Schur complement).
    //
    // R22 = (t+1)I + S_L(sigma1), det(R22) = (t+1)^w (S_L nilpotent).
    // R22^{-1} = sum_{j=0}^{q} (t+1)^{-(j+1)} S_L^j, q = floor((w-1)/sigma1).
    //
    // Clear denominators: Q = (t+1)^{q+1} · R22^{-1} is polynomial.
    // S' = (t+1)^{q+1} · S (polynomial w×w matrix).
    // phi(t) = det(R)/t^r = det(S') / ((t+1)^{wq} · t^r)

    int q = (w - 1) / sigma1_;  // max power of S_L

    // Build (t+1)^j for j = 0..q+1
    std::vector<NTL::GF2X> tp1_pow(q + 2);
    NTL::set(tp1_pow[0]);  // (t+1)^0 = 1
    NTL::GF2X t_plus_1;
    NTL::SetCoeff(t_plus_1, 0);
    NTL::SetCoeff(t_plus_1, 1);
    for (int j = 1; j <= q + 1; j++)
        NTL::mul(tp1_pow[j], tp1_pow[j - 1], t_plus_1);

    // Build Q = sum_{j=0}^{q} (t+1)^{q-j} · S_L^j as a w×w polynomial matrix.
    // Q[i][k] = (t+1)^{q-h} where h = (k-i)/sigma1, if (k-i) >= 0, divisible
    // by sigma1, and h <= q.  Sparse: each row has at most q+1 entries.
    std::vector<std::vector<NTL::GF2X>> Q(w, std::vector<NTL::GF2X>(w));
    for (int i = 0; i < w; i++) {
        for (int h = 0; h <= q; h++) {
            int k = i + h * sigma1_;
            if (k < w)
                Q[i][k] = tp1_pow[q - h];
        }
    }

    // Build R21 as w×w polynomial matrix.
    // R21[i][j] = AU[i][j] + t·AL[i][j] + t^M·δ(i,j)
    NTL::GF2X t_poly; NTL::SetCoeff(t_poly, 1);     // t
    NTL::GF2X t_M;    NTL::SetCoeff(t_M, M);        // t^M
    NTL::GF2X t_1;    NTL::SetCoeff(t_1, 0);        // 1

    std::vector<std::vector<NTL::GF2X>> R21(w, std::vector<NTL::GF2X>(w));
    for (int i = 0; i < w; i++) {
        for (int j = 0; j < w; j++) {
            NTL::GF2X& r = R21[i][j];
            if (bin_bit(AU[i], j, w)) r += t_1;
            if (bin_bit(AL[i], j, w)) r += t_poly;
            if (i == j) r += t_M;
        }
    }

    // Compute P = R12 · Q · R21 (w×w polynomial matrix)
    // R12 = B2B1 (constant binary matrix)
    // First: tmp = Q · R21
    std::vector<std::vector<NTL::GF2X>> QR21(w, std::vector<NTL::GF2X>(w));
    for (int i = 0; i < w; i++) {
        for (int j = 0; j < w; j++) {
            NTL::GF2X& acc = QR21[i][j];
            // Q[i][k] is nonzero only for k = i + h*sigma1 (h=0..q)
            for (int h = 0; h <= q; h++) {
                int k = i + h * sigma1_;
                if (k >= w) break;
                if (!IsZero(Q[i][k]) && !IsZero(R21[k][j])) {
                    NTL::GF2X tmp;
                    NTL::mul(tmp, Q[i][k], R21[k][j]);
                    acc += tmp;
                }
            }
        }
    }

    // P = B2B1 · QR21 (binary matrix × polynomial matrix)
    std::vector<std::vector<NTL::GF2X>> P(w, std::vector<NTL::GF2X>(w));
    for (int i = 0; i < w; i++) {
        for (int j = 0; j < w; j++) {
            NTL::GF2X& acc = P[i][j];
            for (int k = 0; k < w; k++) {
                if (bin_bit(B2B1[i], k, w) && !IsZero(QR21[k][j]))
                    acc += QR21[k][j];
            }
        }
    }

    // Build S' = (t+1)^{q+1} · R11 - P (= (t+1)^{q+1}·R11 + P in GF(2))
    // R11[i][j] = t^s·δ(i,j) + C1U[i][j] + t·C1L[i][j] + t^M·B2[i][j]
    NTL::GF2X t_s; NTL::SetCoeff(t_s, s);

    std::vector<std::vector<NTL::GF2X>> S(w, std::vector<NTL::GF2X>(w));
    for (int i = 0; i < w; i++) {
        for (int j = 0; j < w; j++) {
            // Build R11[i][j]
            NTL::GF2X r11;
            if (i == j) r11 += t_s;
            if (bin_bit(C1U[i], j, w)) r11 += t_1;
            if (bin_bit(C1L[i], j, w)) r11 += t_poly;
            if (bin_bit(B2_mat[i], j, w)) r11 += t_M;

            // S'[i][j] = (t+1)^{q+1} · R11[i][j] + P[i][j]
            NTL::GF2X scaled;
            NTL::mul(scaled, tp1_pow[q + 1], r11);
            S[i][j] = scaled + P[i][j];  // + = XOR in GF(2)
        }
    }

    // ── Compute determinant of w×w matrix S' via Bareiss ─────────────
    NTL::GF2X det = bareiss_det_gf2x(S, w);

    // ── Divide by (t+1)^{wq} · t^r ─────────────────────────────────
    // First divide by t^r
    if (r > 0) {
        NTL::GF2X t_r;
        NTL::SetCoeff(t_r, r);
        NTL::div(det, det, t_r);
    }
    // Then divide by (t+1)^{wq}
    if (q > 0) {
        NTL::GF2X tp1_wq;
        NTL::set(tp1_wq);
        for (int j = 0; j < w * q; j++)
            NTL::mul(tp1_wq, tp1_wq, t_plus_1);
        NTL::div(det, det, tp1_wq);
    }

    // ── Convert to BitVect ──────────────────────────────────────────
    // det is now φ(t) = t^K + c_{K-1}t^{K-1} + ... + c_0
    // BitVect stores: bit j = coefficient of t^j (no leading term)
    BitVect bv(K);
    for (int j = 0; j < K; j++)
        bv.set_bit(j, IsOne(coeff(det, j)) ? 1 : 0);
    return bv;
}

// ── Copy ────────────────────────────────────────────────────────────────

std::unique_ptr<Generator> MELGGen::copy() const {
    auto g = std::make_unique<MELGGen>(w_, N_, M_, r_, sigma1_, sigma2_, a_, L_);
    g->state_ = state_.copy();
    g->i_ = i_;
    return g;
}

// ── Factory methods ────────────────────────────────────────────────────

std::unique_ptr<Generator> MELGGen::from_params(const Params& params, int L) {
    int w = (int)params.get_int("w", 64);
    int N = (int)params.get_int("N");
    int M = (int)params.get_int("M");
    int r = (int)params.get_int("r");
    int sigma1 = (int)params.get_int("sigma1");
    int sigma2 = (int)params.get_int("sigma2");
    uint64_t a = (uint64_t)params.get_int("a");
    return std::make_unique<MELGGen>(w, N, M, r, sigma1, sigma2, a, L);
}

std::vector<ParamSpec> MELGGen::param_specs() {
    return {
        {"w",      "int", true,  true,  64, "",        "", false},
        {"N",      "int", true,  false, 0,  "",        "", false},
        {"r",      "int", true,  false, 0,  "",        "", false},
        {"M",      "int", false, false, 0,  "range",   "1,N-2", false},
        {"sigma1", "int", false, false, 0,  "range",   "1,w-1", false},
        {"sigma2", "int", false, false, 0,  "range",   "1,w-1", false},
        {"a",      "int", false, false, 0,  "bitmask", "w", false},
    };
}
