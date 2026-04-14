#include "generateur.h"
#include <vector>
#include <cstring>
#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/ZZ.h>

// ── Default get_output ──────────────────────────────────────────────────

BitVect Generateur::get_output() const {
    BitVect out(L_);
    out.copy_part_from(state_, L_);
    return out;
}

// ── Default get_transition_state ───────────────────────────────��─────────

void Generateur::get_transition_state(uint64_t* out_words, int out_nwords) const {
    // Copy first k_ bits of state_ into out_words.
    // out_words is MSB-first: out_words[0] holds bits 0..63, etc.
    // state_ already stores bits MSB-first, so we can copy directly.
    BitVect tmp(k_);
    tmp.copy_part_from(state_, k_);
    int n = std::min(out_nwords, tmp.nwords());
    for (int i = 0; i < n; i++)
        out_words[i] = tmp.data()[i];
    for (int i = n; i < out_nwords; i++)
        out_words[i] = 0;
}

// ── Packed Berlekamp-Massey ──────────────────────────────────────────────
// Polynomials C, B and the output sequence are stored as packed uint64_t
// arrays (MSB-first, matching BitVect convention).  This gives ~64×
// speedup over the element-wise int-array version for the inner loop.

namespace {

inline int pk_get(const uint64_t* data, int i) {
    return (data[i >> 6] >> (63 - (i & 63))) & 1;
}

inline void pk_set1(uint64_t* data, int i) {
    data[i >> 6] |= 1ULL << (63 - (i & 63));
}

// Extract 64 contiguous bits starting at bit position pos.
// Caller must ensure data[] has ≥ 1 word of padding at the end.
inline uint64_t pk_extract64(const uint64_t* data, int pos) {
    int w = pos >> 6;
    int b = pos & 63;
    if (b == 0) return data[w];
    return (data[w] << b) | (data[w + 1] >> (64 - b));
}

// C ^= (B right-shifted by 'shift' bits).
// Polynomial semantics: C(z) += z^shift · B(z).
static void pk_shift_xor(uint64_t* C, const uint64_t* B, int shift, int nw) {
    int ws = shift >> 6;
    int bs = shift & 63;
    if (ws >= nw) return;
    if (bs == 0) {
        for (int i = nw - 1; i >= ws; i--)
            C[i] ^= B[i - ws];
    } else {
        int comp = 64 - bs;
        for (int i = nw - 1; i > ws; i--)
            C[i] ^= (B[i - ws] >> bs) | (B[i - ws - 1] << comp);
        C[ws] ^= B[0] >> bs;
    }
}

// Core packed BM.  Runs the generator from init_state for 2K steps.
// Returns linear complexity L.
// If out_poly is non-null, stores the connection polynomial as a BitVect.
static int packed_bm(const Generateur& gen, const BitVect& init_state,
                     int K, BitVect* out_poly)
{
    auto g = gen.copy();
    g->init(init_state);

    int N_total = 2 * K;
    int nw_seq  = (N_total + 63) / 64 + 2;  // +2 padding for pk_extract64
    int nw_poly = (K + 1 + 63) / 64 + 1;    // +1 padding

    // Reversed sequence: rseq bit i stores s[N_total-1-i].
    // This makes the discrepancy window contiguous for word-level AND+popcount.
    std::vector<uint64_t> rseq(nw_seq, 0);

    std::vector<uint64_t> C(nw_poly, 0);
    std::vector<uint64_t> B(nw_poly, 0);
    std::vector<uint64_t> T(nw_poly);  // pre-allocated temp

    C[0] = 1ULL << 63;  // C[0] = 1 (bit 0 = MSB)
    B[0] = 1ULL << 63;  // B[0] = 1

    int Len = 0;
    int x = 1;

    for (int N = 0; N < N_total; N++) {
        g->next();
        int sN = g->get_output().get_bit(0);

        if (sN) pk_set1(rseq.data(), N_total - 1 - N);

        // Discrepancy: d = Σ_{j=0}^{Len} C[j] · s[N-j]
        //              = Σ_{j=0}^{Len} C[j] · rseq[base + j]
        int base = N_total - 1 - N;
        int d = 0;
        int nw_active = (Len + 64) / 64;
        for (int w = 0; w < nw_active; w++) {
            uint64_t cw = C[w];
            if (cw)
                d ^= __builtin_popcountll(cw & pk_extract64(rseq.data(), base + 64 * w));
        }
        d &= 1;

        if (d == 0) {
            x++;
        } else if (2 * Len > N) {
            pk_shift_xor(C.data(), B.data(), x, nw_poly);
            x++;
        } else {
            std::memcpy(T.data(), C.data(), nw_poly * sizeof(uint64_t));
            pk_shift_xor(C.data(), B.data(), x, nw_poly);
            Len = N + 1 - Len;
            std::memcpy(B.data(), T.data(), nw_poly * sizeof(uint64_t));
            x = 1;
        }
    }

    if (out_poly) {
        *out_poly = BitVect(K);
        for (int j = 0; j < K; j++)
            if (pk_get(C.data(), K - j))
                out_poly->set_bit(j, 1);
    }
    return Len;
}

// ── Splitmix64 seed generator ───────────────────────────────────────────

inline uint64_t splitmix64(uint64_t& state) {
    state += 0x9e3779b97f4a7c15ULL;
    uint64_t z = state;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

}  // namespace

// ── char_poly (default: packed Berlekamp-Massey) ────────────────────────

BitVect Generateur::char_poly() const {
    int K = k_;
    BitVect init_bv(K);
    init_bv.set_bit(0, 1);
    BitVect result;
    packed_bm(*this, init_bv, K, &result);
    return result;
}

// ── LCP-based full period test ──────────────────────────────────────────
// For Mersenne prime exponents: irreducible ⟺ primitive.
// Run BM with num_seeds random seeds.  If any gives L < k, the
// characteristic polynomial is reducible (certain).  If all give L = k,
// it is primitive with probability ≥ 1 − 2^{−num_seeds}.

bool is_full_period_lcp(const Generateur& gen, int num_seeds) {
    int K = gen.k();
    for (int seed_idx = 0; seed_idx < num_seeds; seed_idx++) {
        BitVect seed(K);
        uint64_t sm = (uint64_t)(seed_idx + 1) * 0x123456789ABCDEFULL;
        for (int w = 0; w < seed.nwords(); w++)
            seed.data()[w] = splitmix64(sm);
        // Clear tail bits
        int rem = K % 64;
        if (rem != 0)
            seed.data()[seed.nwords() - 1] &= (~0ULL) << (64 - rem);
        // Ensure nonzero
        bool all_zero = true;
        for (int w = 0; w < seed.nwords() && all_zero; w++)
            if (seed.data()[w]) all_zero = false;
        if (all_zero) seed.set_bit(0, 1);

        int L = packed_bm(gen, seed, K, nullptr);
        if (L != K)
            return false;
    }
    return true;
}

// ── Primitivity test ─────────────────────────────────────────────────────

// Known Mersenne prime exponents (2^p - 1 is prime).
// For these, irreducible over GF(2) ↔ primitive.
static bool is_mersenne_prime_exponent(int p) {
    static const int mp[] = {
        2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607,
        1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213,
        19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091,
        756839, 859433, 1257787, 1398269, 2976221, 3021377, 6972593,
        13466917, 20996011, 24036583, 25964951, 30402457, 32582657,
        37156667, 42643801, 43112609, 57885161, 74207281, 77232917,
        82589933, 0
    };
    for (const int* q = mp; *q; q++)
        if (*q == p) return true;
    return false;
}

// Factor 2^k - 1 for small k using trial division.
// Returns the distinct prime factors as NTL ZZ values.
static std::vector<NTL::ZZ> factor_mersenne(int k) {
    NTL::ZZ n = NTL::power(NTL::ZZ(2), k) - 1;
    std::vector<NTL::ZZ> factors;
    NTL::ZZ d(2);
    // Trial division up to sqrt(n) or a reasonable limit
    NTL::ZZ limit = NTL::SqrRoot(n) + 1;
    // Cap trial division for very large numbers
    if (limit > NTL::ZZ(1) << 30)
        limit = NTL::ZZ(1) << 30;
    while (d <= limit && n > 1) {
        if (n % d == 0) {
            factors.push_back(d);
            while (n % d == 0) n /= d;
            limit = NTL::SqrRoot(n) + 1;
            if (limit > NTL::ZZ(1) << 30)
                limit = NTL::ZZ(1) << 30;
        }
        d++;
    }
    if (n > 1) factors.push_back(n);
    return factors;
}

bool Generateur::is_full_period() const {
    int K = k_;
    BitVect cp = char_poly();

    NTL::GF2X f;
    NTL::SetCoeff(f, K);
    for (int j = 0; j < K; j++)
        if (cp.get_bit(j))
            NTL::SetCoeff(f, j);

    // Quick check: constant term must be 1
    if (!IsOne(coeff(f, 0)))
        return false;

    // Check irreducibility
    if (NTL::IterIrredTest(f) == 0)
        return false;

    // For Mersenne prime exponents: irreducible → primitive
    if (is_mersenne_prime_exponent(K))
        return true;

    // General case: check that x has order exactly 2^K - 1 in GF(2)[x]/f(x).
    // For each prime factor p of 2^K - 1, verify x^((2^K-1)/p) ≠ 1 mod f.
    NTL::GF2XModulus F;
    NTL::build(F, f);

    NTL::ZZ order = NTL::power(NTL::ZZ(2), K) - 1;
    auto factors = factor_mersenne(K);

    for (const auto& p : factors) {
        NTL::ZZ exp = order / p;
        NTL::GF2X r;
        NTL::PowerXMod(r, exp, F);
        if (IsOne(r))
            return false;
    }
    return true;
}

// ── Transition matrix ────────────────────────────────────────────────────

std::vector<BitVect> Generateur::transition_matrix() const {
    int K = k_;
    auto gen = copy();

    BitVect bc(K);
    bc.set_bit(0, 1);

    int nw = (K + 63) / 64;

    // Phase 1: build columns (= rows of Aᵀ). Column i is the k-bit
    // state after one step with basis vector e_i.
    std::vector<std::vector<uint64_t>> cols(K, std::vector<uint64_t>(nw, 0));
    for (int i = 0; i < K; i++) {
        gen->init(bc);
        bc.rshift(1);
        gen->next();
        gen->get_transition_state(cols[i].data(), nw);
    }

    // Phase 2: transpose into rows. Row[row] bit i is set iff
    // cols[i] has bit row set.
    std::vector<BitVect> rows(K, BitVect(K));
    for (int i = 0; i < K; i++) {
        for (int row = 0; row < K; row++) {
            // Check bit 'row' in cols[i]
            int w = row / 64;
            int b = 63 - (row % 64);
            if ((cols[i][w] >> b) & 1)
                rows[row].set_bit(i, 1);
        }
    }
    return rows;
}
