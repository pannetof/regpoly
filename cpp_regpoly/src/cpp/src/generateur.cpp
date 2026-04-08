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

// ── Berlekamp-Massey (char_poly) ─────────────────────────────────────────

BitVect Generateur::char_poly() const {
    int K = k_;

    // Create a copy of self
    auto gen = copy();

    // Init with canonical vector (bit 0 set)
    BitVect init_bv(K);
    init_bv.set_bit(0, 1);
    gen->init(init_bv);

    // Collect 2*K output bits: MSB of output each time
    std::vector<int> seq(2 * K, 0);
    for (int i = 0; i < 2 * K; i++) {
        gen->next();
        // Get MSB of output (bit 0) — uses get_output() to handle
        // generators with circular buffers or output rotation
        seq[i] = gen->get_output().get_bit(0);
    }

    // Berlekamp-Massey algorithm (SequenceMinimalPolynomial)
    std::vector<int> C(K + 1, 0);
    std::vector<int> B(K + 1, 0);
    C[0] = B[0] = 1;
    int Len = 0;
    int x = 1;

    for (int N = 0; N < 2 * K; N++) {
        int d = seq[N];
        for (int j = 1; j <= Len; j++) {
            if (C[j])
                d ^= seq[N - j];
        }
        d &= 1;

        if (d == 0) {
            x += 1;
        } else if (2 * Len > N) {
            // C = C - x^x * B (in GF(2): XOR)
            std::vector<int> temp(K + 1, 0);
            for (int i = 0; i + x <= K; i++)
                temp[i + x] = B[i];
            for (int i = 0; i <= K; i++)
                C[i] ^= temp[i];
            x += 1;
        } else {
            std::vector<int> T = C;
            std::vector<int> temp(K + 1, 0);
            for (int i = 0; i + x <= K; i++)
                temp[i + x] = B[i];
            for (int i = 0; i <= K; i++)
                C[i] ^= temp[i];
            Len = N + 1 - Len;
            B = T;
            x = 1;
        }
    }

    // Pack into BitVect: public bit j holds C[K - j]
    BitVect bv(K);
    for (int j = 0; j < K; j++)
        bv.set_bit(j, C[K - j]);
    return bv;
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

    // Build NTL polynomial: f(x) = x^K + sum of cp bits
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
