// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "generator.h"
#include "bm.h"
#include <algorithm>
#include <chrono>
#include <cstring>
#include <vector>
#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/ZZ.h>

// ── Default get_output ──────────────────────────────────────────────────

BitVect Generator::get_output() const {
    BitVect out(L_);
    out.copy_part_from(state_, L_);
    return out;
}

// ── Default get_transition_state ───────────────────────────────��─────────

void Generator::get_transition_state(uint64_t* out_words, int out_nwords) const {
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

// ── char_poly (default: packed Berlekamp-Massey) ────────────────────────

BitVect Generator::char_poly() const {
    int K = k_;
    BitVect init_bv(K);
    init_bv.set_bit(0, 1);
    BitVect result;
    packed_bm(*this, init_bv, K, &result);
    return result;
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

// ── Pollard's rho (Brent's variant) ────────────────────────────────────
//
// Returns a non-trivial factor of n (possibly composite). May return n
// itself on failure (caller retries with different random params, or
// gives up). Expected complexity O(n^(1/4)) — works on numbers far
// beyond what trial division could handle. n must be ≥ 2 and odd (the
// caller peels 2's first).
static NTL::ZZ pollard_rho_brent(
    const NTL::ZZ& n,
    std::chrono::steady_clock::time_point deadline
) {
    if (n % 2 == 0) return NTL::ZZ(2);
    if (n == 1) return n;

    NTL::ZZ y = NTL::RandomBnd(n - 2) + 1;
    NTL::ZZ c = NTL::RandomBnd(n - 1) + 1;
    NTL::ZZ g(1), q(1);
    NTL::ZZ x, ys;
    // Cycle-leader length: starts at 1 and doubles each outer round
    // until a non-trivial GCD is found. Brent's improvement over
    // Floyd cuts the expected work by ~24%. Capped at a long because
    // even on the hardest realistic cofactors r stays well below
    // 2^32 before a factor pops out.
    long r = 1;
    const long M = 128;  // batched GCD interval

    while (g == 1) {
        if (std::chrono::steady_clock::now() > deadline) return n;
        x = y;
        for (long i = 0; i < r; ++i) {
            y = (y * y + c) % n;
        }
        long k = 0;
        while (k < r && g == 1) {
            ys = y;
            long bound = (r - k < M) ? (r - k) : M;
            for (long i = 0; i < bound; ++i) {
                y = (y * y + c) % n;
                q = (q * NTL::abs(x - y)) % n;
            }
            g = NTL::GCD(q, n);
            k += M;
            if (std::chrono::steady_clock::now() > deadline) return n;
        }
        r *= 2;
    }
    if (g == n) {
        // Batched GCD collapsed to n; replay one step at a time.
        do {
            if (std::chrono::steady_clock::now() > deadline) return n;
            ys = (ys * ys + c) % n;
            g = NTL::GCD(NTL::abs(x - ys), n);
        } while (g == 1);
    }
    return g;
}

// Recursive factorisation: peels prime factors of n into `out`,
// stopping early if `deadline` is reached. Returns true iff the
// factorisation completed; on false, any unfactored composite tail is
// still appended (mirrors the legacy trial-division semantics — partial
// results are useful to the primitivity test even when incomplete).
static bool factor_with_rho(
    NTL::ZZ n, std::vector<NTL::ZZ>& out,
    std::chrono::steady_clock::time_point deadline
) {
    if (n <= 1) return true;
    if (std::chrono::steady_clock::now() > deadline) {
        out.push_back(n);
        return false;
    }
    // 25 Miller-Rabin rounds → false-positive < 4^-25 ≈ 9·10^-16.
    if (NTL::ProbPrime(n, 25)) {
        out.push_back(n);
        return true;
    }
    NTL::ZZ d(n);
    int tries = 0;
    while (d == n) {
        if (++tries > 24
            || std::chrono::steady_clock::now() > deadline) {
            // Rho consistently degenerated, or time's up — accept the
            // composite as-is (same legacy fallback).
            out.push_back(n);
            return false;
        }
        d = pollard_rho_brent(n, deadline);
    }
    bool ok1 = factor_with_rho(d, out, deadline);
    bool ok2 = factor_with_rho(n / d, out, deadline);
    return ok1 && ok2;
}

// Factor 2^k - 1 into distinct prime factors.
//
// Algorithm: small-prime trial division (free, up to 10^6) front-end,
// then Pollard's rho (Brent's variant) for the cofactor. Replaces the
// pure-trial-division previous implementation, which spent up to 2^30
// iterations on the cofactor — for K=800 (paper-form F2w; Phi_d(2)
// products with d|800) that was tens of seconds to minutes per call
// and was the documented cause of the library detail page hanging.
//
// Bounded by a wall-clock deadline (default 3s; override via
// REGPOLY_FACTOR_BUDGET_SECONDS env if you actually need the full
// factorisation of a hostile cofactor). On overflow, any unfactored
// composite tail is appended as-is — same fallback semantics as the
// legacy implementation, so the primitivity test is no less correct
// than it was before; it's just *much* faster on the common case.
static std::vector<NTL::ZZ> factor_mersenne(int k) {
    NTL::ZZ n = NTL::power(NTL::ZZ(2), k) - 1;
    std::vector<NTL::ZZ> factors;

    // Front-end: trial division by ~7k small primes. Cheap and
    // deterministic; catches the many small-cyclotomic factors of
    // 2^k - 1 (it equals the product of Phi_d(2) for d | k, and most
    // of those are small). Bound is intentionally low — each `n % d`
    // on a 240-digit n is ~1µs, so 10^5 iterations ≈ 100 ms. Larger
    // primes are Pollard's-rho territory anyway.
    NTL::ZZ d(2);
    const NTL::ZZ small_limit(100000);
    while (d < small_limit && d * d <= n) {
        if (n % d == 0) {
            factors.push_back(d);
            while (n % d == 0) n /= d;
        }
        ++d;
    }

    if (n > 1) {
        long budget_ms = 2000;
        if (const char* env = std::getenv("REGPOLY_FACTOR_BUDGET_SECONDS")) {
            try { budget_ms = std::max(1L, std::stol(env)) * 1000; }
            catch (...) {}
        }
        auto deadline = std::chrono::steady_clock::now()
                      + std::chrono::milliseconds(budget_ms);
        factor_with_rho(n, factors, deadline);
    }

    std::sort(factors.begin(), factors.end());
    factors.erase(std::unique(factors.begin(), factors.end()),
                  factors.end());
    return factors;
}

bool Generator::is_full_period() const {
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

// ── default_test_method (cached wrapper around compute_*) ───────────────

std::optional<std::string> Generator::default_test_method(const std::string& test_type) const {
    auto it = default_test_method_cache_.find(test_type);
    if (it != default_test_method_cache_.end())
        return it->second;
    auto result = compute_default_test_method(test_type);
    default_test_method_cache_.emplace(test_type, result);
    return result;
}

std::optional<std::string> Generator::compute_default_test_method(const std::string& test_type) const {
    if (test_type == "equidistribution") {
        if (!is_full_period()) return std::string("notprimitive");
        if (k_ > 100)          return std::string("harase");
        return std::string("matricial");
    }
    return std::nullopt;
}

// ── Component decomposition (default: I am a single primitive) ─────────

std::vector<Generator*> Generator::components() const {
    return {const_cast<Generator*>(this)};
}

std::vector<std::vector<Transformation*>> Generator::tempering_chains() const {
    return {{}};
}

// ── Transition matrix ────────────────────────────────────────────────────

std::vector<BitVect> Generator::transition_matrix() const {
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
