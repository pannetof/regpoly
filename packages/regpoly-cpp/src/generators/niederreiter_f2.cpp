// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.

#include "niederreiter_f2.h"
#include "primitivity.h"
#include "bitvect.h"

#include <cstdint>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace regpoly::core {

namespace {

// Polynomial multiplication in F_2: (a * b) over GF(2)[x], packed as
// uint64_t (bit i = coefficient of x^i). Assumes the product fits in
// 64 bits (caller is responsible for not overflowing).
uint64_t poly_mul_gf2(uint64_t a, uint64_t b) {
    uint64_t result = 0;
    while (b) {
        if (b & 1ULL) result ^= a;
        a <<= 1;
        b >>= 1;
    }
    return result;
}

// Compute the Laurent-series coefficients a_n of 1/q(x) over F_2((x⁻¹))
// for n in [0, max_n]. By convention a_n = 0 for n < d = deg(q),
// a_d = 1, and a_{n'} = sum_{j=0..d-1} q_j · a_{n' - d + j} for n' ≥ d+1
// (this is the recurrence derived from q(x) · A(x⁻¹) = 1).
//
// `q` is packed as uint64_t with bit i = coefficient of x^i. The degree
// d is given explicitly so the caller controls indexing (q may have
// implicit leading bits in higher positions).
std::vector<uint8_t> laurent_inverse(uint64_t q, int d, int max_n) {
    std::vector<uint8_t> a(max_n + 1, 0);
    if (d <= max_n) a[d] = 1;
    for (int np = d + 1; np <= max_n; ++np) {
        uint8_t acc = 0;
        for (int j = 0; j < d; ++j) {
            if ((q >> j) & 1ULL) {
                acc ^= a[np - d + j];
            }
        }
        a[np] = acc;
    }
    return a;
}

}  // namespace

NiederreiterF2Gen::NiederreiterF2Gen(int m, int s_max)
    : DigitalNet(m, s_max) {
    // Polynomial arithmetic uses uint64_t. Worst case degree we touch
    // is q = p_j^r with d = r*e_j; we use n up to k+m in the Laurent
    // recurrence, which is ≤ e_max + m - 1. Require ≤ 63 to be safe.
    // For typical s_max ≤ 32, e_max ≤ 8 — so m ≤ ~55 fits. Pin the
    // documented v1 cap at 30 to leave generous headroom.
    if (m > 30) {
        throw std::invalid_argument(
            "NiederreiterF2Gen: m > 30 not supported in v1 (got "
            + std::to_string(m) + "). Polynomial arithmetic uses uint64_t.");
    }
    enumerate_irreducibles(s_max);
    build_matrices(m, s_max);
}

void NiederreiterF2Gen::enumerate_irreducibles(int s_max) {
    polys_.reserve(s_max);
    degs_.reserve(s_max);
    int count = 0;
    for (int d = 1; count < s_max; ++d) {
        for (uint64_t lo = 0; lo < (1ULL << d); ++lo) {
            // is_irreducible_gf2w_modM expects "z^w + lower_w_bits"
            // with lo packed LSB-first. It returns false for z itself
            // (constant-term-zero degree-1) by NTL/codebase convention,
            // but every monic degree-1 polynomial over F_2 is irreducible
            // mathematically, so accept all degree-1 candidates manually.
            bool irr = (d == 1) ? true : is_irreducible_gf2w_modM(lo, d);
            if (irr) {
                uint64_t full = (1ULL << d) | lo;
                polys_.push_back(full);
                degs_.push_back(d);
                ++count;
                if (count >= s_max) break;
            }
        }
    }
}

void NiederreiterF2Gen::build_matrices(int m, int s_max) {
    auto mats = std::make_shared<std::vector<std::vector<BitVect>>>();
    mats->reserve(s_max);

    for (int j = 1; j <= s_max; ++j) {
        const uint64_t pj = polys_[j - 1];
        const int ej = degs_[j - 1];
        const int r_max = (m - 1) / ej + 1;     // largest Q+1 we ever read
        if (r_max * ej > 63) {
            std::ostringstream oss;
            oss << "NiederreiterF2Gen: polynomial degree " << (r_max * ej)
                << " (= r_max * e_j) exceeds uint64_t arithmetic limit. "
                   "Try a smaller m / s_max combination.";
            throw std::invalid_argument(oss.str());
        }

        // p_j^r for r = 1..r_max.
        std::vector<uint64_t> p_pow(r_max + 1, 0);
        p_pow[1] = pj;
        for (int r = 2; r <= r_max; ++r) {
            p_pow[r] = poly_mul_gf2(p_pow[r - 1], pj);
        }

        // Per-r Laurent series of 1/p_j(x)^r over F_2((x⁻¹)),
        // computed up to index k_max + m  =  (e_j-1) + m.
        const int laurent_len = (ej - 1) + m;
        std::vector<std::vector<uint8_t>> L(r_max + 1);
        for (int r = 1; r <= r_max; ++r) {
            int d = r * ej;
            L[r] = laurent_inverse(p_pow[r], d, laurent_len);
        }

        // Build C_j: row i (1-based), column ℓ (0-based) =
        //   a^{(j)}(Q+1, k, ℓ) = L[Q+1][k + ℓ + 1]
        // where i - 1 = Q * e_j + k, 0 ≤ k < e_j.
        std::vector<BitVect> Cj(m, BitVect(m));
        for (int i = 1; i <= m; ++i) {
            int Q = (i - 1) / ej;
            int k = (i - 1) % ej;
            int r = Q + 1;
            for (int ell = 0; ell < m; ++ell) {
                int idx = k + ell + 1;
                if (idx < static_cast<int>(L[r].size()) && L[r][idx]) {
                    Cj[i - 1].set_bit(ell, 1);
                }
            }
        }
        mats->push_back(std::move(Cj));
    }
    matrices_ = std::shared_ptr<const std::vector<std::vector<BitVect>>>(
        std::move(mats));
}

const std::vector<BitVect>& NiederreiterF2Gen::generating_matrix(int j) const {
    return (*matrices_)[j - 1];
}

std::string NiederreiterF2Gen::display_str() const {
    std::ostringstream oss;
    oss << "NiederreiterF2Gen(m=" << m() << ", s_max=" << s_max() << ")";
    return oss.str();
}

std::unique_ptr<Generator> NiederreiterF2Gen::copy() const {
    return std::unique_ptr<Generator>(new NiederreiterF2Gen(*this));
}

}  // namespace regpoly::core
