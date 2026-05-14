// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "gen_enumerator.h"

#include <sstream>
#include <stdexcept>

#include <NTL/ZZ.h>

NTL::ZZ binomial_zz(int n, int k) {
    if (k < 0 || k > n) return NTL::ZZ(0);
    if (k == 0 || k == n) return NTL::ZZ(1);
    if (k > n - k) k = n - k;
    NTL::ZZ num(1);
    NTL::ZZ den(1);
    for (int i = 0; i < k; ++i) {
        num *= (n - i);
        den *= (i + 1);
    }
    return num / den;
}

NTL::ZZ parse_zz(const std::string& dec) {
    NTL::ZZ z;
    std::istringstream in(dec);
    in >> z;
    if (in.fail())
        throw std::invalid_argument("parse_zz: not a decimal integer");
    return z;
}

std::string zz_to_dec(const NTL::ZZ& z) {
    std::ostringstream os;
    os << z;
    return os.str();
}

std::vector<int> unrank_combination(int n, int k, const NTL::ZZ& idx) {
    // Lexicographic unranking via the combinatorial number system.
    // Produces a strictly ascending subset of {0, 1, ..., n-1} of
    // size k at rank `idx` (0-based), matching Python's
    // itertools.combinations ordering.
    if (k < 0 || k > n)
        throw std::invalid_argument("unrank_combination: invalid k");
    NTL::ZZ total = binomial_zz(n, k);
    if (idx < 0 || idx >= total)
        throw std::out_of_range("unrank_combination: idx out of range");

    std::vector<int> out;
    out.reserve(k);

    NTL::ZZ remaining = idx;
    int start = 0;
    for (int chosen = 0; chosen < k; ++chosen) {
        // Number of subsets we still need to pick.
        int remaining_k = k - chosen - 1;
        for (int candidate = start; candidate <= n - remaining_k - 1; ++candidate) {
            NTL::ZZ below = binomial_zz(n - candidate - 1, remaining_k);
            if (remaining < below) {
                out.push_back(candidate);
                start = candidate + 1;
                break;
            }
            remaining -= below;
        }
    }
    return out;
}

std::vector<NTL::ZZ> mixed_radix_decode(
    const std::vector<NTL::ZZ>& sizes, const NTL::ZZ& idx)
{
    // The first axis varies slowest.  Compute suffix products then
    // divmod from the outside in.
    const size_t n = sizes.size();
    std::vector<NTL::ZZ> suffix(n + 1);
    suffix[n] = NTL::ZZ(1);
    for (int i = (int)n - 1; i >= 0; --i)
        suffix[i] = suffix[i + 1] * sizes[i];

    if (idx < 0 || idx >= suffix[0])
        throw std::out_of_range("mixed_radix_decode: idx out of range");

    std::vector<NTL::ZZ> out(n);
    NTL::ZZ remaining = idx;
    for (size_t i = 0; i < n; ++i) {
        // Quotient by the product of the axes still to come.
        NTL::ZZ q, r;
        DivRem(q, r, remaining, suffix[i + 1]);
        out[i] = q;
        remaining = r;
    }
    return out;
}
