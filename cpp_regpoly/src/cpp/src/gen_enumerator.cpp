#include "gen_enumerator.h"

#include <stdexcept>

#include <NTL/ZZ.h>

namespace {

// C(n, k) with NTL::ZZ to avoid overflow at large k.  Returns 0 when
// k < 0 or k > n.
NTL::ZZ binomial(int n, int k) {
    if (k < 0 || k > n) return NTL::ZZ(0);
    if (k == 0 || k == n) return NTL::ZZ(1);
    // Choose the smaller half for speed.
    if (k > n - k) k = n - k;
    NTL::ZZ num(1);
    NTL::ZZ den(1);
    for (int i = 0; i < k; ++i) {
        num *= (n - i);
        den *= (i + 1);
    }
    return num / den;
}

} // anonymous namespace

std::vector<int> unrank_combination(int n, int k, const NTL::ZZ& idx) {
    // Lexicographic unranking via the combinatorial number system.
    // Produces a strictly ascending subset of {0, 1, ..., n-1} of
    // size k at rank `idx` (0-based), matching Python's
    // itertools.combinations ordering.
    if (k < 0 || k > n)
        throw std::invalid_argument("unrank_combination: invalid k");
    NTL::ZZ total = binomial(n, k);
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
            NTL::ZZ below = binomial(n - candidate - 1, remaining_k);
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
