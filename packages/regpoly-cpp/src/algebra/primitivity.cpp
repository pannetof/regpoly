#include "primitivity.h"

#include "bitvect.h"
#include "generator.h"

#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/ZZ.h>

#include <algorithm>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace {

const std::unordered_set<int>& mersenne_prime_exponents() {
    static const std::unordered_set<int> set = {
        2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607,
        1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213,
        19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091,
        756839, 859433, 1257787,
    };
    return set;
}

std::vector<int> divisors_of(int n) {
    std::vector<int> out;
    for (int i = 1; static_cast<long long>(i) * i <= n; ++i) {
        if (n % i == 0) {
            out.push_back(i);
            if (i != n / i)
                out.push_back(n / i);
        }
    }
    std::sort(out.begin(), out.end());
    return out;
}

// Build the NTL polynomial f from a BitVect char_poly + leading x^k.
NTL::GF2X char_poly_to_ntl(const BitVect& char_poly_bv, int k) {
    NTL::GF2X f;
    NTL::SetCoeff(f, k);
    for (int j = 0; j < k; ++j) {
        if (char_poly_bv.get_bit(j))
            NTL::SetCoeff(f, j);
    }
    return f;
}

bool ntl_is_irreducible(const NTL::GF2X& f) {
    if (!IsOne(coeff(f, 0))) return false;
    return NTL::IterIrredTest(f) != 0;
}

bool ntl_is_primitive_with_factors(
    const NTL::GF2X& f, int k,
    const std::vector<std::string>& factor_strings)
{
    if (!IsOne(coeff(f, 0))) return false;
    if (NTL::IterIrredTest(f) == 0) return false;

    NTL::GF2XModulus F;
    NTL::build(F, f);
    NTL::ZZ order = NTL::power(NTL::ZZ(2), k) - 1;

    for (const auto& s : factor_strings) {
        NTL::ZZ p = NTL::conv<NTL::ZZ>(s.c_str());
        NTL::ZZ exp = order / p;
        NTL::GF2X r;
        NTL::PowerXMod(r, exp, F);
        if (IsOne(r)) return false;
    }
    return true;
}

}  // namespace

bool is_mersenne_prime_exponent(int k) {
    const auto& s = mersenne_prime_exponents();
    return s.find(k) != s.end();
}

std::optional<std::vector<std::string>> get_primitive_factors_for_k(int k) {
    if (k <= 0) return std::nullopt;

    // Mersenne prime exponent: 2^k - 1 is itself prime (single factor).
    if (is_mersenne_prime_exponent(k)) {
        NTL::ZZ order = NTL::power(NTL::ZZ(2), k) - 1;
        std::ostringstream oss;
        oss << order;
        return std::vector<std::string>{oss.str()};
    }

    // Walk every divisor d > 1 of k and gather Phi_d(2)'s factors.
    std::set<std::string> uniq;  // sorts numerically by string length
                                 // *not* numerically — we sort
                                 // separately at the end.
    for (int d : divisors_of(k)) {
        if (d == 1) continue;
        bool complete = false;
        const auto* facs = regpoly_internal::lookup_factors(d, complete);
        if (facs == nullptr) return std::nullopt;
        if (!complete) return std::nullopt;
        for (const auto& p : *facs)
            uniq.insert(p);
    }

    // Sort numerically. ZZ comparison works for arbitrary precision.
    std::vector<NTL::ZZ> as_zz;
    as_zz.reserve(uniq.size());
    for (const auto& s : uniq)
        as_zz.push_back(NTL::conv<NTL::ZZ>(s.c_str()));
    std::sort(as_zz.begin(), as_zz.end());

    std::vector<std::string> out;
    out.reserve(as_zz.size());
    for (const auto& z : as_zz) {
        std::ostringstream oss;
        oss << z;
        out.push_back(oss.str());
    }
    return out;
}

bool is_full_period(const BitVect& char_poly_bv, int k) {
    NTL::GF2X f = char_poly_to_ntl(char_poly_bv, k);

    if (is_mersenne_prime_exponent(k))
        return ntl_is_irreducible(f);

    auto factors = get_primitive_factors_for_k(k);
    if (!factors)
        throw std::runtime_error(
            "is_full_period: complete factorisation of 2^k - 1 is "
            "unavailable for k = " + std::to_string(k) +
            ". Extend packages/regpoly/src/regpoly/data/"
            "primitive_factors.json and regenerate "
            "primitive_factors_data.cpp.");

    return ntl_is_primitive_with_factors(f, k, *factors);
}

bool is_full_period(const Generator& gen) {
    BitVect cp = gen.char_poly();
    return is_full_period(cp, gen.k());
}
