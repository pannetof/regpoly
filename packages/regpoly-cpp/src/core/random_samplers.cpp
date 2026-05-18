// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "random_samplers.h"

#include "primitivity.h"
#include "tausworthe.h"

#include <algorithm>
#include <cctype>
#include <random>
#include <set>
#include <stdexcept>
#include <string>

using namespace regpoly::core;
using namespace regpoly::internal;
using namespace regpoly::random;


namespace regpoly::random {

namespace {

thread_local std::mt19937_64 tls_rng{std::random_device{}()};

std::string strip(const std::string& s) {
    size_t a = 0;
    while (a < s.size() && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    size_t b = s.size();
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return s.substr(a, b - a);
}

bool parse_literal_int(const std::string& s, int64_t& out) {
    if (s.empty()) return false;
    size_t i = (s[0] == '-') ? 1 : 0;
    if (i == s.size()) return false;
    for (; i < s.size(); ++i)
        if (!std::isdigit(static_cast<unsigned char>(s[i]))) return false;
    try {
        out = std::stoll(s);
        return true;
    } catch (...) {
        return false;
    }
}

int64_t fetch_int_param(const std::string& name, const Params& p) {
    auto it = p.ints().find(name);
    if (it == p.ints().end())
        throw std::invalid_argument(
            "random sampler: parameter '" + name + "' is not in the params bag");
    return it->second;
}

uint64_t random_bits(int nbits) {
    if (nbits <= 0) return 0;
    if (nbits >= 64) {
        return tls_rng();
    }
    uint64_t mask = (nbits == 64) ? ~uint64_t(0) : ((uint64_t(1) << nbits) - 1);
    return tls_rng() & mask;
}

}  // namespace

void seed(uint64_t s) {
    if (s == 0) {
        tls_rng.seed(std::random_device{}());
    } else {
        tls_rng.seed(s);
    }
}

std::mt19937_64& rng() { return tls_rng; }

int64_t eval_param_expr(const std::string& expr_in, const Params& params) {
    std::string expr = strip(expr_in);

    int64_t lit = 0;
    if (parse_literal_int(expr, lit)) return lit;

    // Look for + or -, but ignore a leading '-' (handled above as a
    // literal). Operators inside the expression must split at the
    // first occurrence past index 0.
    for (char op : {'-', '+'}) {
        size_t pos = expr.find(op, 1);
        if (pos != std::string::npos) {
            std::string name = strip(expr.substr(0, pos));
            std::string off_s = strip(expr.substr(pos + 1));
            int64_t off = 0;
            if (!parse_literal_int(off_s, off))
                throw std::invalid_argument(
                    "eval_param_expr: cannot parse offset in '" + expr_in + "'");
            int64_t base = fetch_int_param(name, params);
            return op == '-' ? base - off : base + off;
        }
    }

    return fetch_int_param(expr, params);
}

bool sample_generic_into(const ParamSpec& spec, Params& params) {
    const std::string& rt = spec.rand_type;
    const std::string& ra = spec.rand_args;

    if (rt == "bitmask") {
        int64_t bits = eval_param_expr(ra, params);
        params.set_int(spec.name,
                       static_cast<int64_t>(random_bits(static_cast<int>(bits))));
        return true;
    }

    if (rt == "range") {
        size_t comma = ra.find(',');
        if (comma == std::string::npos)
            throw std::invalid_argument(
                "range sampler: rand_args must be 'lo,hi' (got '" + ra + "')");
        int64_t lo = eval_param_expr(ra.substr(0, comma), params);
        int64_t hi = eval_param_expr(ra.substr(comma + 1), params);
        if (hi < lo)
            throw std::invalid_argument(
                "range sampler: hi < lo for spec '" + spec.name + "'");
        std::uniform_int_distribution<int64_t> dist(lo, hi);
        params.set_int(spec.name, dist(tls_rng));
        return true;
    }

    if (rt == "poly_exponents") {
        int64_t k = eval_param_expr(ra, params);
        if (k < 2)
            throw std::invalid_argument(
                "poly_exponents sampler: k must be >= 2 (got "
                + std::to_string(k) + ")");
        int64_t hi = std::min<int64_t>(k - 1, 10);
        std::uniform_int_distribution<int64_t> n_dist(1, hi);
        int64_t n = n_dist(tls_rng);

        // Sample n distinct integers from [1, k-1] uniformly without
        // replacement. Mirrors Python's random.sample.
        std::set<int64_t> chosen;
        std::uniform_int_distribution<int64_t> dist(1, k - 1);
        while (static_cast<int64_t>(chosen.size()) < n)
            chosen.insert(dist(tls_rng));

        std::vector<int> out;
        out.reserve(n + 1);
        for (auto v : chosen) out.push_back(static_cast<int>(v));
        std::sort(out.begin(), out.end());
        out.push_back(0);
        params.set_int_vec(spec.name, out);
        return true;
    }

    if (rt == "irreducible_gf2") {
        // Sample a degree-w polynomial over F_2 by picking the lower w
        // bits uniformly and asserting `z^w + lower` is irreducible.
        // For w = 32 the density of irreducibles is ~1/w, so ~32 draws
        // expected per success.
        int64_t w_bits = eval_param_expr(ra, params);
        if (w_bits <= 0 || w_bits > 63)
            throw std::invalid_argument(
                "irreducible_gf2 sampler: w must be in [1, 63] (got "
                + std::to_string(w_bits) + ")");
        // Bounded retry budget to surface bugs in misuse rather than
        // hang forever; 64 * w typically suffices by a wide margin.
        const int max_tries = 64 * static_cast<int>(w_bits);
        for (int i = 0; i < max_tries; ++i) {
            uint64_t cand = random_bits(static_cast<int>(w_bits));
            if (is_irreducible_gf2w_modM(cand, static_cast<int>(w_bits))) {
                params.set_int(spec.name, static_cast<int64_t>(cand));
                return true;
            }
        }
        throw std::runtime_error(
            "irreducible_gf2 sampler: failed to draw an irreducible "
            "polynomial in " + std::to_string(max_tries)
            + " tries (w=" + std::to_string(w_bits) + ")");
    }

    if (rt == "bitmask_vec") {
        size_t comma = ra.find(',');
        if (comma == std::string::npos)
            throw std::invalid_argument(
                "bitmask_vec sampler: rand_args must be "
                "'bits_param,length_param' (got '" + ra + "')");
        int64_t bits = eval_param_expr(ra.substr(0, comma), params);
        std::string len_param = strip(ra.substr(comma + 1));
        // Length resolution: try as a vector-param name first (size of
        // the existing int_vec/uint_vec). If that fails, fall back to
        // a scalar expression (parameter or literal). This lets paper-
        // notation specs use `coeff` with rand_args="w,m" — where `m`
        // is the scalar polynomial-term-count (2 or 3).
        size_t length = 0;
        auto it_iv = params.int_vecs().find(len_param);
        auto it_uv = params.uint_vecs().find(len_param);
        if (it_iv != params.int_vecs().end()) {
            length = it_iv->second.size();
        } else if (it_uv != params.uint_vecs().end()) {
            length = it_uv->second.size();
        } else {
            // Fall back to scalar expression resolution.
            try {
                int64_t n = eval_param_expr(len_param, params);
                if (n < 0)
                    throw std::invalid_argument(
                        "bitmask_vec sampler: scalar length '"
                        + len_param + "' resolved to negative value");
                length = static_cast<size_t>(n);
            } catch (const std::invalid_argument&) {
                throw std::invalid_argument(
                    "bitmask_vec sampler: length-source '"
                    + len_param
                    + "' is neither a vector nor a resolvable scalar");
            }
        }

        std::vector<uint64_t> out;
        out.reserve(length);
        for (size_t i = 0; i < length; ++i)
            out.push_back(random_bits(static_cast<int>(bits)));
        params.set_uint_vec(spec.name, out);
        return true;
    }

    return false;
}

bool sample_param_into(const ParamSpec& spec, Params& params, int L) {
    if (sample_generic_into(spec, params)) return true;

    const std::string& rt = spec.rand_type;
    if (rt == "tausworthe_s" || rt == "tausworthe_poly") {
        auto r = TauswortheGen::generate_random(rt, spec.rand_args, params, L);
        if (r.is_vec)
            params.set_int_vec(spec.name, r.vec_val);
        else
            params.set_int(spec.name, r.int_val);
        for (const auto& kv : r.side_ints)
            params.set_int(kv.first, kv.second);
        return true;
    }

    return false;
}

}  // namespace regpoly::random
