#include "random_samplers.h"

#include <algorithm>
#include <cctype>
#include <random>
#include <set>
#include <stdexcept>
#include <string>

namespace regpoly_random {

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

    if (rt == "bitmask_vec") {
        size_t comma = ra.find(',');
        if (comma == std::string::npos)
            throw std::invalid_argument(
                "bitmask_vec sampler: rand_args must be "
                "'bits_param,length_param' (got '" + ra + "')");
        int64_t bits = eval_param_expr(ra.substr(0, comma), params);
        std::string len_param = strip(ra.substr(comma + 1));
        // Length is the .size() of the existing vector parameter named
        // `len_param`. Mirrors Python's `len(params[length_param])`.
        size_t length = 0;
        auto it_iv = params.int_vecs().find(len_param);
        auto it_uv = params.uint_vecs().find(len_param);
        if (it_iv != params.int_vecs().end())
            length = it_iv->second.size();
        else if (it_uv != params.uint_vecs().end())
            length = it_uv->second.size();
        else
            throw std::invalid_argument(
                "bitmask_vec sampler: length-source parameter '"
                + len_param + "' not found in params");

        std::vector<uint64_t> out;
        out.reserve(length);
        for (size_t i = 0; i < length; ++i)
            out.push_back(random_bits(static_cast<int>(bits)));
        params.set_uint_vec(spec.name, out);
        return true;
    }

    return false;
}

}  // namespace regpoly_random
