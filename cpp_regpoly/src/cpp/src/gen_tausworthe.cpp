#include "gen_tausworthe.h"
#include <algorithm>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <cstdio>

Tausworthe::Tausworthe(int k, const std::vector<int>& Q, int s, bool quicktaus, int L)
    : Generateur(k, L), Q_(Q), NbCoeff_((int)Q.size()),
      s_(s), quicktaus_(quicktaus), gen_kms_(k - s)
{
    state_ = BitVect(L);
}

std::string Tausworthe::name() const { return "Tausworthe Generator"; }

std::string Tausworthe::display_str() const {
    // Build polynomial string: iterate Q_ from index 0 to NbCoeff_-1
    // If Q_[i]==1 print " x +", else print " x^{Q[i]} +"
    // Then append " 1 "
    std::string poly_str;
    for (int i = 0; i < NbCoeff_; i++) {
        if (Q_[i] == 1)
            poly_str += " x +";
        else
            poly_str += " x^" + std::to_string(Q_[i]) + " +";
    }
    poly_str += " 1 ";

    // Right-pad to 40 chars
    while ((int)poly_str.size() < 40)
        poly_str += " ";

    // Append "   s={s_}"
    poly_str += "   s=" + std::to_string(s_);

    return poly_str;
}

void Tausworthe::init(const BitVect& init_bv) {
    state_ = BitVect(L_);
    state_.copy_part_from(init_bv, k_);

    for (int j = k_; j < L_; j++) {
        int bit = 0;
        for (int i = 0; i < NbCoeff_ - 1; i++)
            bit ^= state_.get_bit(j - (k_ - Q_[i]));
        if (bit & 1)
            state_.set_bit(j, 1);
    }
}

void Tausworthe::next() {
    if (quicktaus_)
        next_quick();
    else
        next_general();
}

std::unique_ptr<Generateur> Tausworthe::copy() const {
    auto p = std::make_unique<Tausworthe>(k_, Q_, s_, quicktaus_, L_);
    p->state_ = state_.copy();
    p->gen_kms_ = gen_kms_;
    return p;
}

BitVect Tausworthe::char_poly() const {
    // The characteristic polynomial is z^k + z^Q[0] + z^Q[1] + ... + 1
    // Q_ contains the non-leading, non-constant exponents, plus k itself.
    // Q_ is sorted. Q_.back() == k_.
    // The constant term (z^0) is always present.
    // Return k bits where bit j = coefficient of z^j (no leading term).
    BitVect bv(k_);
    bv.set_bit(0, 1);  // constant term z^0
    for (int i = 0; i < NbCoeff_ - 1; i++)
        bv.set_bit(Q_[i], 1);
    return bv;
}

void Tausworthe::get_transition_state(uint64_t* out_words, int out_nwords) const {
    BitVect tmp(k_);
    tmp.copy_part_from(state_, k_);
    int n = std::min(out_nwords, tmp.nwords());
    for (int i = 0; i < n; i++)
        out_words[i] = tmp.data()[i];
    for (int i = n; i < out_nwords; i++)
        out_words[i] = 0;
}

void Tausworthe::next_quick() {
    BitVect gen_B(L_);
    gen_B.zero();

    for (int j = 1; j < NbCoeff_ - 1; j++) {
        BitVect shifted = state_.copy();
        shifted.lshift(Q_[j]);
        gen_B.xor_with(shifted);
    }
    gen_B.xor_with(state_);
    gen_B.rshift(gen_kms_);

    state_.and_mask(k_);
    state_.lshift(s_);
    state_.xor_with(gen_B);
}

void Tausworthe::next_general() {
    int ss = std::min(s_, L_ - k_);
    int m = s_;

    while (m > 0) {
        state_.lshift(ss);
        for (int j = L_ - ss; j < L_; j++) {
            int bit = 0;
            for (int i = 0; i < NbCoeff_ - 1; i++)
                bit ^= state_.get_bit(j - (k_ - Q_[i]));
            if (bit & 1)
                state_.set_bit(j, 1);
        }
        m -= ss;
        if (m < ss)
            ss = m;
    }
}

// ── Factory methods ────────────────────────────────────────────────────

std::unique_ptr<Generateur> Tausworthe::from_params(const Params& params, int L) {
    bool quicktaus = params.get_bool("quicktaus", true);

    std::vector<int> Q;
    int k;

    if (params.has("poly")) {
        // Either legacy (poly only) or shape-based with Python having
        // already sampled poly via the tausworthe_poly rand_type.
        auto poly_list = params.get_int_vec("poly");
        Q.assign(poly_list.begin(), poly_list.end());
        std::sort(Q.begin(), Q.end());
        if (Q.empty())
            throw std::invalid_argument("Tausworthe: poly is empty");
        k = Q.back();
    } else {
        // No poly in hand — caller must have given k (nb_terms defaults
        // to 3).  Sample a polynomial and s here.
        k = (int)params.get_int("k");
        if (k <= 0)
            throw std::invalid_argument(
                "Tausworthe: supply `poly`, or `k` (nb_terms optional, "
                "defaults to 3)");
        int nb_terms = 3;
        if (params.has("nb_terms")) {
            int v = (int)params.get_int("nb_terms");
            if (v > 0) nb_terms = v;
        }
        int s_hint = params.has("s") ? (int)params.get_int("s") : 0;
        Q = Tausworthe::random_poly(k, nb_terms, quicktaus, L, s_hint);
    }

    int s;
    if (params.has("s")) {
        s = (int)params.get_int("s");
    } else if (quicktaus && Q.size() >= 2) {
        s = k - Q[Q.size() - 2];
    } else {
        s = 1;
    }

    if (!Tausworthe::is_admissible(Q, s, quicktaus, L))
        throw std::runtime_error(
            "Tausworthe polynomial violates the quicktaus "
            "admissibility inequality (0 < s <= k - q_{t-2} < k <= L)");

    return std::make_unique<Tausworthe>(k, Q, s, quicktaus, L);
}

std::vector<ParamSpec> Tausworthe::param_specs() {
    // `k` and `nb_terms` are logically structural, but marked
    // has_default=true so Python's fill_params doesn't reject a call
    // that provides `poly` directly (legacy form).  C++ from_params
    // enforces the real contract: supply either `poly`, or both `k`
    // and `nb_terms`.
    return {
        {"k",         "int",     true,  true,  0, "",                "",            false},
        {"nb_terms",  "int",     true,  true,  0, "",                "",            false},
        {"quicktaus", "bool",    true,  true,  1, "",                "",            false},
        {"s",         "int",     false, true,  0, "",                "",            false},
        {"poly",      "int_vec", false, false, 0, "tausworthe_poly", "k,nb_terms",  false},
    };
}

// ── Polynomial admissibility & sampler ─────────────────────────────────

bool Tausworthe::is_admissible(
    const std::vector<int>& Q,
    int s, bool quicktaus, int L)
{
    const int n = (int)Q.size();
    if (n < 3) return false;            // nb_terms must be >= 3
    if (n % 2 == 0) return false;       // nb_terms must be odd
    if (Q.front() != 0) return false;
    const int k = Q.back();
    if (k <= 0) return false;

    // Exponents must be strictly increasing.
    for (int i = 1; i < n; i++)
        if (Q[i] <= Q[i - 1]) return false;

    if (!quicktaus)
        return s >= 1 && s < k;

    // quicktaus inequality: 0 < s <= k - q_{t-2} < k <= L
    const int q = Q[n - 2];
    return s >= 1 && q >= 1 && q < k && s <= k - q && k <= L;
}

std::vector<int> Tausworthe::random_poly(
    int k, int nb_terms, bool quicktaus, int L, int s)
{
    if (nb_terms < 3)
        throw std::invalid_argument(
            "Tausworthe::random_poly: nb_terms must be >= 3");
    if (nb_terms % 2 == 0)
        throw std::invalid_argument(
            "Tausworthe::random_poly: nb_terms must be odd");
    if (k <= 0)
        throw std::invalid_argument(
            "Tausworthe::random_poly: k must be positive");
    if (nb_terms - 2 > k - 1)
        throw std::invalid_argument(
            "Tausworthe::random_poly: nb_terms too large for k");

    if (quicktaus && k > L)
        throw std::invalid_argument(
            "Tausworthe::random_poly: quicktaus requires k <= L");

    const int t_interior = nb_terms - 2;   // exponents strictly in (0, k)

    // Thread-local PRNG — shared seed source across calls but per-thread
    // state for safety inside worker processes.
    thread_local std::mt19937_64 rng(std::random_device{}());

    // Sample exactly `count` distinct ints uniformly from [lo, hi].
    // Requires (hi - lo + 1) >= count (caller's responsibility).
    auto sample_distinct =
        [&](int lo, int hi, int count) -> std::vector<int> {
            std::vector<int> out;
            out.reserve(count);
            std::unordered_set<int> taken;
            taken.reserve(count * 2);
            std::uniform_int_distribution<int> dist(lo, hi);
            while ((int)out.size() < count) {
                int x = dist(rng);
                if (taken.insert(x).second)
                    out.push_back(x);
            }
            return out;
        };

    std::vector<int> interior;
    interior.reserve(t_interior);

    if (quicktaus) {
        // The binding tap q_{t-2} must satisfy q_{t-2} <= k - s.  When
        // s == 0 the caller wants auto-derivation (s = k - q_{t-2}),
        // so q_{t-2} may be anywhere in [1, k-1].
        const int q_max = (s > 0) ? (k - s) : (k - 1);
        // We need at least (t_interior - 1) distinct values strictly
        // below q_{t-2}; so q_{t-2} >= t_interior.  Combined with
        // q_{t-2} <= q_max, the combination is admissible only if
        // q_max >= t_interior.
        if (q_max < t_interior)
            throw std::invalid_argument(
                "Tausworthe::random_poly: nb_terms too large for the "
                "quicktaus bound on q_{t-2} (k - s)");

        std::uniform_int_distribution<int> top_dist(t_interior, q_max);
        int q_top = top_dist(rng);
        interior.push_back(q_top);
        if (t_interior >= 2) {
            auto below = sample_distinct(1, q_top - 1, t_interior - 1);
            interior.insert(interior.end(), below.begin(), below.end());
        }
    } else {
        interior = sample_distinct(1, k - 1, t_interior);
    }

    std::sort(interior.begin(), interior.end());

    std::vector<int> Q;
    Q.reserve(nb_terms);
    Q.push_back(0);
    Q.insert(Q.end(), interior.begin(), interior.end());
    Q.push_back(k);
    return Q;
}
