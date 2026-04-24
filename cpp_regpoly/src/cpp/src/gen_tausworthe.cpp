#include "gen_tausworthe.h"
#include "gen_enumerator.h"

#include <algorithm>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <cstdio>

#include <NTL/ZZ.h>

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

        // When the caller also supplied structural hints, make sure
        // they're consistent with the polynomial itself.
        if (params.has("k")) {
            int k_user = (int)params.get_int("k");
            if (k_user > 0 && k_user != k) {
                throw std::invalid_argument(
                    "Tausworthe: k=" + std::to_string(k_user) +
                    " disagrees with max(poly)=" + std::to_string(k));
            }
        }
        if (params.has("nb_terms")) {
            int nbt_user = (int)params.get_int("nb_terms");
            if (nbt_user > 0 && nbt_user != (int)Q.size()) {
                throw std::invalid_argument(
                    "Tausworthe: nb_terms=" + std::to_string(nbt_user) +
                    " disagrees with len(poly)=" + std::to_string(Q.size()));
            }
        }
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
    //
    // `s` carries rand_type="tausworthe_s" so the UI shows a Random
    // checkbox for it.  When the user fixes `poly` but leaves `s`
    // random, each try draws a fresh admissible s coprime with
    // 2^k - 1.
    return {
        {"k",         "int",     true,  true,  0, "",                "",            false},
        {"nb_terms",  "int",     true,  true,  0, "",                "",            false},
        {"quicktaus", "bool",    true,  true,  1, "",                "",            false},
        {"s",         "int",     false, false, 0, "tausworthe_s",    "k",           false},
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

// Pick an admissible decimation step s in [1, s_max] with
// gcd(s, 2^k - 1) = 1.  Throws std::invalid_argument on s_max < 1 or
// after 256 failed tries.  Used by both rand_types handled in
// Tausworthe::generate_random.
static int pick_admissible_s(int k, int s_max, std::mt19937_64& rng) {
    if (s_max < 1)
        throw std::invalid_argument(
            "Tausworthe::generate_random: no admissible s (s_max < 1)");
    const NTL::ZZ period = NTL::power2_ZZ(k) - 1;
    std::uniform_int_distribution<int> dist(1, s_max);
    for (int tries = 0; tries < 256; ++tries) {
        int s = dist(rng);
        if (NTL::IsOne(NTL::GCD(period, NTL::ZZ(s)))) return s;
    }
    throw std::invalid_argument(
        "Tausworthe::generate_random: could not find s in [1, s_max] "
        "coprime with 2^k - 1 after 256 tries");
}

RandomParamResult Tausworthe::generate_random(
    const std::string& rand_type,
    const std::string& /*rand_args*/,
    const Params& params, int L)
{
    thread_local std::mt19937_64 rng(std::random_device{}());

    if (rand_type == "tausworthe_s") {
        const bool quicktaus = params.get_bool("quicktaus", true);
        const bool has_poly = params.has("poly");

        // Derive k from the explicit param when present, otherwise
        // from the largest exponent of the caller-supplied poly.
        int k = 0;
        if (params.get_int("k", 0) > 0) {
            k = (int)params.get_int("k");
        } else if (has_poly) {
            for (int x : params.get_int_vec("poly"))
                if (x > k) k = x;
        } else {
            throw std::invalid_argument(
                "Tausworthe: randomizing s needs k or poly");
        }

        int s_max;
        if (quicktaus && has_poly) {
            // Largest interior exponent q_{t-2}: s <= k - q_{t-2}
            const auto poly = params.get_int_vec("poly");
            int q_top = 0;
            for (int x : poly)
                if (x > 0 && x < k && x > q_top) q_top = x;
            if (q_top == 0)
                throw std::invalid_argument(
                    "Tausworthe: poly has no interior exponents");
            s_max = k - q_top;
        } else if (quicktaus) {
            int nb_terms = (int)params.get_int("nb_terms", 0);
            if (nb_terms <= 0 && has_poly)
                nb_terms = (int)params.get_int_vec("poly").size();
            if (nb_terms <= 0) nb_terms = 3;
            s_max = k - (nb_terms - 2);
        } else {
            s_max = k - 1;
        }

        RandomParamResult out;
        out.int_val = pick_admissible_s(k, s_max, rng);
        return out;
    }

    if (rand_type == "tausworthe_poly") {
        const int k = (int)params.get_int("k", 0);
        if (k <= 0)
            throw std::invalid_argument(
                "Tausworthe: randomizing poly needs k > 0");
        int nb_terms = (int)params.get_int("nb_terms", 0);
        if (nb_terms <= 0) nb_terms = 3;
        const bool quicktaus = params.get_bool("quicktaus", true);
        int s = (int)params.get_int("s", 0);

        RandomParamResult out;
        if (s <= 0) {
            const int s_max = quicktaus ? (k - (nb_terms - 2)) : (k - 1);
            s = pick_admissible_s(k, s_max, rng);
            out.side_ints["s"] = s;   // propagate to caller
        }

        out.is_vec = true;
        out.vec_val = Tausworthe::random_poly(
            k, nb_terms, quicktaus, L, s);
        return out;
    }

    throw std::invalid_argument(
        "Tausworthe::generate_random: unsupported rand_type '"
        + rand_type + "'");
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

// ═══════════════════════════════════════════════════════════════════════════
// Exhaustive enumerator
// ═══════════════════════════════════════════════════════════════════════════

namespace {

// C(n, k) as NTL::ZZ.  Returns 0 for invalid (k < 0 || k > n).
NTL::ZZ binomial_zz(int n, int k) {
    if (k < 0 || k > n) return NTL::ZZ(0);
    if (k == 0 || k == n) return NTL::ZZ(1);
    if (k > n - k) k = n - k;
    NTL::ZZ num(1), den(1);
    for (int i = 0; i < k; ++i) {
        num *= (n - i);
        den *= (i + 1);
    }
    return num / den;
}

// Admissible decimation values: s in [1, s_max] with gcd(s, 2^k-1) = 1.
std::vector<int> admissible_s_values_(int k, bool quicktaus,
                                      int t_interior) {
    int s_max = quicktaus ? (k - t_interior) : (k - 1);
    if (s_max < 1) return {};

    // period = 2^k - 1
    NTL::ZZ period = NTL::power2_ZZ(k) - 1;

    std::vector<int> out;
    out.reserve(s_max);
    for (int s = 1; s <= s_max; ++s) {
        NTL::ZZ g = NTL::GCD(period, NTL::ZZ(s));
        if (NTL::IsOne(g)) out.push_back(s);
    }
    return out;
}

// Upper bound on q_{t-2} for a given s.
//   quicktaus=true:  q_top <= k - s
//   quicktaus=false: q_top <= k - 1
int q_top_max_(int k, bool quicktaus, int s) {
    return quicktaus ? (k - s) : (k - 1);
}

// Concrete enumerator for Tausworthe (s × poly).  Enumeration order:
//   outer axis = s    (varies slowest)
//   inner axis = poly (q_top, then unrank the interior subset)
class TausworthePolyEnumerator : public GenEnumerator {
public:
    TausworthePolyEnumerator(int k, int nb_terms, bool quicktaus, int L,
                             std::vector<int> s_values)
        : k_(k), nb_terms_(nb_terms), quicktaus_(quicktaus), L_(L),
          s_values_(std::move(s_values))
    {
        const int t_interior = nb_terms_ - 2;
        const int q_top_lo = std::max(1, t_interior);

        // Per-s: cumulative subset counts across q_top values.  An
        // "entry" is (s_index, q_top, local_index into
        // C(q_top-1, t_interior-1) subsets).
        per_s_.reserve(s_values_.size());
        NTL::ZZ total(0);
        s_offsets_.reserve(s_values_.size() + 1);
        s_offsets_.push_back(total);
        for (int s : s_values_) {
            const int q_top_hi = q_top_max_(k_, quicktaus_, s);
            SEntry entry;
            entry.q_top_lo = q_top_lo;
            entry.cum_by_q_top.reserve(
                std::max(0, q_top_hi - q_top_lo + 2));
            NTL::ZZ acc(0);
            entry.cum_by_q_top.push_back(acc);
            for (int q_top = q_top_lo; q_top <= q_top_hi; ++q_top) {
                acc += binomial_zz(q_top - 1, t_interior - 1);
                entry.cum_by_q_top.push_back(acc);
            }
            entry.size = acc;
            per_s_.push_back(std::move(entry));
            total += per_s_.back().size;
            s_offsets_.push_back(total);
        }
        total_ = total;
    }

    std::string size_dec() const override {
        std::ostringstream os;
        os << total_;
        return os.str();
    }

    std::vector<Axis> axes() const override {
        std::vector<Axis> out;
        Axis s_ax;
        s_ax.name = "s";
        {
            std::ostringstream os;
            os << s_values_.size();
            s_ax.size_dec = os.str();
        }
        if (s_values_.empty()) {
            s_ax.describe = "admissible decimation steps";
        } else if (s_values_.size() == 1) {
            std::ostringstream d;
            d << "fixed s=" << s_values_[0];
            s_ax.describe = d.str();
        } else {
            std::ostringstream d;
            d << "s in [1, " << s_values_.back()
              << "] coprime with 2^" << k_ << "-1";
            s_ax.describe = d.str();
        }
        out.push_back(std::move(s_ax));

        Axis poly_ax;
        poly_ax.name = "poly";
        // The poly axis size varies with s; report the total subset
        // count / |s| as a rough describe.  For the UI summary we
        // prefer the accurate product: total / |s_values|.
        NTL::ZZ poly_size(0);
        if (!s_values_.empty()) {
            poly_size = total_ / NTL::ZZ((long)s_values_.size());
        }
        {
            std::ostringstream os;
            os << poly_size;
            poly_ax.size_dec = os.str();
        }
        {
            std::ostringstream d;
            d << "admissible polynomials with k=" << k_
              << ", nb_terms=" << nb_terms_
              << (quicktaus_ ? " (quicktaus)" : "");
            poly_ax.describe = d.str();
        }
        out.push_back(std::move(poly_ax));
        return out;
    }

    Params at(const std::string& idx_dec) const override {
        NTL::ZZ idx;
        {
            std::istringstream is(idx_dec);
            is >> idx;
        }
        if (idx < 0 || idx >= total_)
            throw std::out_of_range("Tausworthe enumerator: idx out of range");

        // Binary-search s_offsets_ to locate the s axis.
        // s_offsets_[i] = cumulative count up to (but excluding) s_index i.
        // Find the greatest i such that s_offsets_[i] <= idx.
        size_t lo = 0, hi = per_s_.size();
        while (lo + 1 < hi) {
            size_t mid = lo + (hi - lo) / 2;
            if (s_offsets_[mid] <= idx) lo = mid;
            else hi = mid;
        }
        const size_t s_index = lo;
        const int s = s_values_[s_index];
        const SEntry& entry = per_s_[s_index];
        NTL::ZZ within_s = idx - s_offsets_[s_index];

        // Binary-search entry.cum_by_q_top to locate q_top.
        size_t q_lo = 0, q_hi = entry.cum_by_q_top.size() - 1;
        while (q_lo + 1 < q_hi) {
            size_t mid = q_lo + (q_hi - q_lo) / 2;
            if (entry.cum_by_q_top[mid] <= within_s) q_lo = mid;
            else q_hi = mid;
        }
        const int q_top_offset = (int)q_lo;
        const int q_top = entry.q_top_lo + q_top_offset;
        NTL::ZZ within_q = within_s - entry.cum_by_q_top[q_top_offset];

        // Unrank the (nb_terms - 3)-subset of {1, ..., q_top - 1}.
        const int t_interior = nb_terms_ - 2;
        std::vector<int> interior_below;
        if (t_interior >= 2) {
            // unrank_combination picks from {0, ..., q_top - 2}; shift
            // to {1, ..., q_top - 1}.
            auto ks = unrank_combination(q_top - 1, t_interior - 1,
                                         within_q);
            interior_below.reserve(ks.size());
            for (int x : ks) interior_below.push_back(x + 1);
        }

        std::vector<int> poly;
        poly.reserve(nb_terms_);
        poly.push_back(0);
        poly.insert(poly.end(), interior_below.begin(), interior_below.end());
        poly.push_back(q_top);
        poly.push_back(k_);
        // Ensure strictly increasing (they are by construction).
        std::sort(poly.begin(), poly.end());

        Params out;
        out.set_int("k", k_);
        out.set_int("nb_terms", nb_terms_);
        out.set_bool("quicktaus", quicktaus_);
        out.set_int("s", s);
        out.set_int_vec("poly", poly);
        return out;
    }

private:
    struct SEntry {
        int q_top_lo;                      // lower bound on q_top (== max(1, nb_terms-2))
        std::vector<NTL::ZZ> cum_by_q_top; // cum_by_q_top[i] = Σ C(q_top-1, t-1) for q_top_lo..q_top_lo+i-1
        NTL::ZZ size;                      // cum_by_q_top.back()
    };

    int k_;
    int nb_terms_;
    bool quicktaus_;
    int L_;
    std::vector<int> s_values_;
    std::vector<SEntry> per_s_;
    std::vector<NTL::ZZ> s_offsets_;       // prefix sums of per_s_.size
    NTL::ZZ total_;
};

} // anonymous namespace

// ── Public factory ─────────────────────────────────────────────────────────

std::unique_ptr<GenEnumerator> Tausworthe::make_enumerator(
    const Params& resolved, int L)
{
    if (!resolved.has("k"))
        throw std::invalid_argument("needs_k");
    if (!resolved.has("nb_terms"))
        throw std::invalid_argument("needs_nb_terms");

    int k = (int)resolved.get_int("k");
    int nb_terms = (int)resolved.get_int("nb_terms");
    bool quicktaus = resolved.get_bool("quicktaus", true);

    if (k <= 0)
        throw std::invalid_argument("needs_positive_k");
    if (nb_terms < 3 || nb_terms % 2 == 0)
        throw std::invalid_argument("needs_odd_nb_terms");
    if (nb_terms - 2 > k - 1)
        throw std::invalid_argument("needs_smaller_nb_terms");
    if (quicktaus && k > L)
        throw std::invalid_argument("needs_k_le_L_for_quicktaus");

    const int t_interior = nb_terms - 2;
    std::vector<int> s_values;

    if (resolved.has("s")) {
        int s_user = (int)resolved.get_int("s");
        int s_max = quicktaus ? (k - t_interior) : (k - 1);
        if (s_user < 1 || s_user > s_max)
            throw std::invalid_argument("needs_admissible_fixed_s");
        NTL::ZZ period = NTL::power2_ZZ(k) - 1;
        NTL::ZZ g = NTL::GCD(period, NTL::ZZ(s_user));
        if (!NTL::IsOne(g))
            throw std::invalid_argument("needs_coprime_fixed_s");
        s_values.push_back(s_user);
    } else {
        s_values = admissible_s_values_(k, quicktaus, t_interior);
        if (s_values.empty())
            throw std::invalid_argument("needs_admissible_s_values");
    }

    return std::make_unique<TausworthePolyEnumerator>(
        k, nb_terms, quicktaus, L, std::move(s_values));
}
