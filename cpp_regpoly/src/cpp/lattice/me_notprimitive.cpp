#include "me_notprimitive.h"
#include "me_helpers.h"   // test_me_lat (full-period fast path)
                             // + find_polys/normalize_polys (slow path)
#include "dual_lattice.h"    // DualLatticeBase — dual-lattice reduction
#include <algorithm>
#include <climits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/vec_GF2.h>
#include <NTL/GF2.h>

namespace {

// ── Mersenne-prime-exponent fast path (mirrors generateur.cpp) ─────────
bool is_mersenne_prime_exponent(int p) {
    static const int mp[] = {
        2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607,
        1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213,
        19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091,
        756839, 859433, 1257787, 1398269, 2976221, 3021377, 6972593,
        13466917, 20996011, 24036583, 25964951, 30402457, 32582657,
        37156667, 42643801, 43112609, 57885161, 74207281, 77232917,
        82589933, 0
    };
    for (const int* q = mp; *q; q++) if (*q == p) return true;
    return false;
}

// ─────────────────────────────────────────────────────────────────────
// CombinedView: wraps J components + their tempering chains as a
// single Generator of dimension kg = Σ k_j.
// ─────────────────────────────────────────────────────────────────────
class CombinedView : public Generator {
public:
    CombinedView(const std::vector<Generator*>& gens,
                 const std::vector<std::vector<Transformation*>>& trans,
                 int kg, int L)
        : Generator(kg, L), trans_(trans)
    {
        gens_.reserve(gens.size());
        k_per_.reserve(gens.size());
        for (auto* g : gens) {
            gens_.push_back(g->copy());
            k_per_.push_back(g->k());
        }
        sync_state_from_components();
    }

    std::string name()        const override { return "CombinedView"; }
    std::string display_str() const override { return ""; }

    int output_phases() const override {
        long long acc = 1;
        for (const auto& g : gens_) {
            int p = g->output_phases();
            if (p <= 0) p = 1;
            long long g_p = std::__gcd<long long>(acc, p);
            acc = (acc / g_p) * p;
        }
        return (int)acc;
    }

    void init(const BitVect& init_bv) override {
        int offset = 0;
        for (size_t j = 0; j < gens_.size(); j++) {
            int k_j = k_per_[j];
            BitVect per_j(k_j);
            for (int i = 0; i < k_j; i++) {
                if (offset + i < init_bv.nbits() && init_bv.get_bit(offset + i))
                    per_j.set_bit(i, 1);
            }
            gens_[j]->init(per_j);
            offset += k_j;
        }
        sync_state_from_components();
    }

    void next() override {
        for (auto& g : gens_) g->next();
        sync_state_from_components();
    }

    std::unique_ptr<Generator> copy() const override {
        std::vector<Generator*> raw_gens;
        raw_gens.reserve(gens_.size());
        for (auto& g : gens_) raw_gens.push_back(g.get());
        auto c = std::make_unique<CombinedView>(raw_gens, trans_, k_, L_);
        c->state_ = state_.copy();
        return c;
    }

    BitVect get_output() const override {
        BitVect out(L_);
        for (size_t j = 0; j < gens_.size(); j++) {
            BitVect tj = gens_[j]->get_output().copy();
            for (auto* t : trans_[j]) t->apply(tj);
            int n = std::min(L_, tj.nbits());
            for (int b = 0; b < n; b++) {
                if (tj.get_bit(b)) out.set_bit(b, out.get_bit(b) ^ 1);
            }
        }
        return out;
    }

private:
    std::vector<std::unique_ptr<Generator>> gens_;
    std::vector<std::vector<Transformation*>> trans_;
    std::vector<int> k_per_;

    void sync_state_from_components() {
        BitVect s(k_);
        int offset = 0;
        for (size_t j = 0; j < gens_.size(); j++) {
            const BitVect& gj = gens_[j]->state();
            int k_j = k_per_[j];
            for (int i = 0; i < k_j; i++) {
                if (i < gj.nbits() && gj.get_bit(i)) s.set_bit(offset + i, 1);
            }
            offset += k_j;
        }
        state_ = std::move(s);
    }
};

NTL::GF2X gf2x_lcm(const NTL::GF2X& a, const NTL::GF2X& b) {
    if (NTL::IsZero(a)) return b;
    if (NTL::IsZero(b)) return a;
    NTL::GF2X g, ab, q;
    NTL::GCD(g, a, b);
    NTL::mul(ab, a, b);
    NTL::div(q, ab, g);
    return q;
}

BitVect step_once(const Generator& proto, const BitVect& r) {
    auto g = proto.copy();
    g->init(r);
    g->next();
    return g->state().copy();
}

BitVect apply_polynomial(const NTL::GF2X& g, const BitVect& s,
                         const Generator& proto)
{
    long d = NTL::deg(g);
    BitVect r(s.nbits());
    for (long i = d; i >= 0; i--) {
        if (i < d) r = step_once(proto, r);
        if (NTL::IsOne(NTL::coeff(g, i))) r.xor_with(s);
    }
    return r;
}

struct PhiPick {
    NTL::GF2X phi;
    int p;
    bool primitivity_certified;
};

PhiPick select_phi(const NTL::GF2X& chi) {
    NTL::vec_pair_GF2X_long fac;
    NTL::CanZass(fac, chi);

    std::vector<long> order(fac.length());
    for (long i = 0; i < fac.length(); i++) order[i] = i;
    std::sort(order.begin(), order.end(),
              [&](long a, long b) { return NTL::deg(fac[a].a) > NTL::deg(fac[b].a); });

    for (long i : order) {
        long d = NTL::deg(fac[i].a);
        if (d <= 0) continue;
        if (is_mersenne_prime_exponent((int)d) &&
            !NTL::IsZero(NTL::coeff(fac[i].a, 0)))
        {
            return { fac[i].a, (int)d, true };
        }
    }
    for (long i : order) {
        long d = NTL::deg(fac[i].a);
        if (d > 0 && !NTL::IsZero(NTL::coeff(fac[i].a, 0)))
            return { fac[i].a, (int)d, false };
    }
    throw std::runtime_error(
        "test_me_notprimitive: chi has no nontrivial irreducible factor");
}

NTL::GF2X recover_char_poly(const CombinedView& view_proto, int kg) {
    auto build_state_seq = [&](const BitVect& seed_state, int coord, int N_total) {
        NTL::vec_GF2 seq;
        seq.SetLength(N_total);
        BitVect cur = seed_state.copy();
        for (int N = 0; N < N_total; N++) {
            cur = step_once(view_proto, cur);
            if (coord < cur.nbits() && cur.get_bit(coord)) seq.put(N, 1);
        }
        return seq;
    };

    NTL::GF2X chi;
    NTL::set(chi);

    std::mt19937_64 rng(0xA110CAB1EE1ULL);
    BitVect seed(kg);
    for (int i = 0; i < kg; i++) if (rng() & 1ULL) seed.set_bit(i, 1);
    bool nonz = false;
    for (int w = 0; w < seed.nwords(); w++)
        if (seed.data()[w]) { nonz = true; break; }
    if (!nonz) seed.set_bit(0, 1);

    int N_total = 2 * kg;

    NTL::vec_GF2 seq = build_state_seq(seed, 0, N_total);
    NTL::GF2X mu;
    NTL::MinPolySeq(mu, seq, kg);
    chi = gf2x_lcm(chi, mu);

    for (int coord = 1; coord < kg && NTL::deg(chi) < kg; coord += std::max(1, kg / 32)) {
        NTL::vec_GF2 seq2 = build_state_seq(seed, coord, N_total);
        NTL::GF2X mu2;
        NTL::MinPolySeq(mu2, seq2, kg);
        chi = gf2x_lcm(chi, mu2);
    }
    return chi;
}

}  // anonymous namespace

// ─────────────────────────────────────────────────────────────────────
MeLatResult test_me_notprimitive(
    const std::vector<Generator*>& gens,
    const std::vector<std::vector<Transformation*>>& trans,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse)
{
    CombinedView view_proto(gens, trans, kg, L);

    // Stage 1 — recover χ_f.  Strategy: try polychar_comb (= product
    // of each component's gen.char_poly()) first.  For generators with
    // analytical char_poly() overrides (MT, MELGGen, TauswortheGen, …) this
    // gives the EXACT χ_f; select_phi then identifies the largest
    // Mersenne primitive factor and the fast path triggers.
    //
    // For generators with no analytical char_poly() (SFMTGen defaults to
    // BM on bit 0), polychar_comb's result may be missing factors.  In
    // that case fall back to a state-coordinate Krylov BM via
    // recover_char_poly(): this captures factors that bit-0 BM misses
    // (e.g. SFMTGen-607's φ_607 lives in lanes other than u[0] bit 31).
    //
    // The two are then LCM'd so we keep the union of factor information.
    NTL::GF2X chi;
    {
        BitVect chi_bv = polychar_comb(gens);
        NTL::SetCoeff(chi, kg);
        for (int j = 0; j < kg; j++)
            if (chi_bv.get_bit(j)) NTL::SetCoeff(chi, j);
    }
    PhiPick pick = select_phi(chi);

    // If polychar_comb didn't yield a "good" φ (low degree relative to
    // kg, or no Mersenne primitive certified), try state-coord BM and
    // LCM.  This catches SFMTGen-style cases where the analytical char
    // poly path is incomplete.
    if (pick.p < kg / 2) {
        NTL::GF2X chi_bm = recover_char_poly(view_proto, kg);
        chi = gf2x_lcm(chi, chi_bm);
        pick = select_phi(chi);
    }
    int p = pick.p;

    if (p == kg) {
        return test_me_lat(gens, trans, kg, L, maxL, delta, mse);
    }

    BitVect M_phi(p + 1);
    for (int i = 0; i <= p; i++) {
        if (NTL::IsOne(NTL::coeff(pick.phi, i))) M_phi.set_bit(i, 1);
    }

    int Deg = p;
    int RES = std::min(maxL, Deg);
    int Phi = std::max(1, view_proto.output_phases());

    NTL::GF2X chi_psi;
    NTL::div(chi_psi, chi, pick.phi);
    long m_psi = NTL::deg(chi_psi);

    auto random_nonzero = [&](int kbits, std::mt19937_64& rng) {
        BitVect s(kbits);
        for (int i = 0; i < kbits; i++) if (rng() & 1ULL) s.set_bit(i, 1);
        bool any = false;
        for (int w = 0; w < s.nwords(); w++) if (s.data()[w]) { any = true; break; }
        if (!any) s.set_bit(0, 1);
        return s;
    };

    std::mt19937_64 rng(0xD15C0FFEDB17EULL);
    std::vector<BitVect> seeds(gens.size());
    for (size_t j = 0; j < gens.size(); j++) {
        std::vector<Generator*> single = {gens[j]};
        std::vector<std::vector<Transformation*>> single_trans = {trans[j]};
        CombinedView view_j(single, single_trans, gens[j]->k(), gens[j]->L());

        for (int attempt = 0; attempt < 16; attempt++) {
            BitVect s_random = random_nonzero(gens[j]->k(), rng);
            BitVect s_in_V = (m_psi >= 0)
                ? apply_polynomial(chi_psi, s_random, view_j)
                : s_random;
            bool any = false;
            for (int w = 0; w < s_in_V.nwords(); w++)
                if (s_in_V.data()[w]) { any = true; break; }
            if (any) { seeds[j] = std::move(s_in_V); break; }
        }
        if (seeds[j].nbits() == 0) {
            throw std::runtime_error(
                "test_me_notprimitive: V-nonzero seed construction failed");
        }
    }

    MeLatResult res;
    res.ecart.assign(maxL + 1, 0);
    res.se = 0;

    auto read_combined_now = [&](std::vector<std::unique_ptr<Generator>>& gen_copies) {
        BitVect combined(L);
        for (size_t j = 0; j < gens.size(); j++) {
            BitVect out = gen_copies[j]->get_output().copy();
            for (auto* t : trans[j]) t->apply(out);
            int n = std::min(L, out.nbits());
            for (int b = 0; b < n; b++)
                if (out.get_bit(b)) combined.set_bit(b, combined.get_bit(b) ^ 1);
        }
        return combined;
    };

    for (int sm = 0; sm < Phi; sm++) {
        std::vector<std::unique_ptr<Generator>> gen_copies;
        gen_copies.reserve(gens.size());
        for (size_t j = 0; j < gens.size(); j++)
            gen_copies.push_back(gens[j]->copy());
        for (size_t j = 0; j < gens.size(); j++)
            gen_copies[j]->init(seeds[j]);
        for (int s = 0; s < sm; s++) {
            for (size_t j = 0; j < gens.size(); j++) gen_copies[j]->next();
        }

        std::vector<BitVect> A(RES, BitVect(p));
        BitVect combined = read_combined_now(gen_copies);
        for (int i = 0; i < RES; i++) A[i].set_bit(0, combined.get_bit(i));
        for (int k = 1; k < p; k++) {
            for (size_t j = 0; j < gens.size(); j++) gen_copies[j]->next();
            combined = read_combined_now(gen_copies);
            for (int i = 0; i < RES; i++) A[i].set_bit(k, combined.get_bit(i));
        }

        BitVect M_sm = M_phi.copy();
        int Deg_sm = p;
        int RES_sm = std::min(maxL, Deg_sm);

        std::vector<BitVect> polys(RES_sm, BitVect(p + 1));
        BitVect temp(p + 1);
        for (int i = 0; i < RES_sm; i++) {
            polys[i] = BitVect(p + 1);
            for (int k = 0; k < p; k++) {
                if (A[i].get_bit(k)) {
                    temp = M_sm.copy();
                    temp.lshift(k + 1);
                    polys[i].xor_with(temp);
                }
            }
            polys[i].and_mask(p);
        }

        Deg_sm = normalize_polys(polys, M_sm, Deg_sm, RES_sm);
        RES_sm = std::min(maxL, Deg_sm);

        DualLatticeBase base(RES_sm, Deg_sm);
        base.dual_base(polys, M_sm, 1);

        for (int l = 1; l <= RES_sm; l++) {
            int length = base.lenstra(l);
            int d = std::min(length, Deg_sm / l);
            int ec = Deg_sm / l - d;
            if (ec > res.ecart[l]) res.ecart[l] = ec;
            if (l != RES_sm) base.dual_base_increase(polys);
        }
    }

    int maxl_overall = maxL;
    for (int l = 1; l <= maxL; l++) {
        res.se += res.ecart[l];
        if (res.ecart[l] > delta[l] || res.se > mse) {
            maxl_overall = l;
            break;
        }
    }
    res.se = 0;
    for (int l = 1; l <= maxl_overall; l++) res.se += res.ecart[l];
    for (int l = maxl_overall + 1; l <= maxL; l++) res.ecart[l] = INT_MAX;
    return res;
}
