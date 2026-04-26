#include "simd_de.h"
#include "lattice_polys.h"   // test_me_lat (full-period fast path) +
                             // polychar_comb (chi recovery)
#include <algorithm>
#include <climits>
#include <memory>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/GF2.h>
#include <NTL/vec_GF2.h>

namespace {

// ─────────────────────────────────────────────────────────────────────
// Helpers shared with notprimitive_de.cpp (chi recovery, V-projection).
// Inlined here to keep this file self-contained.
// ─────────────────────────────────────────────────────────────────────

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

NTL::GF2X gf2x_lcm(const NTL::GF2X& a, const NTL::GF2X& b) {
    if (NTL::IsZero(a)) return b;
    if (NTL::IsZero(b)) return a;
    NTL::GF2X g, ab, q;
    NTL::GCD(g, a, b);
    NTL::mul(ab, a, b);
    NTL::div(q, ab, g);
    return q;
}

struct PhiPick {
    NTL::GF2X phi;
    int p;
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
        if (d > 0 && is_mersenne_prime_exponent((int)d) &&
            !NTL::IsZero(NTL::coeff(fac[i].a, 0)))
            return { fac[i].a, (int)d };
    }
    for (long i : order) {
        long d = NTL::deg(fac[i].a);
        if (d > 0 && !NTL::IsZero(NTL::coeff(fac[i].a, 0)))
            return { fac[i].a, (int)d };
    }
    throw std::runtime_error(
        "test_me_simd_notprimitive: chi has no nontrivial irreducible factor");
}

// step_once / apply_polynomial used to build V-projected seeds.
//
// step_once via init+next applies one full f-step.  CRITICAL: SFMT's
// init() runs rebuild_work_from_state which has an all-zero guard
// (inserts a synthetic 1-bit if state is zero) — that breaks
// F2-linearity for apply_polynomial when intermediate r = 0.  We use
// set_raw_state + simd_advance_one_word * N instead, which faithfully
// replicates one full step without the guard.  pword_idx_ is set to
// N-1 before the N-cycle so the advance order matches generate_all
// (state slots 0, 1, ..., N-1 in that order, matching the matrix
// recovered by recover_chi_state_coord_bm via init+next).
BitVect step_once(const Generateur& proto, const BitVect& r) {
    auto g = proto.copy();
    g->set_raw_state(r);
    // Use the generator's declared full-step word count so dSFMT (which
    // carries an extra "lung" word in its F2-linear state) advances by
    // N words per step rather than N+1.  SFMT and any future
    // lung-less SIMD gen has simd_full_step_words() == k/128 = N.
    int n_per_word = g->simd_full_step_words();
    if (n_per_word <= 0) {
        // Non-SIMD fallback: one next() = one full f-step.
        g->simd_reset_word_index();
        g->next();
        return g->state().copy();
    }
    // SFMT path: replicate generate_all via N per-word advances,
    // starting at pword_idx_ = N-1 so the first advance updates state[0]
    // (matching generate_all's update order [0, 1, ..., N-1]).
    g->simd_reset_word_index();   // sets pword_idx_ = 0
    // Advance N-1 extra times so first useful advance starts at slot 0.
    // Cleaner: just advance N times — for SFMT, simd_advance_one_word
    // increments-then-updates, so starting at pword_idx_=0 we update
    // slots [1, 2, ..., N-1, 0] in N advances; that's equivalent to
    // generate_all up to a cyclic shift, which doesn't affect linearity
    // (M_full is the same composition either way modulo rotation).
    for (int i = 0; i < n_per_word; i++) {
        g->simd_advance_one_word();
    }
    return g->state().copy();
}

BitVect apply_polynomial(const NTL::GF2X& g, const BitVect& s,
                         const Generateur& proto)
{
    long d = NTL::deg(g);
    BitVect r(s.nbits());
    for (long i = d; i >= 0; i--) {
        if (i < d) r = step_once(proto, r);
        if (NTL::IsOne(NTL::coeff(g, i))) r.xor_with(s);
    }
    return r;
}


// ─────────────────────────────────────────────────────────────────────
// SimdLinVec — port of MTToolBox simd_linear_generator_vector
// (AlgorithmSIMDEquidistribution.hpp:86–222).
//
// Represents one PIS basis vector: the underlying generator clone
// (carrying its full F2-linear state), the lane-interleaved super-word
// `next` (PIS leading vector), the cached previous super-word for
// weight-mode blending, and the polynomial-degree counter `count`.
// ─────────────────────────────────────────────────────────────────────
struct SimdLinVec {
    std::unique_ptr<Generateur> gen;
    int count;
    bool zero;
    BitVect next;        // 128-bit lane-interleaved super-word
    BitVect previous;    // 128-bit cache for weight-mode blending
    int lane_count;
    int L_bits;
    int start_mode;
    int weight_mode;

    // "Rand" basis-vector constructor: clones g (state and pword_idx_
    // preserved); previous, next initialized to zero; count = 0.
    SimdLinVec(const Generateur& g, int lane_count_, int L_bits_,
               int sm, int wm)
        : gen(g.copy()),
          count(0), zero(false),
          next(128), previous(128),
          lane_count(lane_count_), L_bits(L_bits_),
          start_mode(sm), weight_mode(wm)
    {}

    // Standard-basis constructor: gen state cleared to zero,
    // pword_idx_ reset; previous = 0, next = 1 at MSB-first position
    // bit_pos.
    SimdLinVec(const Generateur& g, int bit_pos,
               int lane_count_, int L_bits_, int sm, int wm)
        : gen(g.copy()),
          count(0), zero(false),
          next(128), previous(128),
          lane_count(lane_count_), L_bits(L_bits_),
          start_mode(sm), weight_mode(wm)
    {
        gen->set_raw_state(BitVect(gen->k()));   // zero the gen state
        gen->simd_reset_word_index();
        next.set_bit(bit_pos, 1);
    }

    // Find largest set BitVect bit index in `next` within [0, bit_size).
    // Matches MTToolBox calc_1pos which returns the MSB-position of the
    // trailing (LSB-most) set 1 bit.  In our MSB-first BitVect encoding,
    // bit 0 is the highest, so the LSB-most set bit is the largest index.
    int trailing_bit_pos(int bit_size) const {
        // Words of the 128-bit BitVect: data[0] = bits [0, 64),
        // data[1] = bits [64, 128).  In each word, BitVect bit i within
        // the word maps to LSB-pos (63 - i_in_word) of that word.
        // Largest set BitVect bit = lowest LSB-pos set bit in the
        // highest-index word that has any bits in [0, bit_size).
        if (bit_size <= 0) return -1;
        // Word 1 (bits [64, 128)) — only consider if bit_size > 64.
        if (bit_size > 64) {
            uint64_t v = next.data()[1];
            int unused_lsb = 128 - bit_size;     // lowest LSB bits beyond bit_size
            if (unused_lsb > 0) {
                v &= ~((uint64_t(1) << unused_lsb) - 1);
            }
            if (v) {
                int ctz = __builtin_ctzll(v);    // lowest LSB-pos set bit
                return 127 - ctz;                // BitVect index
            }
        }
        // Word 0 (bits [0, 64)).
        uint64_t v0 = next.data()[0];
        int eff = std::min(64, bit_size);
        int unused_lsb = 64 - eff;
        if (unused_lsb > 0) {
            v0 &= ~((uint64_t(1) << unused_lsb) - 1);
        }
        if (v0) {
            int ctz = __builtin_ctzll(v0);
            return 63 - ctz;
        }
        return -1;
    }

    static bool is_zero_bitvect(const BitVect& bv) {
        for (int w = 0; w < bv.nwords(); w++)
            if (bv.data()[w]) return false;
        return true;
    }

    // GF(2) addition: XOR generator state (index-aligned via SFMT
    // override), XOR next, XOR previous cache.  All three are part of
    // the augmented linear state.
    void add(const SimdLinVec& other) {
        gen->simd_add_state(*other.gen);
        next.xor_with(other.next);
        previous.xor_with(other.previous);
    }

    // Read fresh super-word from gen, blend with `previous` per
    // weight_mode (lanes [0, wm) from fresh, [wm, lane_count) from
    // previous), update previous = fresh, then pack the v MSBs of each
    // blended lane into `next` at MSB-first positions [i*v, (i+1)*v).
    void get_next(int v) {
        BitVect fresh = gen->simd_read_super_word(start_mode);
        BitVect blended(128);
        for (int i = 0; i < lane_count; i++) {
            const BitVect& src = (i < weight_mode) ? fresh : previous;
            int base = i * L_bits;
            for (int b = 0; b < L_bits; b++) {
                if (src.get_bit(base + b)) blended.set_bit(base + b, 1);
            }
        }
        previous = fresh.copy();
        next = BitVect(128);
        for (int i = 0; i < lane_count; i++) {
            int src_base = i * L_bits;
            int dst_base = i * v;
            for (int b = 0; b < v; b++) {
                if (blended.get_bit(src_base + b))
                    next.set_bit(dst_base + b, 1);
            }
        }
    }

    // Advance one f-application step then build a fresh `next`.  If
    // `next` lands at zero, keep advancing (counts as additional
    // polynomial-degree increments) until either a non-zero `next` is
    // produced or 2*k consecutive zero outputs trip the dead-state
    // safeguard (mirrors MTToolBox's bitSize()*2 check).
    void next_state(int v) {
        if (zero) return;
        gen->simd_advance_one_word();
        get_next(v);
        count++;
        int zero_count = 0;
        while (is_zero_bitvect(next)) {
            zero_count++;
            if (zero_count > 2 * gen->k()) {
                zero = true;
                break;
            }
            gen->simd_advance_one_word();
            get_next(v);
            count++;
        }
    }
};

// ─────────────────────────────────────────────────────────────────────
// SIMD-PIS reduction — port of
// AlgorithmSIMDEquidistribution::get_equidist_main (lines 587–729).
//
// `warmed_gen` and `warmed_previous` carry the post-warmup state of the
// rand basis vector (one simd_advance_one_word + simd_read_super_word
// has been done), matching MTToolBox's `work.generate()` warmup before
// constructing AlgorithmSIMDEquidistribution.
// ─────────────────────────────────────────────────────────────────────
int simd_pis_get_equidist(
    const Generateur& warmed_gen, const BitVect& warmed_previous,
    int v, int lane_count, int L_bits, int sm, int wm)
{
    int bit_size = v * lane_count;
    if (bit_size <= 0 || bit_size > 128) return 0;

    std::vector<std::unique_ptr<SimdLinVec>> basis;
    basis.reserve(bit_size + 1);

    // basis[0..bit_size-1]: standard basis.  Each clones warmed_gen but
    // immediately zeros the generator state and resets pword_idx_, so
    // the gen carries no "rand info" yet.  The `next` super-word holds
    // a single 1 at MSB-first position i.
    for (int i = 0; i < bit_size; i++) {
        basis.emplace_back(std::make_unique<SimdLinVec>(
            warmed_gen, i, lane_count, L_bits, sm, wm));
    }

    // basis[bit_size]: rand basis.  Clone warmed_gen (preserving
    // gen state and pword_idx_), seed previous cache from warmup,
    // then advance one step via next_state to populate `next`.
    auto rand_lv = std::make_unique<SimdLinVec>(
        warmed_gen, lane_count, L_bits, sm, wm);
    rand_lv->previous = warmed_previous.copy();
    rand_lv->next_state(v);
    basis.emplace_back(std::move(rand_lv));

    int pivot_index = basis[bit_size]->trailing_bit_pos(bit_size);

    while (!basis[bit_size]->zero) {
        if (pivot_index < 0 || pivot_index >= bit_size) {
            // Defensive: should not happen given the algorithm's
            // monotonic-pivot invariant.  Bail out cleanly.
            break;
        }
        // Keep counts uniform: swap so basis[bit_size] has the lower
        // count (= more "room" to advance in subsequent next_state).
        if (basis[bit_size]->count > basis[pivot_index]->count) {
            std::swap(basis[bit_size], basis[pivot_index]);
        }
        basis[bit_size]->add(*basis[pivot_index]);
        if (SimdLinVec::is_zero_bitvect(basis[bit_size]->next)) {
            basis[bit_size]->next_state(v);
        }
        pivot_index = basis[bit_size]->trailing_bit_pos(bit_size);
    }

    int min_count = INT_MAX;
    for (int i = 0; i < bit_size; i++) {
        if (basis[i]->zero) continue;
        if (basis[i]->count < min_count) min_count = basis[i]->count;
    }
    if (min_count == INT_MAX) return 0;
    return min_count;
}

// ─────────────────────────────────────────────────────────────────────
// Chi recovery via Krylov BM with random linear functionals.
//
// State-coordinate observations (bit i of state) often only "see" the
// primitive component of SFMT's char poly — the non-V dim-31 subspace
// may project trivially onto specific state bits, so LCM of state-coord
// min polys = primitive_factor (degree mexp), missing the degree-31
// quotient factor needed for V-projection.
//
// Random linear functionals (random XOR-combinations of state bits) hit
// both the primitive AND the non-primitive subspaces of the char poly's
// Jordan decomposition, so their min polys equal the FULL min poly of
// the transition matrix (degree = state_dim with high probability).
//
// LCM'd with polychar_comb's analytical χ to ensure the primitive
// factor is found even for tricky configurations.
// ─────────────────────────────────────────────────────────────────────
NTL::GF2X recover_chi_state_coord_bm(
    const std::vector<Generateur*>& gens, int kg, int /*L*/)
{
    // Use the SAME step as apply_polynomial (= step_once = guard-free
    // full step via N per-word advances).  Consistency is critical:
    // chi must annihilate the operator that V-projection uses.
    auto step = [&](const BitVect& r) {
        return step_once(*gens[0], r);
    };

    auto build_func_seq = [&](const BitVect& seed, const BitVect& func, int N_total) {
        NTL::vec_GF2 seq;
        seq.SetLength(N_total);
        BitVect cur = seed.copy();
        for (int N = 0; N < N_total; N++) {
            cur = step(cur);
            uint64_t parity = 0;
            int nw = std::min(func.nwords(), cur.nwords());
            for (int w = 0; w < nw; w++) {
                parity ^= func.data()[w] & cur.data()[w];
            }
            if (__builtin_popcountll(parity) & 1) seq.put(N, 1);
        }
        return seq;
    };

    auto random_bv = [&](int kbits, std::mt19937_64& rng, bool ensure_nonzero) {
        BitVect s(kbits);
        for (int i = 0; i < kbits; i++) if (rng() & 1ULL) s.set_bit(i, 1);
        if (ensure_nonzero) {
            bool any = false;
            for (int w = 0; w < s.nwords(); w++)
                if (s.data()[w]) { any = true; break; }
            if (!any) s.set_bit(0, 1);
        }
        return s;
    };

    NTL::GF2X chi;
    NTL::set(chi);
    std::mt19937_64 rng(0xA110CAB1EE1ULL);
    BitVect seed = random_bv(kg, rng, true);
    int N_total = 2 * kg;

    // Try up to 32 random functionals; stop early if chi reaches degree kg.
    for (int trial = 0; trial < 32 && NTL::deg(chi) < kg; trial++) {
        BitVect func = random_bv(kg, rng, true);
        NTL::vec_GF2 seq = build_func_seq(seed, func, N_total);
        NTL::GF2X mu;
        NTL::MinPolySeq(mu, seq, kg);
        chi = gf2x_lcm(chi, mu);
    }
    return chi;
}

}  // anonymous namespace

// ─────────────────────────────────────────────────────────────────────
// Public entry point: test_me_simd_notprimitive.
// Stages 1–3: chi recovery, factor, select_phi (= primitive component).
// Fast path: full-period non-SIMD generators delegate to test_me_lat.
// Slow path: per (sm, wm) warmup + per-v PIS reduction; aggregate
// per-sm MAX over wm, then MIN over sm; rescale per Saito–Matsumoto
// 2008 formula  e = e * lane_count - (lane_count - wm).
// ─────────────────────────────────────────────────────────────────────
MeLatResult test_me_simd_notprimitive(
    const std::vector<Generateur*>& gens,
    const std::vector<std::vector<Transformation*>>& trans,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse)
{
    int lane_count = 1;
    for (auto* g : gens) {
        int n = g->simd_lane_count();
        if (n > lane_count) lane_count = n;
    }

    // Fast path: pure non-SIMD generator (lane_count == 1).  Reuse the
    // existing notprimitive chi recovery + dual-lattice pipeline via
    // test_me_lat for the full-period case (gives identical k(v) to
    // METHOD_DUALLATTICE).  This collapses the SIMD method to the
    // standard lattice method for non-SIMD generators.
    if (lane_count == 1) {
        BitVect chi_bv = polychar_comb(gens);
        NTL::GF2X chi;
        NTL::SetCoeff(chi, kg);
        for (int j = 0; j < kg; j++)
            if (chi_bv.get_bit(j)) NTL::SetCoeff(chi, j);
        PhiPick pick = select_phi(chi);
        if (pick.p == kg) {
            return test_me_lat(gens, trans, kg, L, maxL, delta, mse);
        }
    }

    // SIMD path: single-component only.
    if (gens.size() != 1) {
        throw std::runtime_error(
            "test_me_simd_notprimitive: multi-component SIMD configurations "
            "are not yet supported.  Use single-component SFMT/dSFMT/MTGP.");
    }

    // V-PROJECTION: per MTToolBox `anni()`, project the seed into
    // V = Ker(phi(f_full)) where phi is the largest mersenne prime
    // irreducible factor of chi_full (the char poly under the FULL
    // step f_full = generate_all = N consecutive per-word advances).
    //
    // V_full is invariant under per-word advance M too: M and M^N
    // commute, so for x in Ker(phi(M^N)), phi(M^N)·M·x = M·phi(M^N)·x = 0.
    // Hence M·x is also in Ker(phi(M^N)) = V_full.  This means
    // V-projecting via the full step is sufficient for the SIMD-PIS,
    // which advances by per-word M.
    //
    // BM via per-word M (simd_advance_one_word) does NOT work for chi
    // recovery — the per-word recurrence is non-stationary (operator
    // depends on pword_idx_), so the state-coord sequence at a single
    // coord is constant within blocks of N steps, giving BM a degenerate
    // sequence and a tiny min poly.  The full step IS stationary.
    // Chi recovery: polychar_comb's "monic char_poly" form (deg=kg) +
    // state-coord random-functional BM (LCM'd if polychar_comb misses
    // the primitive factor, as it does for SFMT).  The combination
    // matches MTToolBox's chi for the full-step operator (= step_once).
    // V_full = Ker(primitive(M_full)) is invariant under per-word M_pw,
    // so V_full-projecting the seed via apply_polynomial(quotient, s, proto)
    // keeps PIS basis vectors in V_full.
    BitVect chi_bv = polychar_comb(gens);
    NTL::GF2X chi;
    NTL::SetCoeff(chi, kg);
    for (int j = 0; j < kg; j++)
        if (chi_bv.get_bit(j)) NTL::SetCoeff(chi, j);
    PhiPick pick = select_phi(chi);
    int target_p = gens[0]->period_exponent();
    if (pick.p < target_p) {
        NTL::GF2X chi_bm = recover_chi_state_coord_bm(gens, kg, L);
        chi = gf2x_lcm(chi, chi_bm);
        pick = select_phi(chi);
    }
    int p = pick.p;
    (void)target_p;

    // V-projection.  In MTToolBox this is `anni(sf)`: apply quotient(f)·s
    // where quotient = chi/primitive.  For SFMT, the per-word transition
    // matrix M_pw has min poly = primitive (degree mexp = p); the char
    // poly factors as primitive * (degree-(kg-p) factor) where (kg-p) =
    // 31 for SFMT-19937, 33 for SFMT-607.  Computing the explicit
    // (kg-p)-degree factor requires the full char poly which our BM
    // can't recover.
    //
    // Workaround: apply M_pw^(kg-p) to a random seed.  Since chi_pw
    // factors as primitive * (something of degree kg-p), and the
    // "something" annihilates a (kg-p)-dimensional subspace W, we have
    // M_pw^(kg-p) maps state into V_pw = Ker(primitive(M_pw)) — the
    // primitive-orbit subspace of dimension p — provided W is killed
    // by M_pw^(kg-p) (which holds when W is the Jordan/invariant block
    // for the non-primitive factor; for SFMT this is typical: the
    // 31-dim non-V subspace is non-invariant under per-word advance and
    // collapses after ≤ kg-p applications of M_pw).
    // V-project via apply_polynomial (init+next = full step).
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
    std::mt19937_64 rng(0x51DABEEF42ULL);
    BitVect seed_v;
    for (int attempt = 0; attempt < 16; attempt++) {
        BitVect s_random = random_nonzero(gens[0]->k(), rng);
        BitVect s_in_V = (m_psi >= 0)
            ? apply_polynomial(chi_psi, s_random, *gens[0])
            : s_random;
        bool any = false;
        for (int w = 0; w < s_in_V.nwords(); w++)
            if (s_in_V.data()[w]) { any = true; break; }
        if (any) { seed_v = std::move(s_in_V); break; }
    }
    if (seed_v.nbits() == 0) {
        throw std::runtime_error(
            "test_me_simd_notprimitive: V-nonzero seed construction failed");
    }

    int RES = std::min(maxL, p);
    std::vector<int> veq(RES + 1, INT_MAX);

    for (int sm = 0; sm < lane_count; sm++) {
        std::vector<int> veq_per_sm(RES + 1, -1);
        for (int wm = 1; wm <= lane_count; wm++) {
            // Build the warmup template per (sm, wm).  Use set_raw_state
            // to avoid init's spurious generate_all (which would advance
            // the state by N words and offset pword_idx_).  After
            // simd_reset_word_index, pword_idx_ = 0; one warmup advance
            // brings it to 1 (matching MTToolBox's `work.generate()`
            // post-setWeightMode).
            auto warmed_gen = gens[0]->copy();
            warmed_gen->set_raw_state(seed_v);
            warmed_gen->simd_reset_word_index();
            warmed_gen->simd_advance_one_word();
            BitVect warmed_previous = warmed_gen->simd_read_super_word(sm);

            for (int v = 1; v <= RES; v++) {
                int e_simd = simd_pis_get_equidist(
                    *warmed_gen, warmed_previous,
                    v, lane_count, L, sm, wm);
                int e = e_simd * lane_count - (lane_count - wm);
                if (e < 0) e = 0;
                if (e > p / v) e = p / v;
                if (e > veq_per_sm[v]) veq_per_sm[v] = e;
            }
        }
        for (int v = 1; v <= RES; v++) {
            if (veq_per_sm[v] >= 0 && veq_per_sm[v] < veq[v])
                veq[v] = veq_per_sm[v];
        }
    }

    MeLatResult res;
    res.ecart.assign(maxL + 1, 0);
    res.se = 0;
    for (int v = 1; v <= RES; v++) {
        int upper = p / v;
        int kv = (veq[v] == INT_MAX) ? 0 : veq[v];
        if (kv > upper) kv = upper;
        res.ecart[v] = upper - kv;
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
