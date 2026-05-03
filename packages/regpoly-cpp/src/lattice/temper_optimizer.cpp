#include "temper_optimizer.h"
#include <climits>
#include <algorithm>

// ── NTL helpers (same as lattice_polys.cpp) ─────────────────────────────

static NTL::GF2X bv_to_gf2x_local(const BitVect& bv, int deg) {
    NTL::GF2X p;
    for (int i = 0; i < deg; i++)
        if (bv.get_bit(i)) NTL::SetCoeff(p, i);
    return p;
}

static void gf2x_to_bv_local(BitVect& bv, const NTL::GF2X& poly) {
    bv.zero();
    for (int i = 0; i <= NTL::deg(poly); i++)
        if (IsOne(coeff(poly, i))) bv.set_bit(i, 1);
}

// Get tempered output from a generator
static BitVect get_transformed_output_local(
    Generator* gen, const std::vector<Transformation*>& chain)
{
    if (chain.empty())
        return gen->get_output();
    BitVect out = gen->get_output();
    for (auto* t : chain) t->apply(out);
    BitVect result(gen->L());
    result.copy_part_from(out, gen->L());
    return result;
}

// ══════════════════════════════════════════════════════════════════════════
// TemperOptCache implementation
// ══════════════════════════════════════════════════════════════════════════

TemperOptCache::TemperOptCache(
    const std::vector<Generator*>& gens,
    const std::vector<std::vector<Transformation*>>& trans,
    int kg, int L)
    : gens_(gens), trans_(trans), kg_orig_(kg), Deg_(kg), L_(L), RES_(std::min(L, kg))
{
    // 1. Combined characteristic polynomial
    M_ = polychar_comb(gens);
    M_orig_ = M_.copy();

    // 2. Find ALL generating polynomials
    find_polys(gens, trans, Deg_, M_, polys_, RES_);

    // 3. Normalize (may reduce Deg_ if g_0 shares factor with M)
    Deg_ = normalize_polys(polys_, M_, Deg_, RES_);
    RES_ = std::min(L_, Deg_);

    // 4. Cache the NTL representations for single-poly recomputation
    M_ntl_ = bv_to_gf2x_local(M_, Deg_ + 1);
    NTL::GF2X g0 = bv_to_gf2x_local(polys_[0], Deg_);
    // polys_[0] is already normalized to 1, but we need inv_g0 from
    // the ORIGINAL g_0 (before normalization). Since polys_ are already
    // normalized, h_0 = 1, so inv_g0 was used to produce them.
    // We can recover inv_g0 from polys_[0] = 1 → inv_g0 was g_0^{-1}.
    // For recomputation of a single poly: h_v = g_v_raw * inv_g0 mod M.
    // We need to store inv_g0. Recompute it from the raw g_0.

    // Recompute raw g_0 to get inv_g0
    std::vector<BitVect> raw_polys;
    find_polys(gens, trans, Deg_, M_, raw_polys, 1);

    // Handle GCD reduction
    NTL::GF2X raw_g0 = bv_to_gf2x_local(raw_polys[0], Deg_);
    NTL::GF2X gcd_val = NTL::GCD(raw_g0, M_ntl_);
    if (NTL::deg(gcd_val) > 0)
        raw_g0 /= gcd_val;
    inv_g0_ = NTL::InvMod(raw_g0, M_ntl_);

    // 5. Allocate stack (empty, populated by compute_all)
    stack_.resize(RES_ + 1);
}

std::vector<int> TemperOptCache::compute_all() {
    std::vector<int> ecart(L_ + 1, -1);

    auto base = std::make_unique<DualLatticeBase>(RES_, Deg_);
    base->dual_base(polys_, M_, 1);

    for (int l = 1; l <= RES_; l++) {
        // Save snapshot BEFORE dual_base_increase + lenstra at l.
        // For l=1: snapshot has 1 row (just P).
        // For l>1: snapshot has l-1 rows (after previous lenstra).
        // compute_gap(l) will restore this, add row l via
        // dual_base_increase, then run lenstra(l).
        stack_[l] = std::make_unique<DualLatticeBase>(*base);

        if (l > 1)
            base->dual_base_increase(polys_);

        int length = base->lenstra(l);
        int d = std::min(length, Deg_ / l);
        ecart[l] = Deg_ / l - d;
    }

    for (int l = 1; l <= L_; l++)
        if (ecart[l] == -1) ecart[l] = 0;

    return ecart;
}

std::vector<int> TemperOptCache::rebuild() {
    // Restore original M and degree
    M_ = M_orig_.copy();
    Deg_ = kg_orig_;
    RES_ = std::min(L_, Deg_);

    // Recompute all generating polynomials from current bitmask
    find_polys(gens_, trans_, Deg_, M_, polys_, RES_);

    // Normalize (may reduce Deg_)
    Deg_ = normalize_polys(polys_, M_, Deg_, RES_);
    RES_ = std::min(L_, Deg_);

    // Recompute NTL representation and inv_g0
    M_ntl_ = bv_to_gf2x_local(M_, Deg_ + 1);
    refresh_inv_g0();

    // Rebuild all snapshots
    stack_.clear();
    stack_.resize(RES_ + 1);
    return compute_all();
}

// ── Incremental optimization (mirrors temperMK.c) ───────────────────

void TemperOptCache::reset_step() {
    previous_v_ = 0;
    current_base_.reset();
    opt_stack_.clear();
    opt_stack_.resize(RES_);  // indices 0..RES_-2 used for v=2..RES_
}

int TemperOptCache::step(int v) {
    if (v < 1 || v > RES_) return 0;

    if (v == 1) {
        // Always rebuild from scratch at v=1 (perturbation may change g_0)
        refresh_inv_g0();

        BitVect new_poly(Deg_ + 1);
        find_single_poly(0, new_poly);
        normalize_single(new_poly);
        polys_[0] = new_poly;

        current_base_ = std::make_unique<DualLatticeBase>(RES_, Deg_);
        current_base_->dual_base(polys_, M_, 1);

    } else if (v > previous_v_) {
        // Forward: save snapshot, then DualBaseIncrease
        BitVect new_poly(Deg_ + 1);
        find_single_poly(v - 1, new_poly);
        normalize_single(new_poly);
        polys_[v - 1] = new_poly;

        opt_stack_[v - 2] = std::make_unique<DualLatticeBase>(*current_base_);
        current_base_->dual_base_increase(polys_);

    } else {
        // Retry (v <= previous_v_): restore snapshot, then DualBaseIncrease
        BitVect new_poly(Deg_ + 1);
        find_single_poly(v - 1, new_poly);
        normalize_single(new_poly);
        polys_[v - 1] = new_poly;

        *current_base_ = *opt_stack_[v - 2];
        current_base_->dual_base_increase(polys_);
    }

    previous_v_ = v;
    int length = current_base_->lenstra(v);
    int d = std::min(length, Deg_ / v);
    return Deg_ / v - d;
}

// ── Single-poly helpers ─────────────────────────────────────────────

void TemperOptCache::find_single_poly(int idx, BitVect& out_poly) {
    // Collect output bit 'idx' from the generator for K steps
    int K = Deg_;
    int J = (int)gens_.size();

    std::vector<std::unique_ptr<Generator>> gen_copies;
    for (int j = 0; j < J; j++)
        gen_copies.push_back(gens_[j]->copy());

    // Init with canonical e_0
    for (int j = 0; j < J; j++) {
        BitVect bc(gens_[j]->k());
        bc.set_bit(0, 1);
        gen_copies[j]->init(bc);
    }

    // Collect one bit per step
    BitVect A(K);
    for (int k = 0; k < K; k++) {
        if (k > 0)
            for (int j = 0; j < J; j++)
                gen_copies[j]->next();

        BitVect combined(gens_[0]->L());
        combined.zero();
        for (int j = 0; j < J; j++) {
            BitVect out = get_transformed_output_local(
                gen_copies[j].get(), trans_[j]);
            combined.xor_with(out);
        }
        A.set_bit(k, combined.get_bit(idx));
    }

    // g_idx = sum_{k: A[k]=1} M >> (k+1), masked to K bits
    out_poly = BitVect(K + 1);
    out_poly.zero();
    BitVect temp(K + 1);
    for (int k = 0; k < K; k++) {
        if (A.get_bit(k)) {
            temp = M_.copy();
            temp.lshift(k + 1);
            out_poly.xor_with(temp);
        }
    }
    out_poly.and_mask(K);
}

void TemperOptCache::normalize_single(BitVect& poly) {
    NTL::GF2X g = bv_to_gf2x_local(poly, Deg_);
    // Handle GCD reduction (same factor as during init)
    NTL::GF2X gcd_val = NTL::GCD(g, M_ntl_);
    if (NTL::deg(gcd_val) > 0)
        g /= gcd_val;
    NTL::GF2X h = NTL::MulMod(g, inv_g0_, M_ntl_);
    gf2x_to_bv_local(poly, h);
}

void TemperOptCache::refresh_inv_g0() {
    // Recompute raw g_0 for the current bitmask
    std::vector<BitVect> raw_polys;
    find_polys(gens_, trans_, Deg_, M_, raw_polys, 1);

    NTL::GF2X raw_g0 = bv_to_gf2x_local(raw_polys[0], Deg_);
    NTL::GF2X gcd_val = NTL::GCD(raw_g0, M_ntl_);
    if (NTL::deg(gcd_val) > 0)
        raw_g0 /= gcd_val;
    inv_g0_ = NTL::InvMod(raw_g0, M_ntl_);
}

int TemperOptCache::compute_gap(int v) {
    if (v < 1 || v > RES_) return 0;
    if (!stack_[v]) return 0;

    // 1. Restore the cached basis.
    //    stack_[v] was saved BEFORE dual_base_increase + lenstra at v.
    //    It has l-1 rows for v > 1, or 1 row (P) for v = 1.
    DualLatticeBase base(*stack_[v]);

    // 2. Recompute the generating polynomial for output bit v-1
    //    using the CURRENT bitmask, normalize it.
    if (v > 1) {
        BitVect new_poly(Deg_ + 1);
        find_single_poly(v - 1, new_poly);
        normalize_single(new_poly);
        polys_[v - 1] = new_poly;

        // Add the new polynomial to the basis (row v-1)
        base.dual_base_increase(polys_);
    }

    // 3. Run Lenstra at resolution v
    int length = base.lenstra(v);
    int d = std::min(length, Deg_ / v);
    return Deg_ / v - d;
}
