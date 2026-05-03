#include "me_helpers.h"
#include "dual_lattice.h"
#include <NTL/GF2X.h>
#include <climits>
#include <algorithm>
#include <stdexcept>

// ══════════════════════════════════════════════════════════════════════════
// NTL conversion helpers
// ══════════════════════════════════════════════════════════════════════════

// BitVect bit i = coefficient of z^i → NTL GF2X
static NTL::GF2X bv_to_gf2x(const BitVect& bv, int deg) {
    NTL::GF2X p;
    for (int i = 0; i < deg; i++)
        if (bv.get_bit(i))
            NTL::SetCoeff(p, i);
    return p;
}

// NTL GF2X → BitVect bit i = coefficient of z^i
static void gf2x_to_bv(BitVect& bv, const NTL::GF2X& poly) {
    bv.zero();
    for (int i = 0; i <= NTL::deg(poly); i++) {
        if (IsOne(coeff(poly, i)))
            bv.set_bit(i, 1);
    }
}

// Convert a single generator's char_poly() to NTL GF2X.
// char_poly() returns K bits: bit j = coefficient of z^j (no leading z^K term).
// The full polynomial is z^K + sum_{j=0}^{K-1} bv[j]*z^j.
static NTL::GF2X char_poly_to_gf2x(Generator* gen) {
    BitVect bv = gen->char_poly();
    int K = gen->k();
    NTL::GF2X p;
    NTL::SetCoeff(p, K);  // leading z^K
    for (int j = 0; j < K; j++)
        if (bv.get_bit(j))
            NTL::SetCoeff(p, j);
    return p;
}

// Get L-bit output from a generator, applying transformation chain.
// Uses get_output() which handles circular buffer rotation (e.g. MT).
static BitVect get_transformed_output(
    Generator* gen,
    const std::vector<Transformation*>& chain)
{
    if (chain.empty())
        return gen->get_output();
    BitVect out = gen->get_output();
    for (auto* t : chain)
        t->apply(out);
    BitVect result(gen->L());
    result.copy_part_from(out, gen->L());
    return result;
}

// ══════════════════════════════════════════════════════════════════════════
// polychar_comb — product of individual characteristic polynomials
// ══════════════════════════════════════════════════════════════════════════

BitVect polychar_comb(const std::vector<Generator*>& gens) {
    NTL::GF2X product;
    NTL::SetCoeff(product, 0);  // start with 1

    for (auto* gen : gens) {
        NTL::GF2X pi = char_poly_to_gf2x(gen);
        product *= pi;
    }

    // Total degree = sum of individual k values
    int K_total = NTL::deg(product);
    BitVect M(K_total + 1);
    gf2x_to_bv(M, product);
    return M;
}

// ══════════════════════════════════════════════════════════════════════════
// find_polys — build generating polynomials g_i(z) for each output bit
// ══════════════════════════════════════════════════════════════════════════

void find_polys(
    const std::vector<Generator*>& gens,
    const std::vector<std::vector<Transformation*>>& trans,
    int K, const BitVect& M,
    std::vector<BitVect>& polys, int resolution)
{
    int J = (int)gens.size();

    // Create working copies of all generators
    std::vector<std::unique_ptr<Generator>> gen_copies;
    gen_copies.reserve(J);
    for (int j = 0; j < J; j++)
        gen_copies.push_back(gens[j]->copy());

    // Build matrix A[resolution][K]: A[i][k] = bit i of combined output at step k
    // Store as resolution BitVects of K bits each
    std::vector<BitVect> A(resolution, BitVect(K));

    // Step 0: initialize each component with canonical vector e_0
    for (int j = 0; j < J; j++) {
        BitVect bc(gens[j]->k());
        bc.set_bit(0, 1);
        gen_copies[j]->init(bc);
    }

    // Collect combined output at step 0
    {
        BitVect combined(gens[0]->L());
        combined.zero();
        for (int j = 0; j < J; j++) {
            BitVect out = get_transformed_output(gen_copies[j].get(), trans[j]);
            combined.xor_with(out);
        }
        for (int i = 0; i < resolution; i++)
            A[i].set_bit(0, combined.get_bit(i));
    }

    // Steps 1..K-1: iterate all generators, collect combined output
    for (int k = 1; k < K; k++) {
        for (int j = 0; j < J; j++)
            gen_copies[j]->next();

        BitVect combined(gens[0]->L());
        combined.zero();
        for (int j = 0; j < J; j++) {
            BitVect out = get_transformed_output(gen_copies[j].get(), trans[j]);
            combined.xor_with(out);
        }
        for (int i = 0; i < resolution; i++)
            A[i].set_bit(k, combined.get_bit(i));
    }

    // Compute polynomials: g_i(z) = sum_{k: A[i][k]=1} z^{k+1} * M(z), masked to K bits
    polys.resize(resolution);
    BitVect temp(K + 1);
    for (int i = 0; i < resolution; i++) {
        polys[i] = BitVect(K + 1);
        polys[i].zero();
        for (int k = 0; k < K; k++) {
            if (A[i].get_bit(k)) {
                // temp = M << (k+1) (multiply M by z^{k+1})
                temp = M.copy();
                temp.lshift(k + 1);
                polys[i].xor_with(temp);
            }
        }
        // Mask to K bits (keep only coefficients of z^0 through z^{K-1})
        polys[i].and_mask(K);
    }
}

// ══════════════════════════════════════════════════════════════════════════
// normalize_polys — multiply each g_i by g_1^{-1} mod M(z) using NTL
// ══════════════════════════════════════════════════════════════════════════

int normalize_polys(
    std::vector<BitVect>& polys, BitVect& M, int K, int resolution)
{
    NTL::GF2X M_ntl = bv_to_gf2x(M, K + 1);
    NTL::GF2X g1_ntl = bv_to_gf2x(polys[0], K);

    // If g_0 shares a common factor with M, divide it out from all polynomials
    // and from M before normalizing. This handles generators where the output
    // function produces a degenerate sequence from the canonical initial state.
    NTL::GF2X g = NTL::GCD(g1_ntl, M_ntl);
    if (NTL::deg(g) > 0) {
        M_ntl /= g;
        for (int i = 0; i < resolution; i++) {
            NTL::GF2X gi = bv_to_gf2x(polys[i], K);
            gi /= g;
            gf2x_to_bv(polys[i], gi);
        }
        g1_ntl /= g;
        K = NTL::deg(M_ntl);
        // Update M to the reduced polynomial
        M = BitVect(K + 1);
        gf2x_to_bv(M, M_ntl);
    }

    NTL::GF2X inv_g1 = NTL::InvMod(g1_ntl, M_ntl);

    for (int i = 0; i < resolution; i++) {
        NTL::GF2X gi = bv_to_gf2x(polys[i], K);
        NTL::GF2X hi = NTL::MulMod(gi, inv_g1, M_ntl);
        gf2x_to_bv(polys[i], hi);
    }
    return K;
}

// ══════════════════════════════════════════════════════════════════════════
// test_me_lat — full TestMELat algorithm
// ══════════════════════════════════════════════════════════════════════════

MeLatResult test_me_lat(
    const std::vector<Generator*>& gens,
    const std::vector<std::vector<Transformation*>>& trans,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse)
{
    int Deg = kg;
    int RES = std::min(maxL, Deg);

    MeLatResult result;
    result.ecart.resize(maxL + 1, -1);
    result.se = 0;

    // 1. Compute combined characteristic polynomial
    BitVect M = polychar_comb(gens);

    // 2-3. Find and normalize generating polynomials
    // normalize_polys may reduce Deg if g_0 shares a factor with M
    std::vector<BitVect> polys;
    find_polys(gens, trans, Deg, M, polys, RES);
    Deg = normalize_polys(polys, M, Deg, RES);
    RES = std::min(maxL, Deg);

    // 4. Build dual lattice basis and run Lenstra
    DualLatticeBase base(RES, Deg);
    base.dual_base(polys, M, 1);

    int maxl = maxL;
    for (int l = 1; l <= RES; l++) {
        int length = base.lenstra(l);
        int d = std::min(length, Deg / l);
        result.ecart[l] = Deg / l - d;
        result.se += result.ecart[l];

        if (result.ecart[l] > delta[l] || result.se > mse) {
            maxl = l;
            break;
        }
        if (l != RES)
            base.dual_base_increase(polys);
    }

    // Finalize: set remaining ecarts
    result.se = 0;
    for (int l = 1; l <= maxl; l++) {
        if (result.ecart[l] == -1)
            result.ecart[l] = 0;
        result.se += result.ecart[l];
    }
    for (int l = maxl + 1; l <= maxL; l++) {
        if (result.ecart[l] == -1)
            result.ecart[l] = INT_MAX;
    }

    return result;
}
