#pragma once
#include "bitvect.h"
#include "generateur.h"
#include "transformation.h"
#include <vector>
#include <memory>

struct MeLatResult {
    std::vector<int> ecart;  // indexed 0..maxL, ecart[0] unused
    int se;
};

// Compute the product of the individual characteristic polynomials.
// Returns a BitVect of (K_total+1) bits where bit i = coefficient of z^i.
BitVect polychar_comb(const std::vector<Generateur*>& gens);

// Build the raw generating polynomials g_i(z) for i=0..resolution-1.
// Each polynomial is stored in polys[i] as a BitVect of (K+1) bits,
// bit j = coefficient of z^j.
void find_polys(
    const std::vector<Generateur*>& gens,
    const std::vector<std::vector<Transformation*>>& trans,
    int K, const BitVect& M,
    std::vector<BitVect>& polys, int resolution);

// Normalize the polynomials: multiply each by g_1^{-1} mod M(z).
// If g_0 shares a common factor with M, divides it out from M and all polys.
// Updates M and returns the effective degree (may be less than K).
int normalize_polys(
    std::vector<BitVect>& polys, BitVect& M, int K, int resolution);

// Full TestMELat: compute ecart[l] for l=1..min(maxL, kg).
MeLatResult test_me_lat(
    const std::vector<Generateur*>& gens,
    const std::vector<std::vector<Transformation*>>& trans,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse);
