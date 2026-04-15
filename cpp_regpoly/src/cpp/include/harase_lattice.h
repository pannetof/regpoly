#pragma once
#include "bitvect.h"
#include "generateur.h"
#include "transformation.h"
#include "lattice_polys.h"  // MeLatResult
#include <vector>

// Harase-Matsumoto-Saito (2011) fast lattice reduction for equidistribution.
//
// Works with the PRIMAL lattice (not dual).  Uses Mulders-Storjohann
// weak reduction: simpler inner loop than Lenstra, no linear system solve.
// Avoids the polynomial inversion step (no NTL InvMod needed for the
// lattice reduction itself).
//
// The generating polynomials are still computed the same way as in the
// Couture-L'Ecuyer method (run the generator from canonical state for k
// steps, collect output bits, form g_i(z)).
//
// Reference:
//   Harase, Matsumoto, Saito (2011). "Fast Lattice Reduction for
//   F2-Linear Pseudorandom Number Generators." Math. Comp. 80(273).

MeLatResult test_me_lat_harase(
    const std::vector<Generateur*>& gens,
    const std::vector<std::vector<Transformation*>>& trans,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse);

// Compute k(v) for a SINGLE resolution v using the PIS method.
// Returns the equidistribution dimension (not the gap).
// Cost: O(k/v) generator steps — much cheaper than full test_me_lat_harase.
int compute_kv(
    const std::vector<Generateur*>& gens,
    const std::vector<std::vector<Transformation*>>& trans,
    int kg, int v);

// ── PIS basis with StackBase caching for tempering optimization ─────────
// Computes all k(v) from v=L down to 1, caching the basis at each v.
// After a bitmask perturbation, restore_and_reduce(v) restores the
// cached basis at v and re-reduces to get the new k(v).  This avoids
// recomputing resolutions v+1..L.

class PISCache {
public:
    PISCache(const std::vector<Generateur*>& gens,
             const std::vector<std::vector<Transformation*>>& trans,
             int kg, int L);

    // Compute all k(v) for v = L..1.  Populates the cache.
    // Returns ecart[v] = kg/v - k(v) for v = 1..L.
    std::vector<int> compute_all();

    // After a bitmask perturbation: restore cached basis at resolution v,
    // rebuild the generator vector, reduce, return the NEW k(v).
    // Resolutions v+1..L are unchanged (guaranteed by safe masks).
    // Cost: O(k/v) — same as a single PIS reduction.
    int restore_and_reduce(int v);

    int kg() const { return kg_; }
    int L() const { return L_; }

private:
    std::vector<Generateur*> gens_;
    std::vector<std::vector<Transformation*>> trans_;
    int kg_;
    int L_;

    // The implementation details are in the .cpp file (opaque structs).
    struct Impl;
    std::shared_ptr<Impl> impl_;
};
