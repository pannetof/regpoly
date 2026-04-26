#pragma once
#include "bitvect.h"
#include "generateur.h"
#include "transformation.h"
#include "lattice_polys.h"  // MeLatResult
#include <vector>

// Equidistribution test that does NOT assume the combined generator is
// full-period (i.e. does not assume the characteristic polynomial χ_f is
// primitive of degree kg).
//
// Pipeline (per docs/C.md):
//   1. Recover χ_f via Berlekamp–Massey on a scalar functional of the
//      combined output; LCM across several initial states / functionals.
//   2. Factor χ_f over F_2 (NTL CanZass on GF2X).
//   3. Pick the largest-degree irreducible factor φ whose period 2^p − 1
//      can be certified maximal (Mersenne fast path; full Knuth check
//      otherwise; fall back to "highest certifiable" if intractable).
//      Set V = Ker φ(f), p = deg φ.
//   4. Build a basis B of V (Route B: pick random s, compute s' = ψ(f)·s
//      via Horner where ψ = χ_f / φ, collect orbit {f^v · s' : v=0..p-1};
//      Route A fallback if rank < p).
//   5. Run the matricial DE core on V: maintain p "virtual register"
//      clones of the combined generator, each initialised to one B[v];
//      step them in lockstep; read L-bit outputs; insert dual rows into
//      a single F_2 echelon (one per output-bit phase) and read off
//      k(v) for v = 1..maxL.
//
// The reported k(v) is on the invariant subspace V, so the corresponding
// equidistribution gap is ecart[v] = floor(p/v) − k(v) (NOT kg/v − k(v));
// the upper bound is p, not kg.  This is the "honest" k(v) for any
// F_2-linear generator regardless of full-period status.

MeLatResult test_me_notprimitive(
    const std::vector<Generateur*>& gens,
    const std::vector<std::vector<Transformation*>>& trans,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse);
