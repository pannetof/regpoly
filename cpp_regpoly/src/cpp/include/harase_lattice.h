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
