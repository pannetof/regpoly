#pragma once
#include "bitvect.h"
#include "generateur.h"
#include "transformation.h"
#include "poly_lattice.h"
#include "lattice_polys.h"
#include <NTL/GF2X.h>
#include <vector>
#include <memory>

// ══════════════════════════════════════════════════════════════════════════
// LatticeOptCache — dual lattice with StackBase for tempering optimization.
//
// Mirrors the C code's strategy in temperMK.c:
//   - compute_all(): run Lenstra for v = 1..L, caching the basis at each v
//   - compute_gap(v): recompute ONLY resolution v after a bitmask change:
//       1. Restore the cached basis at v-1
//       2. Recompute generating polynomial g_v for the current bitmask
//       3. Normalize: h_v = g_v * g0_inv mod M
//       4. DualBaseIncrease + Lenstra at resolution v
//       5. Return the gap
//
// Cost per compute_gap(v): O(k) for FindPolys at one resolution + one Lenstra.
// ══════════════════════════════════════════════════════════════════════════

class LatticeOptCache {
public:
    LatticeOptCache(
        const std::vector<Generateur*>& gens,
        const std::vector<std::vector<Transformation*>>& trans,
        int kg, int L);

    // Compute all gaps v = 1..L.  Populates the StackBase cache.
    // Returns ecart[v] for v = 1..L (ecart[0] unused).
    std::vector<int> compute_all();

    // After a bitmask perturbation: recompute gap at resolution v only.
    // Restores cached basis at v-1, recomputes g_v for current bitmask,
    // runs Lenstra at v.  Cost: O(k) for sequence collection + one Lenstra.
    int compute_gap(int v);

    // Recompute g_0 and inv_g0 from the current bitmask.
    // Must be called after any perturbation that changes output bit 0.
    void refresh_inv_g0();

    // Recompute all polys, inv_g0, and snapshots from current bitmask.
    std::vector<int> rebuild();

    // ── Incremental optimization (mirrors temperMK.c OptimizeTemper) ──
    //
    // step(v) computes gap at resolution v with incremental StackBase:
    //   v == 1:          rebuild inv_g0 + DualBase from scratch
    //   v > previous_v:  save snapshot, DualBaseIncrease, Lenstra (forward)
    //   v <= previous_v: restore snapshot, DualBaseIncrease, Lenstra (retry)
    //
    // Call reset_step() once before starting the optimization loop.
    void reset_step();
    int step(int v);

    int kg() const { return Deg_; }
    int L() const { return L_; }

private:
    std::vector<Generateur*> gens_;
    std::vector<std::vector<Transformation*>> trans_;
    int kg_orig_;  // original degree (before GCD reduction)
    int Deg_;      // effective degree (after GCD reduction)
    int L_;
    int RES_;

    // Combined char poly and normalization
    BitVect M_orig_;  // original M (before GCD reduction)
    BitVect M_;
    NTL::GF2X M_ntl_;
    NTL::GF2X inv_g0_;   // g_0^{-1} mod M

    // All normalized generating polynomials
    std::vector<BitVect> polys_;

    // StackBase for compute_all / compute_gap
    std::vector<std::unique_ptr<PolyLatBase>> stack_;

    // Incremental state for step()
    std::unique_ptr<PolyLatBase> current_base_;
    std::vector<std::unique_ptr<PolyLatBase>> opt_stack_;
    int previous_v_ = 0;

    // Compute generating polynomial for a SINGLE resolution index
    void find_single_poly(int idx, BitVect& out_poly);

    // Normalize a single polynomial: h = g * inv_g0 mod M
    void normalize_single(BitVect& poly);
};
