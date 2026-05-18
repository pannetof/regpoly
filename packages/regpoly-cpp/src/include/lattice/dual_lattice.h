// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "bitvect.h"
#include <vector>
#include <climits>
#include <algorithm>

// PolVect: a basis vector = an l-tuple of polynomials over GF(2).
//
// Storage: a single flat uint64_t array of size (maxresolution * nwords),
// laid out as
//   [coord_0_word_0, ..., coord_0_word_{nw-1},
//    coord_1_word_0, ..., coord_{maxresolution-1}_word_{nw-1}].
// This avoids per-coordinate heap allocations and improves cache behavior.
//
// `resolution` is the number of *currently active* coordinates l.  The
// data array is always sized for a fixed `maxresolution`; `resolution`
// is the prefix the algorithm walks at any given step (it grows from
// 1 to maxresolution as the dual lattice is built up).
//
// `deg` is the maximum degree across coordinates 0..resolution-1
// (INT_MIN if the vector is zero).  `indicemaxdeg` is the coordinate
// that achieves that max degree.

namespace regpoly::core {

struct PolVect {
    std::vector<uint64_t> data;   // flat storage: maxresolution * nwords
    int deg = INT_MIN;
    int indicemaxdeg = INT_MIN;
    int nwords = 0;               // words per coordinate
    int resolution = 0;           // current number of active coordinates

    static constexpr int WL = 64;

    PolVect() = default;
    // resolution_ = the *initial* active resolution. Storage is allocated
    // for resolution_ coordinates; if you need to grow beyond that, you
    // have to construct a wider PolVect. DualLatticeBase pre-sizes its
    // PolVects at maxresolution and then sets each one's `resolution`
    // field as the active basis grows.
    PolVect(int resolution_, int degmax);

    // Coordinate access.
    uint64_t*       coord(int j)       { return data.data() + j * nwords; }
    const uint64_t* coord(int j) const { return data.data() + j * nwords; }

    void swap(PolVect& other);
    // Mark the vector as conceptually zero without clearing the data
    // bytes. Cheaper than a full memset; per-coord setters
    // (set_coord_zero / set_coord_one / set_coord_from_bv) are then
    // responsible for cleaning the stale words on the slots they touch.
    void zero_lazy() { deg = INT_MIN; indicemaxdeg = INT_MIN; }

    // ── Per-coordinate setters ─────────────────────────────────────────
    // Clear coordinate j (zero its words from word 0 up to deg/WL).
    // No-op if the whole vector is currently lazy-zero.
    void set_coord_zero(int j);
    // Set coordinate j to the polynomial "1" (z^0). Updates deg if the
    // vector was previously zero.
    void set_coord_one(int j);
    // Copy a BitVect into coordinate j, then update deg if this
    // coordinate now holds a higher-degree polynomial than before.
    void set_coord_from_bv(int j, const BitVect& A);

    // ── Polynomial operations across all coordinates ───────────────────
    // Recompute deg / indicemaxdeg from the data words. Call after any
    // manual mutation that may have changed the leading degree.
    void update_deg();

    // *this = a * z^s. Requires *this and a to share resolution and
    // nwords (i.e. same maxresolution and same degmax). The result's
    // resolution matches *this's existing resolution.
    void assign_shifted(const PolVect& a, int s);

    // *this ^= a * z^s. Same shape requirement as assign_shifted.
    void add_shifted(const PolVect& a, int s);

    // Convenience: *this = a*z^s if I am zero, else *this ^= a*z^s.
    // The dual-lattice algorithm calls this as it accumulates the
    // running sum of shifted basis vectors.
    void add_self_mult(const PolVect& a, int s);

    // *this = a XOR b. As a side effect, b's tail bits below b.deg may
    // be masked off (preserves the original optimisation in
    // add_polvect that avoided clearing the unused suffix).
    void assign_xor(const PolVect& a, PolVect& b);
};

// DualLatticeBase: polynomial lattice basis for the dual lattice method.
//
// Owns:
//   - the basis (`vect_`),
//   - the column permutation maintained by Lenstra (`perm_`/`invperm_`),
//   - the algorithmic state (current resolution, max resolution, degmax).
//
// Per-PolVect operations (degree update, shift, XOR, etc.) live on
// PolVect itself. The only method here that touches a single PolVect
// is `element_norme_equal`, which indexes through `perm_` and is
// therefore basis-state-bound.
class DualLatticeBase {
public:
    DualLatticeBase(int maxresolution, int degmax);

    void dual_base(const std::vector<BitVect>& polys, const BitVect& M, int res);
    void dual_base_increase(const std::vector<BitVect>& polys);
    int lenstra(int res);

    int resolution() const { return resolution_; }
    int degmax() const { return degmax_; }

private:
    std::vector<PolVect> vect_;
    std::vector<int> perm_;
    std::vector<int> invperm_;
    int resolution_;
    int maxresolution_;
    int degmax_;
    int nwords_;                  // words per coordinate = ceil((degmax+1)/64)

    PolVect& vb(int i) { return vect_[i]; }
    const PolVect& vb(int i) const { return vect_[i]; }

    // Set the active resolution on every basis PolVect. Called whenever
    // resolution_ changes (dual_base / dual_base_increase) so each
    // PolVect knows how many of its coordinates are live.
    void set_all_polvect_resolutions(int r);

    // Read the bit at position `norme` of coordinate `perm_[j]` of P.
    // Couples to the basis's column permutation, so it stays here
    // rather than on PolVect.
    int element_norme_equal(const PolVect& P, int j, int norme) const;

    // Lenstra helpers
    void renumber(int m, int resolution);
    void solve_axb(int m, std::vector<int>& x);
    void permute_coord(int m);
};

}  // namespace regpoly::core
