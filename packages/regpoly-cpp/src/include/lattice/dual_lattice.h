// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "bitvect.h"
#include <vector>
#include <climits>
#include <algorithm>

/**
 * @file dual_lattice.h
 * @brief Couture-L'Ecuyer dual-lattice basis and Lenstra reduction.
 * @ingroup core
 *
 * Defines the polynomial-vector primitive `PolVect` and the
 * `DualLatticeBase` that drives the dual-lattice equidistribution
 * method. Selected via `MethodRegistry` by the name string
 * `"lattice"` in the YAML config.
 *
 * @verbatim
 *   PolVect: a basis vector = an l-tuple of polynomials over GF(2).
 *
 *   Storage: a single flat uint64_t array of size
 *     (maxresolution * nwords), laid out as
 *     [coord_0_word_0, ..., coord_0_word_{nw-1},
 *      coord_1_word_0, ..., coord_{maxresolution-1}_word_{nw-1}].
 *   This avoids per-coordinate heap allocations and improves cache
 *   behaviour.
 *
 *   `resolution` is the number of currently active coordinates l.
 *   The data array is always sized for a fixed `maxresolution`;
 *   `resolution` is the prefix the algorithm walks at any given step
 *   (grows from 1 to maxresolution as the dual lattice is built up).
 *
 *   `deg` is the maximum degree across coordinates 0..resolution-1
 *   (INT_MIN if the vector is zero). `indicemaxdeg` is the
 *   coordinate that achieves that max degree.
 * @endverbatim
 */

namespace regpoly::core {

/**
 * @brief Polynomial vector — one basis row of the dual lattice.
 *
 * Each `PolVect` is an `l`-tuple of polynomials over `GF(2)`, stored
 * as a flat `uint64_t` array of `maxresolution * nwords` words.
 * Coordinate access is via `coord(j)` (a span of `nwords` words).
 * `deg` and `indicemaxdeg` cache the leading-degree state for the
 * Lenstra reduction.
 *
 * @ingroup core
 */
struct PolVect {
    std::vector<uint64_t> data;   ///< Flat storage: `maxresolution * nwords` words.
    int deg = INT_MIN;            ///< Max degree across active coordinates (`INT_MIN` if zero).
    int indicemaxdeg = INT_MIN;   ///< Coordinate achieving `deg`.
    int nwords = 0;               ///< Words per coordinate.
    int resolution = 0;           ///< Current number of active coordinates.

    static constexpr int WL = 64;   ///< Word length in bits.

    PolVect() = default;
    /**
     * @brief Construct a `PolVect` with `resolution_` active coordinates
     *        and per-coordinate capacity for degree up to `degmax`.
     *
     * Storage is allocated for `resolution_` coordinates; if you need
     * to grow beyond that, construct a wider `PolVect`. `DualLatticeBase`
     * pre-sizes its `PolVect`s at `maxresolution` and then sets each
     * one's `resolution` field as the active basis grows.
     *
     * @param resolution_  Initial active resolution.
     * @param degmax       Maximum polynomial degree per coordinate.
     */
    PolVect(int resolution_, int degmax);

    /** @brief Coordinate accessor (mutable). */
    uint64_t*       coord(int j)       { return data.data() + j * nwords; }
    /** @brief Coordinate accessor (const). */
    const uint64_t* coord(int j) const { return data.data() + j * nwords; }

    /** @brief Swap contents with `other`. */
    void swap(PolVect& other);
    /**
     * @brief Mark the vector as conceptually zero without clearing the data bytes.
     *
     * Cheaper than a full `memset`. Per-coordinate setters
     * (`set_coord_zero` / `set_coord_one` / `set_coord_from_bv`) are
     * then responsible for cleaning the stale words on the slots
     * they touch.
     */
    void zero_lazy() { deg = INT_MIN; indicemaxdeg = INT_MIN; }

    // ── Per-coordinate setters ─────────────────────────────────────────
    /**
     * @brief Clear coordinate `j` (zero its words up to `deg/WL`).
     *
     * No-op if the whole vector is currently lazy-zero.
     *
     * @param j  Coordinate index.
     */
    void set_coord_zero(int j);
    /**
     * @brief Set coordinate `j` to the constant polynomial `1` (`z^0`).
     *
     * Updates `deg` if the vector was previously zero.
     *
     * @param j  Coordinate index.
     */
    void set_coord_one(int j);
    /**
     * @brief Copy a `BitVect` into coordinate `j`.
     *
     * Updates `deg` if this coordinate now holds a higher-degree
     * polynomial than before.
     *
     * @param j  Coordinate index.
     * @param A  Polynomial to write into the coordinate.
     */
    void set_coord_from_bv(int j, const BitVect& A);

    // ── Polynomial operations across all coordinates ───────────────────
    /**
     * @brief Recompute `deg` / `indicemaxdeg` from the data words.
     *
     * Call after any manual mutation that may have changed the
     * leading degree.
     */
    void update_deg();

    /**
     * @brief Assign `*this = a * z^s`.
     *
     * Requires `*this` and `a` to share `resolution` and `nwords`
     * (i.e. same `maxresolution` and same `degmax`). The result's
     * resolution matches `*this`'s existing resolution.
     *
     * @param a  Source polynomial vector.
     * @param s  Shift amount (in powers of `z`).
     */
    void assign_shifted(const PolVect& a, int s);

    /**
     * @brief Apply `*this ^= a * z^s` in place.
     *
     * Same shape requirement as `assign_shifted`.
     *
     * @param a  Source polynomial vector.
     * @param s  Shift amount.
     */
    void add_shifted(const PolVect& a, int s);

    /**
     * @brief Convenience: `*this = a*z^s` if I am zero, else `*this ^= a*z^s`.
     *
     * The dual-lattice algorithm calls this as it accumulates the
     * running sum of shifted basis vectors.
     *
     * @param a  Source polynomial vector.
     * @param s  Shift amount.
     */
    void add_self_mult(const PolVect& a, int s);

    /**
     * @brief Assign `*this = a XOR b`.
     *
     * As a side effect, `b`'s tail bits below `b.deg` may be masked
     * off (preserves the original optimisation in `add_polvect` that
     * avoided clearing the unused suffix).
     *
     * @param a  Left operand.
     * @param b  Right operand (may be mutated).
     */
    void assign_xor(const PolVect& a, PolVect& b);
};

/**
 * @brief Polynomial lattice basis for the dual-lattice method.
 *
 * Owns:
 *  - the basis (`vect_`),
 *  - the column permutation maintained by Lenstra
 *    (`perm_` / `invperm_`),
 *  - the algorithmic state (current resolution, max resolution,
 *    `degmax`).
 *
 * Per-`PolVect` operations (degree update, shift, XOR, etc.) live on
 * `PolVect` itself. The only method here that touches a single
 * `PolVect` is `element_norme_equal`, which indexes through `perm_`
 * and is therefore basis-state-bound.
 *
 * Selected via `MethodRegistry` by the name string `"lattice"`.
 *
 * @ingroup core
 */
class DualLatticeBase {
public:
    /**
     * @brief Construct a basis sized for `maxresolution` and `degmax`.
     *
     * @param maxresolution  Maximum resolution the basis can grow to.
     * @param degmax         Maximum polynomial degree per coordinate.
     */
    DualLatticeBase(int maxresolution, int degmax);

    /**
     * @brief Build the dual basis from a fresh set of generating polynomials.
     *
     * @param polys  Generating polynomials.
     * @param M      Combined characteristic polynomial.
     * @param res    Target resolution.
     */
    void dual_base(const std::vector<BitVect>& polys, const BitVect& M, int res);

    /**
     * @brief Extend an already-built basis by one resolution.
     *
     * Cheaper than `dual_base` from scratch when the algorithm walks
     * resolutions incrementally.
     *
     * @param polys  Generating polynomials for the new resolution.
     */
    void dual_base_increase(const std::vector<BitVect>& polys);

    /**
     * @brief Run Lenstra's reduction at the given resolution.
     *
     * @param res  Resolution to reduce.
     * @return     The dimension of the reduced basis at this resolution.
     */
    int lenstra(int res);

    /** @brief Current resolution. */
    int resolution() const { return resolution_; }
    /** @brief Per-coordinate `degmax` the basis was sized for. */
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
