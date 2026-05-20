// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.

#pragma once
#include "digital_net.h"
#include "bitvect.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

/**
 * @file niederreiter_f2.h
 * @brief Niederreiter digital sequence in F_2 (Niederreiter 1988).
 * @ingroup core
 *
 * `NiederreiterF2Gen` is a concrete `DigitalNet` whose j-th generating
 * matrix C_j is built from the j-th monic irreducible polynomial
 * p_j ∈ F_2[x] via the Laurent-series formula (Niederreiter 1988;
 * verified against Dick–Pillichshammer 2010 *Digital Nets and
 * Sequences* §8.1):
 *
 *   For e_j = deg(p_j), r ≥ 1, 0 ≤ k < e_j, take the expansion
 *     x^k / p_j(x)^r = Σ_{ℓ=0..∞} a^{(j)}(r, k, ℓ) · x^{-ℓ-1}
 *   in F_2((x⁻¹)). Then row i (1-based) of C_j has entry
 *     c_{j,i,ℓ} = a^{(j)}(Q+1, k, ℓ)
 *   where i - 1 = Q·e_j + k, 0 ≤ k < e_j.
 *
 * Irreducibles are enumerated on the fly at construction time in the
 * canonical order: p_1 = x, p_2 = x+1, then ascending degree, lex
 * within each degree by binary value with the leading 1 as the MSB.
 * Cost is milliseconds for any realistic `s_max`.
 *
 * **v1 limit**: m ≤ 63 (polynomial arithmetic uses `uint64_t`; the
 * highest polynomial degree handled is `m + e_max - 1`, which fits in
 * 64 bits when e_max ≤ ~9 and m ≤ 55 — generous for any v1 use). The
 * DigitalNet base allows m up to 64; the extra restriction lives here.
 */

namespace regpoly::core {

/**
 * @brief Niederreiter digital sequence in F_2.
 *
 * @ingroup core
 */
class NiederreiterF2Gen : public DigitalNet {
public:
    /**
     * @brief Construct a Niederreiter F_2 net with on-the-fly irreducibles.
     *
     * @param m      Input-index width in bits. Must satisfy `m <= 30`
     *               in v1 to keep polynomial arithmetic safely in
     *               `uint64_t` (polynomial degrees can reach `m + e_max`,
     *               where `e_max` grows with `s_max`).
     * @param s_max  Coordinate count, must satisfy `s_max >= m`.
     */
    NiederreiterF2Gen(int m, int s_max);

    std::string name() const override { return "NiederreiterF2Gen"; }
    std::string display_str() const override;
    std::unique_ptr<Generator> copy() const override;

    /** @brief j-th monic irreducible, packed as a uint64_t (bit i = coeff of x^i). */
    uint64_t irreducible(int j) const { return polys_[j - 1]; }
    /** @brief Degree of the j-th monic irreducible. */
    int irreducible_degree(int j) const { return degs_[j - 1]; }

protected:
    const std::vector<BitVect>& generating_matrix(int j) const override;

private:
    std::shared_ptr<const std::vector<std::vector<BitVect>>> matrices_;
    std::vector<uint64_t> polys_;   // first s_max monic irreducibles
    std::vector<int> degs_;          // their degrees

    void enumerate_irreducibles(int s_max);
    void build_matrices(int m, int s_max);
};

}  // namespace regpoly::core
