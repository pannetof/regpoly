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
 * @file sobol.h
 * @brief Sobol digital net in base 2 (Sobol 1967; Antonov–Saleev recurrence;
 *        Joe–Kuo 2008 direction numbers).
 * @ingroup core
 *
 * `SobolNet` is a concrete `DigitalNet` whose j-th generating matrix C_j
 * is constructed by the Antonov–Saleev recurrence on direction numbers
 * driven by a primitive polynomial over F_2:
 *
 *  - Coordinate j=1 is the van der Corput / radical-inverse coordinate:
 *    C_1 is the m×m identity (column k = the standard basis vector e_k).
 *  - For j ≥ 2 the direction numbers `m_{j,1}, ..., m_{j,m}` satisfy
 *
 *       m_{j,k} = 2·a_1·m_{j,k-1} ⊕ 4·a_2·m_{j,k-2} ⊕ … ⊕
 *                 2^{s-1}·a_{s-1}·m_{j,k-s+1} ⊕
 *                 2^s·m_{j,k-s} ⊕ m_{j,k-s}    for k > s
 *
 *    where s = `deg(p_j)`, the `a_i` are p_j's interior coefficients
 *    encoded by the integer `a_j` (bit i of `a_j` = `a_i`, MSB first
 *    after the leading 1), and `m_{j,1..s}` are the initial direction
 *    numbers seeded from the table. The k-th column of C_j is the
 *    m-bit MSB-first encoding of `m_{j,k} << (m - k)`.
 *
 * v1 ships an embedded trusted table covering dimensions j=2..8
 * (Joe & Kuo 2003 Table 1; values reproduced in numerous downstream
 * implementations — QMCPy, Numerical Recipes §7.7, Wikipedia "Sobol").
 * Extending to the full 32-dimensional `new-joe-kuo-6.21201` table is
 * a follow-up commit requiring user-authorized fetch of that file.
 */

namespace regpoly::core {

/**
 * @brief Direction-number table row for one Sobol coordinate (j ≥ 2).
 *
 * `deg = 0` and `a = 0` with empty `m_init` is the sentinel that
 * coordinate j=1 (radical inverse) is handled in-class; callers do
 * not put j=1 in the table.
 */
struct SobolDirNumbers {
    int deg = 0;                          ///< Degree of the primitive polynomial.
    uint32_t a = 0;                       ///< Interior-coefficient bits, MSB-first
                                          ///< (bit `deg-1-i` set ⇔ a_i = 1).
    std::vector<uint32_t> m_init;         ///< Initial direction numbers m_{j,1..deg}.
};

/**
 * @brief Sobol digital net (base 2).
 *
 * @ingroup core
 */
class SobolNet : public DigitalNet {
public:
    /**
     * @brief Construct a Sobol net with explicit direction-number data.
     *
     * @param m      Input-index width (m-bit precision).
     * @param s_max  Coordinate count. Must satisfy `s_max >= m` (the
     *               `DigitalNet` invariant) and `s_max ≤ table.size() + 1`
     *               (coord 1 is built-in identity; coord j≥2 reads
     *               `table[j-2]`).
     * @param table  Direction-number rows for coordinates 2..s_max.
     *               `table[k]` carries the data for coordinate j = k+2.
     */
    SobolNet(int m, int s_max, const std::vector<SobolDirNumbers>& table);

    /**
     * @brief Construct a Sobol net from the embedded v1 table (j=2..8).
     *
     * @param m      Input-index width.
     * @param s_max  Coordinate count. Must satisfy `m ≤ s_max ≤ 8`
     *               for the embedded table; pass an explicit table
     *               via the other constructor to go beyond.
     */
    SobolNet(int m, int s_max);

    std::string name() const override { return "SobolNet"; }
    std::string display_str() const override;
    std::unique_ptr<Generator> copy() const override;

    /** @brief Return the embedded direction-number table size (v1 = 7 entries → j=2..8). */
    static int embedded_table_size();

protected:
    const std::vector<BitVect>& generating_matrix(int j) const override;

private:
    /// Shared, immutable list of m generating matrices (C_1..C_{s_max}),
    /// each stored as a vector of `m` row-`BitVect`s of width `m`.
    std::shared_ptr<const std::vector<std::vector<BitVect>>> matrices_;

    void build_matrices(int m, int s_max,
                        const std::vector<SobolDirNumbers>& table);
};

}  // namespace regpoly::core
