// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "generator.h"
#include "transformation.h"
#include <vector>
#include <cstdint>

/**
 * @file gauss.h
 * @brief Gaussian-elimination matrix over F_2 for equidistribution kernels.
 * @ingroup core
 *
 * `GaussMatrix` is the incremental F_2 Gaussian-elimination
 * accumulator used by the matricial equidistribution method and the
 * collision-free rank walks. Rows are stored as `nwords`
 * `uint64_t` words each; bit ordering matches `BitVect` (MSB-first
 * within each word).
 */

namespace regpoly::core {

/**
 * @brief Dense F_2 matrix with row XOR, pivot search, and incremental rank tracking.
 *
 * Used by the matricial equidistribution kernel (`dimension_equid`)
 * and the collision-free rank deficit walks (`rang_cf`). Construct
 * from `prepare(...)` (which collects sample outputs from the
 * generators) or directly via the constructor + row mutators.
 *
 * @code{.cpp}
 *   using regpoly::core::GaussMatrix;
 *   GaussMatrix mat = GaussMatrix::prepare(gens, gen_k, trans, kg, indice_max, L);
 *   int dim = mat.dimension_equid(kg, 1, L);   // l = 1
 *   // dim is the F_2 rank seen at resolution l.
 * @endcode
 *
 * @see :py:class:`regpoly.core.matrix.BitMatrix`
 *
 * @ingroup core
 */
class GaussMatrix {
public:
    /**
     * @brief Construct a zero matrix with `nrows` x `ncols` bits.
     * @param nrows  Number of rows.
     * @param ncols  Number of columns (bits).
     */
    GaussMatrix(int nrows, int ncols);

    /** @brief Return an independent deep copy. */
    GaussMatrix copy() const;
    /**
     * @brief Overwrite row `row` with raw words.
     * @param row    Row index.
     * @param words  Source word buffer.
     * @param n      Number of words to copy.
     */
    void set_row_from_words(int row, const uint64_t* words, int n);
    /**
     * @brief Read the bit at `(row, col)`.
     * @param row  Row index.
     * @param col  Column index (bit).
     * @return     True iff the bit is set.
     */
    bool bit_test(int row, int col) const;
    /**
     * @brief Swap two rows in place.
     * @param a  First row index.
     * @param b  Second row index.
     */
    void swap_rows(int a, int b);
    /**
     * @brief XOR row `src` into row `dst`.
     * @param dst  Destination row.
     * @param src  Source row.
     */
    void row_xor(int dst, int src);
    /**
     * @brief Eliminate column `col` below `pivot_row` (full row width).
     *
     * @param pivot_row  Pivot row index.
     * @param col        Column to eliminate.
     * @param start_row  First row to fold against (typically `pivot_row + 1`).
     */
    void eliminate_column(int pivot_row, int col, int start_row);
    /**
     * @brief Same as `eliminate_column` but restricted to `[col_lo, col_hi)`.
     *
     * Avoids touching bits outside the current resolution window.
     *
     * @param pivot_row  Pivot row index.
     * @param col        Column to eliminate.
     * @param start_row  First row to fold against.
     * @param col_lo     Inclusive low column.
     * @param col_hi     Exclusive high column.
     */
    void eliminate_column_masked(int pivot_row, int col, int start_row,
                                 int col_lo, int col_hi);
    /**
     * @brief Find a row with a 1 in column `col`, starting from `start_row`.
     *
     * @param col        Column index.
     * @param start_row  First row to inspect.
     * @return           Pivot row index, or `-1` if no pivot exists.
     */
    int find_pivot(int col, int start_row) const;

    /**
     * @brief Build the generator matrix from generators and tempering chains.
     *
     * Collects sample outputs from each component generator (taking
     * the tempering chain into account) and packs them as rows of an
     * F_2 matrix sized for the matricial equidistribution kernel.
     *
     * @param gens         Component generators.
     * @param gen_k        Per-component state widths.
     * @param trans        Per-component tempering chains.
     * @param kg           Combined state size.
     * @param indice_max   Maximum row index to prepare.
     * @param L            Output word width.
     * @return             A prepared `GaussMatrix`.
     */
    static GaussMatrix prepare(
        const std::vector<Generator*>& gens,
        const std::vector<int>& gen_k,
        const std::vector<std::vector<Transformation*>>& trans,
        int kg, int indice_max, int L);

    /**
     * @brief Matricial dimension at resolution `l`.
     *
     * Runs Gaussian elimination on the current row set restricted to
     * the first `l` output bits per row, returning the F_2 rank.
     *
     * @param kg  Combined state size.
     * @param l   Resolution.
     * @param L   Output word width.
     * @return    Equidistribution dimension `k(l)`.
     */
    int dimension_equid(int kg, int l, int L);
    /**
     * @brief Resolution-aware variant for the collision-free walk.
     *
     * @param kg       Combined state size.
     * @param t        Tuple dimension.
     * @param L        Output word width.
     * @param indices  Row indices to consider.
     * @return         F_2 rank.
     */
    int resolution_equid(int kg, int t, int L, const std::vector<int>& indices);
    /**
     * @brief Rank check for the collision-free test.
     *
     * @param kg  Combined state size.
     * @param t   Tuple dimension.
     * @param l   Resolution.
     * @param L   Output word width.
     * @return    Rank deficit `kg - rank`.
     */
    int rang_cf(int kg, int t, int l, int L);

    /** @brief Row count. */
    int nrows() const { return nrows_; }
    /** @brief Column count. */
    int ncols() const { return ncols_; }
    /** @brief Words per row. */
    int nwords() const { return nwords_; }
    /**
     * @brief Pointer to the raw words of `row`.
     * @param row  Row index.
     * @return     Pointer to the first word of the row.
     */
    const uint64_t* row_data(int row) const { return &data_[row * nwords_]; }

private:
    int nrows_, ncols_, nwords_;
    std::vector<uint64_t> data_;
    static constexpr int WL = 64;
};

}  // namespace regpoly::core
