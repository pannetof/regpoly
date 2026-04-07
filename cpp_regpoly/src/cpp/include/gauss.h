#pragma once
#include "generateur.h"
#include "transformation.h"
#include <vector>
#include <cstdint>

class GaussMatrix {
public:
    GaussMatrix(int nrows, int ncols);

    GaussMatrix copy() const;
    void set_row_from_words(int row, const uint64_t* words, int n);
    bool bit_test(int row, int col) const;
    void swap_rows(int a, int b);
    void row_xor(int dst, int src);
    void eliminate_column(int pivot_row, int col, int start_row);
    void eliminate_column_masked(int pivot_row, int col, int start_row,
                                 int col_lo, int col_hi);
    int find_pivot(int col, int start_row) const;

    // Build the generator matrix from generators and transformations
    static GaussMatrix prepare(
        const std::vector<Generateur*>& gens,
        const std::vector<int>& gen_k,
        const std::vector<std::vector<Transformation*>>& trans,
        int kg, int indice_max, int L);

    // Gaussian elimination algorithms
    int dimension_equid(int kg, int l, int L);
    int resolution_equid(int kg, int t, int L, const std::vector<int>& indices);
    int rang_cf(int kg, int t, int l, int L);

    int nrows() const { return nrows_; }
    int ncols() const { return ncols_; }
    int nwords() const { return nwords_; }
    const uint64_t* row_data(int row) const { return &data_[row * nwords_]; }

private:
    int nrows_, ncols_, nwords_;
    std::vector<uint64_t> data_;
    static constexpr int WL = 64;
};
