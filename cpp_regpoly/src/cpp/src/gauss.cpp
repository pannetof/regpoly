#include "gauss.h"
#include <algorithm>
#include <cstring>

// ═══════════════════════════════════════════════════════════════════════════
// GaussMatrix — construction and basic operations
// ═══════════════════════════════════════════════════════════════════════════

GaussMatrix::GaussMatrix(int nrows, int ncols)
    : nrows_(nrows), ncols_(ncols),
      nwords_((ncols + WL - 1) / WL),
      data_(nrows * ((ncols + WL - 1) / WL), 0ULL) {}

GaussMatrix GaussMatrix::copy() const {
    GaussMatrix m(nrows_, ncols_);
    m.data_ = data_;
    return m;
}

void GaussMatrix::set_row_from_words(int row, const uint64_t* words, int n) {
    int base = row * nwords_;
    int cnt = std::min(n, nwords_);
    for (int i = 0; i < cnt; i++)
        data_[base + i] = words[i];
    for (int i = cnt; i < nwords_; i++)
        data_[base + i] = 0;
}

bool GaussMatrix::bit_test(int row, int col) const {
    int base = row * nwords_;
    int w = col / WL;
    int b = WL - 1 - (col % WL);
    return (data_[base + w] >> b) & 1;
}

void GaussMatrix::swap_rows(int a, int b) {
    int base_a = a * nwords_;
    int base_b = b * nwords_;
    for (int w = 0; w < nwords_; w++)
        std::swap(data_[base_a + w], data_[base_b + w]);
}

void GaussMatrix::row_xor(int dst, int src) {
    int base_d = dst * nwords_;
    int base_s = src * nwords_;
    for (int w = 0; w < nwords_; w++)
        data_[base_d + w] ^= data_[base_s + w];
}

void GaussMatrix::eliminate_column(int pivot_row, int col, int start_row) {
    int w_col = col / WL;
    uint64_t b_mask = 1ULL << (WL - 1 - (col % WL));
    int base_p = pivot_row * nwords_;

    for (int i = start_row; i < nrows_; i++) {
        int base_i = i * nwords_;
        if (data_[base_i + w_col] & b_mask) {
            for (int w = 0; w < nwords_; w++)
                data_[base_i + w] ^= data_[base_p + w];
        }
    }
}

void GaussMatrix::eliminate_column_masked(int pivot_row, int col, int start_row,
                                           int col_lo, int col_hi) {
    int w_col = col / WL;
    uint64_t b_mask = 1ULL << (WL - 1 - (col % WL));
    int base_p = pivot_row * nwords_;
    int w_lo = col_lo / WL;
    int w_hi = (col_hi - 1) / WL;

    int span = w_hi - w_lo + 1;
    std::vector<uint64_t> masks(span);
    std::vector<uint64_t> pwords(span);

    for (int w = w_lo; w <= w_hi; w++) {
        uint64_t m = ~0ULL;
        if (w == w_lo) {
            int bit_in_word = col_lo % WL;
            if (bit_in_word != 0)
                m &= (1ULL << (WL - bit_in_word)) - 1;
        }
        if (w == w_hi) {
            int tail = (w_hi + 1) * WL - col_hi;
            if (tail != 0)
                m &= ~((1ULL << tail) - 1);
        }
        int idx = w - w_lo;
        masks[idx] = m;
        pwords[idx] = data_[base_p + w] & m;
    }

    for (int i = start_row; i < nrows_; i++) {
        int base_i = i * nwords_;
        if (data_[base_i + w_col] & b_mask) {
            for (int idx = 0; idx < span; idx++)
                data_[base_i + w_lo + idx] ^= pwords[idx];
        }
    }
}

int GaussMatrix::find_pivot(int col, int start_row) const {
    int w = col / WL;
    uint64_t b_mask = 1ULL << (WL - 1 - (col % WL));
    for (int i = start_row; i < nrows_; i++) {
        if (data_[i * nwords_ + w] & b_mask)
            return i;
    }
    return -1;
}

// ═══════════════════════════════════════════════════════════════════════════
// GaussMatrix::prepare — Build the generator matrix
// ═══════════════════════════════════════════════════════════════════════════

static void extract_L_bits_into_row(
    const BitVect& state,
    const std::vector<Transformation*>& trans_chain,
    int L, int bit_offset,
    uint64_t* row_words, int row_nwords)
{
    BitVect out = state.copy();
    for (auto* t : trans_chain)
        t->apply(out);

    for (int b = 0; b < L; b++) {
        int src_word = b / 64;
        int src_bit = 63 - (b % 64);
        if (src_word < out.nwords() && ((out.data()[src_word] >> src_bit) & 1)) {
            int dst_pos = bit_offset + b;
            int dst_word = dst_pos / 64;
            int dst_bit = 63 - (dst_pos % 64);
            if (dst_word < row_nwords)
                row_words[dst_word] |= 1ULL << dst_bit;
        }
    }
}

GaussMatrix GaussMatrix::prepare(
    const std::vector<Generateur*>& gens,
    const std::vector<int>& gen_k,
    const std::vector<std::vector<Transformation*>>& trans,
    int kg, int indice_max, int L)
{
    int total_cols = indice_max * L;
    GaussMatrix mat(kg, total_cols);
    int row_nwords = (total_cols + 63) / 64;
    std::vector<uint64_t> row_words(row_nwords, 0);

    int ml = 0;
    for (int j = 0; j < (int)gens.size(); j++) {
        Generateur* gen = gens[j];
        int k = gen_k[j];
        const auto& chain = trans[j];
        bool has_trans = !chain.empty();

        BitVect bc(k);
        bc.set_bit(0, 1);

        for (int gl = 0; gl < k; gl++) {
            auto gen_copy = gen->copy();
            gen_copy->init(bc);
            bc.rshift(1);

            std::fill(row_words.begin(), row_words.end(), 0ULL);

            if (has_trans) {
                extract_L_bits_into_row(gen_copy->get_output(), chain,
                                        L, 0, row_words.data(), row_nwords);
            } else {
                BitVect out_bv = gen_copy->get_output();
                for (int w = 0; w < out_bv.nwords() && w < row_nwords; w++)
                    row_words[w] = out_bv.data()[w];
            }

            for (int mc = 1; mc < indice_max; mc++) {
                gen_copy->next();

                int bit_offset = mc * L;
                if (has_trans) {
                    extract_L_bits_into_row(gen_copy->get_output(), chain,
                                            L, bit_offset, row_words.data(), row_nwords);
                } else {
                    BitVect out_bv = gen_copy->get_output();
                    for (int b = 0; b < L; b++) {
                        int src_word = b / 64;
                        int src_bit = 63 - (b % 64);
                        if ((out_bv.data()[src_word] >> src_bit) & 1) {
                            int dst_pos = bit_offset + b;
                            int dst_word = dst_pos / 64;
                            int dst_bit = 63 - (dst_pos % 64);
                            if (dst_word < row_nwords)
                                row_words[dst_word] |= 1ULL << dst_bit;
                        }
                    }
                }
            }

            mat.set_row_from_words(ml, row_words.data(), row_nwords);
            ml++;
        }
    }
    return mat;
}

// ═══════════════════════════════════════════════════════════════════════════
// GaussMatrix::dimension_equid
// ═══════════════════════════════════════════════════════════════════════════

int GaussMatrix::dimension_equid(int kg, int l, int L) {
    int t = kg / l;
    int rang = 0;

    for (int j = 0; j < t; j++) {
        for (int cl = 0; cl < l; cl++) {
            int col = j * L + cl;
            int i = find_pivot(col, rang);
            if (i >= 0) {
                if (i != rang)
                    swap_rows(rang, i);
                eliminate_column_masked(rang, col, rang + 1, j * L, t * L);
                rang++;
            } else {
                return j;
            }
        }
    }
    return t;
}

// ═══════════════════════════════════════════════════════════════════════════
// GaussMatrix::resolution_equid
// ═══════════════════════════════════════════════════════════════════════════

int GaussMatrix::resolution_equid(int kg, int t, int L,
                                   const std::vector<int>& indices) {
    int l = std::min(kg / t, L);
    int rang = 0;

    for (int cl = 0; cl < l; cl++) {
        for (int j = 0; j < t; j++) {
            int col = indices[j] * L + cl;
            int i = find_pivot(col, rang);
            if (i >= 0) {
                if (i != rang)
                    swap_rows(rang, i);
                eliminate_column(rang, col, rang + 1);
                rang++;
            } else {
                return cl;
            }
        }
    }
    return l;
}

// ═══════════════════════════════════════════════════════════════════════════
// GaussMatrix::rang_cf
// ═══════════════════════════════════════════════════════════════════════════

int GaussMatrix::rang_cf(int /*kg*/, int t, int l, int L) {
    int rang = 0;

    for (int cl = 0; cl < l + 1; cl++) {
        for (int j = 0; j < t; j++) {
            int col = j * L + cl;
            int i = find_pivot(col, rang);
            if (i >= 0) {
                if (i != rang)
                    swap_rows(rang, i);
                eliminate_column(rang, col, rang + 1);
                rang++;
            }
        }
    }
    return rang;
}
