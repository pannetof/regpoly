#pragma once
#include "generateur.h"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

/**
 * MELG — Maximally Equidistributed F2-Linear Generator.
 *
 * Based on Harase & Kimoto (2017).  The state consists of an array
 * of (N-1) w-bit words plus one extra w-bit word v (double feedback).
 *
 * Recurrence (Algorithm 1 from the paper):
 *   x  = (w[i] & upper_mask) ^ (w[(i+1)%(N-1)] & lower_mask)
 *   v  = xA ^ w[(i+M)%(N-1)] ^ vB
 *   w[i] = x ^ vC
 *   i  = (i+1) % (N-1)
 *
 * Where:
 *   xA = (x >> 1) ^ (MSB(x) ? a : 0)      (twist matrix)
 *   vB = v ^ (v << sigma1)                   (left shift feedback)
 *   vC = v ^ (v >> sigma2)                   (right shift feedback)
 *
 * Period: 2^p - 1 where p = N*w - r.
 * Output: w[i] (without tempering).
 */
class MELG : public Generateur {
public:
    MELG(int w, int N, int M, int r, int sigma1, int sigma2,
         uint64_t a, int L);

    std::string name() const override;
    std::string display_str() const override;
    void init(const BitVect& init_bv) override;
    void next() override;
    std::unique_ptr<Generateur> copy() const override;
    BitVect get_output() const override;
    void get_transition_state(uint64_t* out_words, int out_nwords) const override;
    BitVect char_poly() const override;

private:
    int w_;
    int N_;
    int M_;
    int r_;
    int sigma1_;
    int sigma2_;
    uint64_t a_;
    int i_;                    // circular buffer pointer
    uint64_t upper_mask_;
    uint64_t lower_mask_;
    uint64_t word_mask_;
    int arr_size_;             // N-1
    int state_bits_;           // N*w (full storage including unused r bits)

    // Word access into state_: logical index 0..arr_size_-1 for the array,
    // index arr_size_ for v.  Stored in reversed order like MT.
    uint64_t V(int idx) const;
    void SetV(int idx, uint64_t val);

    // Rotate state to canonical form (pointer-independent).
    // Array part rotated by i_*w_, v appended at the end.
    BitVect rotated_state() const;
};
