// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "generator.h"
#include "param_spec.h"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

/**
 * @file melg.h
 * @brief Maximally Equidistributed Long-period Generator (Harase & Kimoto 2018).
 *
 * `MELGGen` implements the F_2-linear word recurrence of Harase &
 * Kimoto (2017/2018). The state consists of an array of `N - 1`
 * `w`-bit words plus one extra `w`-bit word `v` carrying a double
 * feedback. Period is `2^p - 1` where `p = N*w - r`. Catalog file:
 * `docs/library/harase-kimoto-2018.yaml`.
 *
 * Recurrence (Algorithm 1 from the paper):
 *
 * @verbatim
 * x  = (w[i] & upper_mask) ^ (w[(i+1)%(N-1)] & lower_mask)
 * v  = xA ^ w[(i+M)%(N-1)] ^ vB
 * w[i] = x ^ vC
 * i  = (i+1) % (N-1)
 *
 * xA = (x >> 1) ^ (MSB(x) ? a : 0)      (twist matrix)
 * vB = v ^ (v << sigma1)                 (left shift feedback)
 * vC = v ^ (v >> sigma2)                 (right shift feedback)
 * @endverbatim
 *
 * Output: `w[i]` (without tempering).
 *
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief Maximally Equidistributed Long-period generator (Harase & Kimoto 2018).
 *
 * Structural parameters: `w` (word width), `N` (state length in
 * words, including the `v` carry), `M` (mid-state offset), `r`
 * (lower-bit count), `sigma1` / `sigma2` (left / right shift counts
 * for the double feedback), and `a` (twist matrix bottom row).
 * Registered as `"MELGGen"` (alias `"MELG"`). The flagship
 * configuration is `melg19937-64`.
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   Params p;
 *   p.set_int("w", 64);
 *   p.set_int("N", 312);
 *   p.set_int("M", 81);
 *   p.set_int("r", 31);
 *   p.set_int("sigma1", 23);
 *   p.set_int("sigma2", 33);
 *   p.set_int("a", 0x5C32E06DF730FC42);
 *   auto gen = create_generator("MELGGen", p, 64);
 * @endcode
 *
 * @see :py:class:`regpoly.core.generator.Generator`
 * @ingroup core
 */
class MELGGen : public Generator {
public:
    /**
     * @brief Construct a MELGGen with explicit structural parameters.
     *
     * Most callers should go through `create_generator("MELGGen", ...)`
     * rather than this constructor directly.
     *
     * @param w       Word width in bits.
     * @param N       State length in words (including the v carry).
     * @param M       Mid-state offset.
     * @param r       Lower-bit count.
     * @param sigma1  Left-shift count for the v double-feedback.
     * @param sigma2  Right-shift count for the v double-feedback.
     * @param a       Twist matrix bottom row.
     * @param L       Output resolution in bits.
     */
    MELGGen(int w, int N, int M, int r, int sigma1, int sigma2,
         uint64_t a, int L);

    /**
     * @brief Build a MELGGen from a Params dict (registry factory hook).
     * @param params  Parameter dict with keys w, N, M, r, sigma1, sigma2, a.
     * @param L       Output resolution in bits.
     * @return        A constructed MELGGen as a polymorphic Generator pointer.
     * @throws std::runtime_error  If a required parameter is missing or invalid.
     */
    static std::unique_ptr<Generator> from_params(const Params& params, int L);

    /**
     * @brief Parameter specs declared by this family.
     * @return  Vector of ParamSpec records consumed by the factory and the
     *          Python introspection helpers.
     */
    static std::vector<ParamSpec> param_specs();

    /** @brief Family display name — returns the canonical "MELGGen" string. */
    std::string name() const override;
    /** @brief Human-readable parameter summary. */
    std::string display_str() const override;
    /** @brief Seed the state from `init_bv` (the first `N * w` bits are used). */
    void init(const BitVect& init_bv) override;
    /** @brief Advance the MELG recurrence by one step (rotates the ring head). */
    void next() override;
    /** @brief Deep copy this generator (state included). */
    std::unique_ptr<Generator> copy() const override;
    /** @brief Current output word `w[i]` truncated to `L()` bits. */
    BitVect get_output() const override;
    /** @brief Pack the canonical state into `out_words` (raw uint64_t form). */
    void get_transition_state(uint64_t* out_words, int out_nwords) const override;
    /** @brief Characteristic polynomial of the F_2-linear recurrence. */
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

}  // namespace regpoly::core
