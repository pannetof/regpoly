// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "generator.h"
#include "param_spec.h"
#include <cstdint>
#include <memory>
#include <string>

/**
 * @file mt.h
 * @brief Mersenne Twister generator family (Matsumoto & Nishimura 1998).
 *
 * MT is the canonical F_2-linear word-recurrence generator. The state
 * consists of `r` words of `w` bits; the recurrence twists the upper
 * `w - p` bits of `state[i]` with the lower `p` bits of `state[i+1]`
 * and XORs against `state[(i+m) % r]` through a twist matrix whose
 * bottom row is `a`. The factory name in the registry is `"MTGen"`
 * (alias: `"MersenneTwister"`).
 *
 * For the canonical MT19937 parameter set (`w=32, r=624, m=397, p=31,
 * a=0x9908B0DF`), the state has `k = w*r - p = 19937` effective bits
 * and the recurrence has Mersenne-prime period `2^19937 - 1`. See
 * the catalog file `docs/library/matsumoto-nishimura-1998.yaml` for
 * the published parameter set.
 *
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief Mersenne Twister recurrence over GF(2) (Matsumoto & Nishimura 1998).
 *
 * `MTGen` implements the raw F_2-linear state transition; tempering
 * (the standard MK output map) is applied as a separate
 * `Transformation` stacked on top of the generator when the catalog
 * entry includes a `tempering:` block. The catalog YAML
 * `matsumoto-nishimura-1998.yaml` carries the published `MT19937`
 * configuration including its `tempMK2` tempering parameters.
 *
 * Structural parameters: `w` (word width), `r` (state size in words),
 * `m` (mid-state offset used by the recurrence), `p` (number of
 * lower bits taken from `state[i+1]`), and `a` (bottom row of the
 * twist matrix). Output width is set by `L` (typically equal to `w`).
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   Params p;
 *   p.set_int("w", 32);
 *   p.set_int("r", 624);
 *   p.set_int("m", 397);
 *   p.set_int("p", 31);
 *   p.set_int("a", 0x9908B0DF);
 *   auto gen = create_generator("MTGen", p, 32);
 * @endcode
 *
 * @see :py:class:`regpoly.core.generator.Generator`
 * @ingroup core
 */
class MTGen : public Generator {
public:
    /**
     * @brief Construct an MTGen with explicit structural parameters.
     *
     * Most callers should go through `create_generator("MTGen", ...)`
     * rather than this constructor directly — the factory validates
     * parameters and integrates with the catalog / YAML loader.
     *
     * @param w  Word width in bits (typically 32).
     * @param r  State size in words.
     * @param m  Mid-state offset.
     * @param p  Lower-bit count for the upper-bit split.
     * @param a  Twist matrix bottom row.
     * @param L  Output resolution in bits.
     */
    MTGen(int w, int r, int m, int p, uint64_t a, int L);

    /**
     * @brief Build an MTGen from a Params dict (registry factory hook).
     * @param params  Parameter dict with keys w, r, m, p, a.
     * @param L       Output resolution in bits.
     * @return        A constructed MTGen as a polymorphic Generator pointer.
     * @throws std::runtime_error  If a required parameter is missing or invalid.
     */
    static std::unique_ptr<Generator> from_params(const Params& params, int L);

    /**
     * @brief Parameter specs declared by this family.
     * @return  Vector of ParamSpec records consumed by the factory and the
     *          Python introspection helpers.
     */
    static std::vector<ParamSpec> param_specs();

    /** @brief Family display name — returns the canonical "MTGen" string. */
    std::string name() const override;
    /** @brief Human-readable parameter summary for logs and result tables. */
    std::string display_str() const override;
    /** @brief Seed the state from `init_bv`; the first `w*r - p` bits are used. */
    void init(const BitVect& init_bv) override;
    /** @brief Apply one MT recurrence step (advance the ring head by one word). */
    void next() override;
    /** @brief Deep copy this generator (state included). */
    std::unique_ptr<Generator> copy() const override;
    /** @brief Characteristic polynomial of the F_2-linear recurrence. */
    BitVect char_poly() const override;
    /** @brief Current output word (the freshly-twisted state[i], pre-tempering). */
    BitVect get_output() const override;
    /** @brief Pack the canonical state into `out_words` (raw uint64_t form). */
    void get_transition_state(uint64_t* out_words, int out_nwords) const override;

    /** @brief Word width in bits. */
    int w() const { return w_; }
    /** @brief State size in words. */
    int r() const { return r_; }
    /** @brief Mid-state offset used by the recurrence. */
    int m() const { return m_; }
    /** @brief Lower-bit count for the upper-bit split. */
    int p() const { return p_; }
    /** @brief Twist matrix bottom row. */
    uint64_t a() const { return a_; }
    /** @brief Current ring-head index inside the state buffer. */
    int i_val() const { return i_; }

private:
    int w_, r_, m_, p_;
    uint64_t a_;
    int i_;
    uint64_t uu_, ll_;
    uint64_t maskw_;
    int state_bits_;

    uint64_t V(int idx) const;
    void SetV(int idx, uint64_t val);
    BitVect rotated_state() const;
};

}  // namespace regpoly::core
