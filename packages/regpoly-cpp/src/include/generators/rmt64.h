// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "generator.h"
#include "param_spec.h"
#include "params.h"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

/**
 * @file rmt64.h
 * @brief 64-bit Reducible Mersenne Twister (Saito 2011, Harase 2014).
 *
 * `RMT64Gen` is the 64-bit Reducible Mersenne Twister shipped as a
 * sample in MTToolBox (`samples/RMT/rmt64.hpp`). The "Reducible"
 * name signals that the characteristic polynomial of the F_2-linear
 * state is *not* primitive in general — the period is `2^mexp - 1`
 * within an invariant subspace of the `(mexp/64 + 1) * 64` F_2-linear
 * state. Equidistribution analysis must therefore use
 * `METHOD_NOTPRIMITIVE`. Catalog: `docs/library/rmt_params.yaml`
 * (cross-check parameter sets at multiple `mexp` values).
 *
 * State layout in `state_` (BitVect, MSB-first):
 *
 * @verbatim
 * bits   0..  63  = status[0]
 * bits  64.. 127  = status[1]
 *  ...
 * bits  64*(N-1)..64*N-1 = status[N-1]
 *   where N = mexp/64 + 1.
 * @endverbatim
 *
 * We do not roll the index into the BitVect; instead, like
 * MTToolBox, we keep an `index_` integer that points to the
 * most-recently-updated word.
 *
 * Recurrence (from `rmt64.hpp::generate`, sh1=29, sh2=17, sh3=37, sh4=43):
 *
 * @verbatim
 * index = (index + 1) % size
 * x = state[index]
 * y = state[(index + pos) % size]
 * y ^= (y << 17)
 * x = y ^ (x >> 1) ^ ((x & 1) ? mata : 0)
 * state[index] = x
 * // tempering:
 * x ^= (x >> 29) & 0x5555_5555_5555_5555
 * x ^= (x << 17) & maskb
 * x ^= (x << 37) & maskc
 * x ^= (x >> 43)
 * return x      // 64-bit output
 * @endverbatim
 *
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief 64-bit Reducible Mersenne Twister (Saito 2011, Harase 2014).
 *
 * Structural parameters: `mexp` (Mersenne exponent — the active
 * invariant-subspace dimension), `pos` (pick-up offset), `mata`
 * (twist matrix bottom row), and `maskb` / `maskc` (tempering
 * masks). Registered as `"RMT64Gen"` (alias `"RMT64"`).
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   Params p;
 *   p.set_int("mexp", 1279);
 *   p.set_int("pos", 1);
 *   p.set_int("mata",  0x121b78491deedbc8);
 *   p.set_int("maskb", 0xe3beff7ffbae0000);
 *   p.set_int("maskc", 0xfaf9fae000000000);
 *   auto gen = create_generator("RMT64Gen", p, 64);
 * @endcode
 *
 * @see :py:class:`regpoly.core.generator.Generator`
 * @ingroup core
 */
class RMT64Gen : public Generator {
public:
    /**
     * @brief Construct an RMT64Gen with explicit structural parameters.
     *
     * Most callers should go through `create_generator("RMT64Gen", ...)`
     * rather than this constructor directly.
     *
     * @param mexp   Mersenne exponent / invariant-subspace dimension.
     * @param pos    Pick-up offset used in the recurrence.
     * @param mata   Twist matrix bottom row.
     * @param maskb  First tempering mask.
     * @param maskc  Second tempering mask.
     * @param L      Output resolution in bits.
     */
    RMT64Gen(int mexp, int pos, uint64_t mata,
          uint64_t maskb, uint64_t maskc, int L);

    /**
     * @brief Build an RMT64Gen from a Params dict (registry factory hook).
     * @param params  Parameter dict with keys mexp, pos, mata, maskb, maskc.
     * @param L       Output resolution in bits.
     * @return        A constructed RMT64Gen as a polymorphic Generator pointer.
     * @throws std::runtime_error  If a required parameter is missing or invalid.
     */
    static std::unique_ptr<Generator> from_params(const Params& params, int L);

    /**
     * @brief Parameter specs declared by this family.
     * @return  Vector of ParamSpec records consumed by the factory and the
     *          Python introspection helpers.
     */
    static std::vector<ParamSpec> param_specs();

    /** @brief Family display name — returns the canonical "RMT64Gen" string. */
    std::string name() const override;
    /** @brief Human-readable parameter summary (mexp, pos, mata, masks). */
    std::string display_str() const override;
    /** @brief Seed the state from the leading `64 * N` bits of `init_bv`. */
    void init(const BitVect& init_bv) override;
    /** @brief Apply one RMT64 step (recurrence + tempering). */
    void next() override;
    /** @brief Deep copy this generator (state included). */
    std::unique_ptr<Generator> copy() const override;
    /** @brief Most recently produced tempered output truncated to `L()` bits. */
    BitVect get_output() const override;

    /** @brief Mersenne exponent (invariant-subspace dimension). */
    int mexp() const { return mexp_; }
    /** @brief State length in 64-bit words. */
    int N() const { return N_; }
    /** @brief Period exponent — equal to `mexp` for RMT64. */
    int period_exponent() const override { return mexp_; }

protected:
    /**
     * @brief Default test method — `notprimitive` (char poly is reducible).
     */
    std::optional<std::string> compute_default_test_method(const std::string& test_type) const override;

private:
    int mexp_;
    int N_;
    int pos_;
    uint64_t mata_, maskb_, maskc_;
    int index_;          // most-recently-updated word
    uint64_t last_output_;

    uint64_t get_word_u64(int i) const {
        return state_.get_word(i, 64);
    }
    void set_word_u64(int i, uint64_t v) {
        state_.set_word(i, 64, v);
    }
};

}  // namespace regpoly::core
