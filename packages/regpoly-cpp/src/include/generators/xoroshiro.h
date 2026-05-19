// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "generator.h"
#include "param_spec.h"
#include "params.h"

#include <cstdint>
#include <memory>
#include <string>

/**
 * @file xoroshiro.h
 * @brief Blackman & Vigna 2022 xoroshiro linear engine.
 *
 * `XoroshiroGen` implements the raw linear engine underlying the
 * xoroshiro family of generators (Blackman & Vigna 2022). The state
 * is `r` words of `w` bits each (`r >= 2`); parameters `(A, B, C)`
 * each in `[1, w-1]`. The catalog file
 * `docs/library/blackman-vigna-2022-scrambled.yaml` carries the
 * published `(A, B, C)` triples for each scrambler variant
 * (`+`, `*`, `++`, `**`).
 *
 * The update follows the `rw × rw` matrix from §3.1:
 *
 * @verbatim
 * new state = (old s_1, old s_2, ..., old s_{r-2},
 *              rotl(s_0, A) ^ t ^ (t << B),
 *              rotl(t, C))
 * where t = s_0 ^ s_{r-1}.
 * @endverbatim
 *
 * The scrambler (`+` / `++` / `*` / `**`) lives outside this class — this
 * is the raw linear engine, suitable for primitivity and
 * equidistribution analysis. The paper writes "k" for the word
 * count; this codebase reserves `k` for the state size in bits, so
 * we use `r` here (matching the variable name the paper uses
 * elsewhere, e.g. in Panneton & L'Ecuyer 2005).
 *
 * @ingroup core
 */

namespace regpoly::core {

class GenEnumerator;

/**
 * @brief xoroshiro linear engine (Blackman & Vigna 2022).
 *
 * Registered as `"XoroshiroGen"` (alias `"Xoroshiro"`). The
 * paper-recommended `xoroshiro128` engine uses `w = 64, r = 2,
 * A = 24, B = 16, C = 37`.
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   Params p;
 *   p.set_int("w", 64);
 *   p.set_int("r", 2);
 *   p.set_int("A", 24);
 *   p.set_int("B", 16);
 *   p.set_int("C", 37);
 *   auto gen = create_generator("XoroshiroGen", p, 64);
 * @endcode
 *
 * @see :py:class:`regpoly.core.generator.Generator`
 * @ingroup core
 */
class XoroshiroGen : public Generator {
public:
    /**
     * @brief Construct a XoroshiroGen with explicit structural parameters.
     *
     * Most callers should go through `create_generator("XoroshiroGen", ...)`
     * rather than this constructor directly.
     *
     * @param w  Word width in bits.
     * @param r  Number of words in the state (r >= 2).
     * @param A  Left-rotation amount in [1, w-1].
     * @param B  Left-shift amount in [1, w-1].
     * @param C  Left-rotation amount in [1, w-1].
     * @param L  Output resolution in bits.
     */
    XoroshiroGen(int w, int r, int A, int B, int C, int L);

    /**
     * @brief Build a XoroshiroGen from a Params dict (registry factory hook).
     * @param params  Parameter dict with keys w, r, A, B, C.
     * @param L       Output resolution in bits.
     * @return        A constructed XoroshiroGen as a polymorphic Generator pointer.
     * @throws std::runtime_error  If a required parameter is missing or invalid.
     */
    static std::unique_ptr<Generator> from_params(const Params& params, int L);

    /**
     * @brief Parameter specs declared by this family.
     * @return  Vector of ParamSpec records consumed by the factory and the
     *          Python introspection helpers.
     */
    static std::vector<ParamSpec> param_specs();

    /**
     * @brief Build the exhaustive-search enumerator.
     *
     * Requires `w` and `r` in `resolved` (the structural axes the
     * user pins). Throws `std::invalid_argument` with a
     * `"needs_<axis>"` message when an input is missing or out of
     * range. Enumeration order (outer first): `A × B × C`, each in
     * `[1, w-1]`.
     *
     * @param resolved  Params bag containing the pinned structural axes.
     * @param L         Output resolution in bits.
     * @return          A heap-allocated enumerator over the (A, B, C) cube.
     */
    static std::unique_ptr<GenEnumerator> make_enumerator(
        const Params& resolved, int L);

    /** @brief Family display name — returns the canonical "XoroshiroGen" string. */
    std::string name() const override;
    /** @brief Human-readable parameter summary (w, r, A, B, C). */
    std::string display_str() const override;
    /** @brief Seed the state from the leading `w * r` bits of `init_bv`. */
    void init(const BitVect& init_bv) override;
    /** @brief Advance the xoroshiro linear engine by one step. */
    void next() override;
    /** @brief Deep copy this generator (state included). */
    std::unique_ptr<Generator> copy() const override;

    /** @brief Word width in bits. */
    int w() const { return w_; }
    /** @brief Number of words in the state. */
    int r() const { return r_; }
    /** @brief Left-rotation amount A. */
    int A() const { return A_; }
    /** @brief Left-shift amount B. */
    int B() const { return B_; }
    /** @brief Left-rotation amount C. */
    int C() const { return C_; }

private:
    const int w_;
    const int r_;
    const int A_;
    const int B_;
    const int C_;
    const uint64_t wmask_;

    uint64_t rotl_(uint64_t x, int rot) const;
};

}  // namespace regpoly::core
