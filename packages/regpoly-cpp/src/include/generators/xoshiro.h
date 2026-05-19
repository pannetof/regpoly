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
 * @file xoshiro.h
 * @brief Blackman & Vigna 2022 xoshiro linear engine.
 *
 * `XoshiroGen` implements the raw linear engine underlying the
 * xoshiro family of generators (Blackman & Vigna 2022). The state is
 * `r` words of `w` bits each, with `r ∈ {4, 8}`; parameters
 * `(A, B)` each in `[1, w-1]`. The catalog file
 * `docs/library/blackman-vigna-2022-scrambled.yaml` carries the
 * published `(A, B)` pairs for each scrambler variant.
 *
 * The update follows the `4w × 4w` (resp. `8w × 8w`) matrix from §3.2:
 *
 * @verbatim
 * r = 4 (S_4w):
 *   t = s[1] << A;  s[2] ^= s[0];  s[3] ^= s[1];
 *   s[1] ^= s[2];   s[0] ^= s[3];  s[2] ^= t;
 *   s[3] = rotl(s[3], B);
 *
 * r = 8 (S_8w):
 *   t = s[1] << A;
 *   s[2] ^= s[0];   s[5] ^= s[1];  s[1] ^= s[2];
 *   s[7] ^= s[3];   s[3] ^= s[4];  s[4] ^= s[5];
 *   s[0] ^= s[6];   s[6] ^= s[7];  s[6] ^= t;
 *   s[7] = rotl(s[7], B);
 * @endverbatim
 *
 * Scramblers (`+`, `++`, `**`) are NOT implemented — this is the
 * raw linear engine, suitable for primitivity and equidistribution
 * analysis. The paper writes "k" for the word count; this codebase
 * reserves `k` for the state size in bits, so we use `r`.
 *
 * @ingroup core
 */

namespace regpoly::core {

class GenEnumerator;

/**
 * @brief xoshiro linear engine (Blackman & Vigna 2022).
 *
 * Registered as `"XoshiroGen"` (alias `"Xoshiro"`). The paper's
 * `xoshiro256` engine uses `w = 64, r = 4, A = 17, B = 45`.
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   Params p;
 *   p.set_int("w", 64);
 *   p.set_int("r", 4);
 *   p.set_int("A", 17);
 *   p.set_int("B", 45);
 *   auto gen = create_generator("XoshiroGen", p, 64);
 * @endcode
 *
 * @see :py:class:`regpoly.core.generator.Generator`
 * @ingroup core
 */
class XoshiroGen : public Generator {
public:
    /**
     * @brief Construct a XoshiroGen with explicit structural parameters.
     *
     * Most callers should go through `create_generator("XoshiroGen", ...)`
     * rather than this constructor directly.
     *
     * @param w  Word width in bits.
     * @param r  Number of words in the state (must be 4 or 8).
     * @param A  Left-shift amount in [1, w-1].
     * @param B  Left-rotation amount in [1, w-1].
     * @param L  Output resolution in bits.
     */
    XoshiroGen(int w, int r, int A, int B, int L);

    /**
     * @brief Build a XoshiroGen from a Params dict (registry factory hook).
     * @param params  Parameter dict with keys w, r, A, B.
     * @param L       Output resolution in bits.
     * @return        A constructed XoshiroGen as a polymorphic Generator pointer.
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
     * Requires `w` and `r` in `resolved` (r must be 4 or 8). Throws
     * `std::invalid_argument` with a `"needs_<axis>"` message on bad
     * inputs. Enumeration order (outer first): `A × B`, each in
     * `[1, w-1]`.
     *
     * @param resolved  Params bag containing the pinned structural axes.
     * @param L         Output resolution in bits.
     * @return          A heap-allocated enumerator over the (A, B) plane.
     */
    static std::unique_ptr<GenEnumerator> make_enumerator(
        const Params& resolved, int L);

    /** @brief Family display name — returns the canonical "XoshiroGen" string. */
    std::string name() const override;
    /** @brief Human-readable parameter summary (w, r, A, B). */
    std::string display_str() const override;
    /** @brief Seed the state from the leading `w * r` bits of `init_bv`. */
    void init(const BitVect& init_bv) override;
    /** @brief Advance the xoshiro linear engine by one step. */
    void next() override;
    /** @brief Deep copy this generator (state included). */
    std::unique_ptr<Generator> copy() const override;

    /** @brief Word width in bits. */
    int w() const { return w_; }
    /** @brief Number of words in the state (4 or 8). */
    int r() const { return r_; }
    /** @brief Left-shift amount A. */
    int A() const { return A_; }
    /** @brief Left-rotation amount B. */
    int B() const { return B_; }

private:
    const int w_;
    const int r_;           // 4 or 8
    const int A_;
    const int B_;
    const uint64_t wmask_;

    uint64_t rotl_(uint64_t x, int rot) const;

    void next_4_();
    void next_8_();
};

}  // namespace regpoly::core
