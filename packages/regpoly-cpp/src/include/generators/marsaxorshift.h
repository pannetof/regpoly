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
 * @file marsaxorshift.h
 * @brief Marsaglia xorshift family (Marsaglia 2003; Panneton & L'Ecuyer 2005).
 *
 * `MarsaXorshiftGen` is the unified Marsaglia-family xorshift
 * generator. It subsumes the older `XorShift128Gen` class and adds
 * paper-notation type tags plus an exhaustive-search enumerator.
 * Catalog file: `docs/library/marsaglia-2003-xorshift.yaml`.
 *
 * @ingroup core
 */

namespace regpoly::core {

// Forward declaration; full definition pulled in only by the .cpp.
// This keeps NTL out of every translation unit that includes this header.
class GenEnumerator;

/**
 * @brief Marsaglia xorshift family — multiple notation types and an enumerator.
 *
 * The `type` selector chooses one of the paper notations:
 *
 *   - **type 1**:   single-word xorshift, shifts `[a, b, c]`.
 *   - **type 2/2x**: two-component, three-shift triples
 *     `p = [p0, p1, p2], q = [q0, q1, q2]`.
 *   - **type 3**:   multi-tap, three taps with shifts.
 *   - **type 4**:   two-component, two-shift pairs `p, q`.
 *   - **type 100**: general multi-component (paper §3.1 "4-xorshift"
 *     form and beyond).
 *
 * Registered as `"MarsaXorshiftGen"`. The catalog's `xor()` triple
 * `(13, 17, 5)` is encoded as `type = 1, w = 32, r = 1,
 * shifts = [-13, 17, -5]`; the named `xor128()` is
 * `type = 2, w = 32, r = 4, m = 1, p = [-11, 8, 0], q = [19, 0, 0]`.
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   Params p;
 *   p.set_int("type", 1);
 *   p.set_int("w", 32);
 *   p.set_int("r", 1);
 *   p.set_int_vec("shifts", {-13, 17, -5});   // Marsaglia xor()
 *   auto gen = create_generator("MarsaXorshiftGen", p, 32);
 * @endcode
 *
 * @see :py:class:`regpoly.core.generator.Generator`
 * @ingroup core
 */
class MarsaXorshiftGen : public Generator {
public:
    /** @brief Type-1 parameter triple (single-word three-shift form). */
    struct Type1Params {
        int a, b, c;
    };

    /** @brief Type-2x parameter pair (two-component, 3 shifts each). */
    struct Type2xParams {
        std::vector<int> p;  ///< 3 values
        std::vector<int> q;  ///< 3 values
    };

    /** @brief One tap in a type-3 multi-tap configuration. */
    struct Tap {
        int position;   ///< 0-indexed slot the tap reads from.
        int shift;      ///< Signed shift amount (`< 0` = left, `> 0` = right).
    };

    /** @brief Type-4 parameter pair (two-component, 2 shifts each). */
    struct Type4Params {
        std::vector<int> p;  ///< 2 values
        std::vector<int> q;  ///< 2 values
    };

    /** @brief One entry in the type-100 general multi-component form. */
    struct MiEntry {
        int position;             ///< 0-indexed slot the entry reads from.
        std::vector<int> shifts;  ///< Per-entry shift list.
    };

    /**
     * @brief Construct a MarsaXorshiftGen with explicit per-type fields.
     *
     * Most callers should go through `create_generator("MarsaXorshiftGen", ...)`
     * rather than this constructor directly — the factory normalises
     * the type-specific keys (`shifts`, `p`, `q`, `mi_positions`, …)
     * into the right field set.
     *
     * @param type  Marsaglia notation type (1, 2, 3, 4, or 100).
     * @param w     Word width in bits.
     * @param r     Number of words in the state.
     * @param m     Pin for type-2/4 (active only for those types).
     * @param t1    Type-1 triple (active when `type == 1`).
     * @param t2x   Type-2/2x triples (active when `type == 2`).
     * @param taps  Type-3 tap list (active when `type == 3`).
     * @param t4    Type-4 pair (active when `type == 4`).
     * @param mi    Type-100 entry list (active when `type == 100`).
     * @param L     Output resolution in bits.
     */
    MarsaXorshiftGen(int type, int w, int r, int m,
                     const Type1Params& t1,
                     const Type2xParams& t2x,
                     const std::vector<Tap>& taps,
                     const Type4Params& t4,
                     const std::vector<MiEntry>& mi,
                     int L);

    /**
     * @brief Build a MarsaXorshiftGen from a Params dict (registry factory hook).
     * @param params  Parameter dict with the type-specific keys (see class docs).
     * @param L       Output resolution in bits.
     * @return        A constructed MarsaXorshiftGen as a polymorphic Generator pointer.
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
     * Requires `type` and `w` in `resolved`; types 2, 3, 4, 100
     * additionally require `r`; types 2 and 4 also accept an
     * optional `m` pin; type 100 accepts optional `nb_taps`
     * (default 3) and `shifts_per_tap` (per-tap vector, default
     * `[1, 1, ..., 1]` of length `nb_taps`). Throws
     * `std::invalid_argument` with a `"needs_<axis>"` message when
     * a required input is missing or inconsistent.
     *
     * Axes per type (outermost first; all coordinates are encoded
     * as mixed-radix digits over `NTL::ZZ` so the enumerator scales
     * to very large totals):
     *
     *   - type=1:    pattern (4) × a (w-1) × b (w-1) × c (w-1)
     *   - type=2:    m (r-1) × p[0..2] (2w-1) × q[0..2] (2w-1)
     *   - type=3:    tap-subset (C(r,3)) × sh[0..2] (2(w-1))
     *   - type=4:    m (r-1) × p[0..1] (2w-1) × q[0..1] (2w-1)
     *   - type=100:  tap-subset (C(r, nb_taps))
     *                × per-tap shift slots (2(w-1) each, total
     *                  sum(shifts_per_tap[i]) slots)
     *
     * @param resolved  Params with `type` and `w` pinned.
     * @param L         Output resolution in bits.
     * @return          A heap-allocated enumerator.
     */
    static std::unique_ptr<GenEnumerator> make_enumerator(
        const Params& resolved, int L);

    /** @brief Family display name — returns the canonical "MarsaXorshiftGen" string. */
    std::string name() const override;
    /** @brief Human-readable parameter summary (type, sizes, active fields). */
    std::string display_str() const override;
    /** @brief Seed the state from the leading `w * r` bits of `init_bv`. */
    void init(const BitVect& init_bv) override;
    /** @brief Apply one xorshift step under the active `type`. */
    void next() override;
    /** @brief Deep copy this generator (state included). */
    std::unique_ptr<Generator> copy() const override;

private:
    // All structural and configuration fields are set in the ctor and
    // never mutated afterwards.  Marking them `const` makes that
    // invariant compiler-enforced.
    const int type_;
    const int w_;
    const int r_;
    const int m_;
    const Type1Params t1_;
    const Type2xParams t2x_;
    const std::vector<Tap> taps_;
    const Type4Params t4_;
    const std::vector<MiEntry> mi_;
    const uint64_t wmask_;            // pre-computed (1 << w_) - 1; 0xFFF..FFFF when w_ >= 64

    static uint64_t ShiftR(uint64_t v, int s);
    uint64_t V(int idx) const;
    void SetV(int idx, uint64_t val);

    // BitVect-backed slow path used when w > 64 (the uint64_t kernel
    // cannot represent words wider than 64 bits).  Same recurrence
    // semantics as the fast path.  Defined in marsaxorshift.cpp.
    void next_wide_();
};

}  // namespace regpoly::core
