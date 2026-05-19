// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "generator.h"
#include "gen_enumerator.h"
#include "param_spec.h"
#include "params.h"
#include <vector>
#include <memory>
#include <string>
#include <unordered_map>
#include <cstdint>

/**
 * @file tausworthe.h
 * @brief Tausworthe LFSR family (Tausworthe 1965, L'Ecuyer 1996/1999).
 *
 * Defines `TauswortheGen`, the Tausworthe LFSR generator family
 * used in L'Ecuyer's combined LFSR113 / LFSR258 constructions
 * (`docs/library/lecuyer-1996.yaml`, `docs/library/lecuyer-1999.yaml`).
 * Also exposes the helper struct `RandomParamResult` carrying
 * sampled parameter values plus side-effect params the sampler
 * wants to splice back into the caller's params bag (used today by
 * TauswortheGen's `poly` sampler, which also chooses `s`).
 *
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief Sampled parameter value plus side-effect params (sampler helper).
 *
 * Returned by `TauswortheGen::generate_random` and similar sampler
 * methods. Either `int_val` (scalar) or `vec_val` (vector) is the
 * primary payload (selected by `is_vec`); `side_ints` carries
 * companion parameters the sampler implicitly chose so the caller
 * can propagate them to the generator constructor.
 *
 * Kept adjacent to the first family that needs it; hoist to a
 * shared header when a second family registers a sampler.
 *
 * @ingroup core
 */
struct RandomParamResult {
    bool is_vec = false;                                ///< True iff the primary payload is in `vec_val`.
    int64_t int_val = 0;                                ///< Sampled scalar value (when `is_vec` is false).
    std::vector<int> vec_val;                           ///< Sampled vector value (when `is_vec` is true).
    std::unordered_map<std::string, int64_t> side_ints; ///< Companion params the sampler also chose.
};

/**
 * @brief Tausworthe LFSR generator (Tausworthe 1965; L'Ecuyer 1996/1999).
 *
 * Structural parameters: `k` (state size in bits), `Q` (sorted
 * exponent list of the connection polynomial, including 0 and k),
 * `s` (decimation step), and `quicktaus` (uses the trinomial-table
 * fast path when admissible). Registered as `"TauswortheGen"`
 * (alias `"Tausworthe"`). LFSR113 wires four Tausworthe components
 * together; see `lecuyer-1999.yaml`.
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   Params p;
 *   p.set_int("k", 31);
 *   p.set_int("nb_terms", 3);
 *   p.set_bool("quicktaus", true);
 *   p.set_int_vec("poly", {0, 6, 31});
 *   p.set_int("s", 18);
 *   auto gen = create_generator("TauswortheGen", p, 32);
 * @endcode
 *
 * @see :py:class:`regpoly.core.generator.Generator`
 * @ingroup core
 */
class TauswortheGen : public Generator {
public:
    /**
     * @brief Construct a TauswortheGen with explicit polynomial and step.
     *
     * Most callers should go through `create_generator("TauswortheGen", ...)`
     * rather than this constructor directly.
     *
     * @param k          State size in bits (degree of the connection polynomial).
     * @param Q          Sorted exponent list `[0, q_1, ..., q_{t-2}, k]`.
     * @param s          Decimation step (number of state advances per `next()`).
     * @param quicktaus  Use the trinomial-table fast path when admissible.
     * @param L          Output resolution in bits.
     */
    TauswortheGen(int k, const std::vector<int>& Q, int s, bool quicktaus, int L);

    /**
     * @brief Build a TauswortheGen from a Params dict (registry factory hook).
     * @param params  Parameter dict with keys k, nb_terms, quicktaus, poly, s.
     * @param L       Output resolution in bits.
     * @return        A constructed TauswortheGen as a polymorphic Generator pointer.
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
     * @brief Sample a random admissible polynomial for a TauswortheGen.
     *
     * Returns the sorted exponent list `[0, q_1, ..., q_{t-2}, k]`.
     *
     * `s`: if 0, the caller wants auto-derivation
     * (`s = k - q_{t-2}`) — the top interior exponent may be
     * anywhere in `[1, k-1]`. If > 0, the poly must satisfy
     * `q_{t-2} <= k - s`.
     *
     * @param k          State size in bits.
     * @param nb_terms   Number of polynomial terms (must be odd, >= 3).
     * @param quicktaus  Whether the quicktaus fast path is requested.
     * @param L          Output resolution in bits.
     * @param s          Decimation step (0 to auto-derive).
     * @return           A sorted, admissible exponent list.
     * @throws std::invalid_argument  When the combination is inadmissible
     *         (nb_terms even or < 3, k <= 0, or quicktaus=true with k > L).
     */
    static std::vector<int> random_poly(
        int k, int nb_terms, bool quicktaus, int L, int s);

    /**
     * @brief Handle `rand_type` values owned by TauswortheGen.
     *
     * Supports `"tausworthe_s"` and `"tausworthe_poly"`. Reads `k`,
     * `nb_terms`, `quicktaus`, `poly`, `s` from `params` as needed.
     * For `"tausworthe_poly"` when `s` is unset, picks an admissible
     * `s` and returns it via `side_ints` so the caller can propagate
     * that exact value to the generator constructor (the poly was
     * sampled to match it).
     *
     * @param rand_type  Sampler name (`"tausworthe_s"` or `"tausworthe_poly"`).
     * @param rand_args  Optional sampler arguments (currently unused).
     * @param params     The current params bag.
     * @param L          Output resolution in bits.
     * @return           Sampled value plus companion `side_ints`.
     * @throws std::invalid_argument  For unsupported `rand_type` or invalid inputs.
     */
    static RandomParamResult generate_random(
        const std::string& rand_type,
        const std::string& rand_args,
        const Params& params, int L);

    /**
     * @brief Admissibility predicate for a (poly, s) pair.
     *
     * True iff the sorted polynomial `Q_sorted` together with
     * `(s, quicktaus, L)` satisfies the quicktaus inequality
     * `0 < s <= k - q_{t-2} < k <= L` (with `q_{t-2}` =
     * `Q_sorted[size-2]` = largest interior exponent), or, when
     * `quicktaus = false`, simply `1 <= s < k`. Also enforces the
     * odd-`nb_terms` requirement and the `Q[0] == 0 / Q[-1] == k`
     * conventions.
     *
     * @param Q_sorted   Sorted exponent list.
     * @param s          Decimation step.
     * @param quicktaus  Whether quicktaus is requested.
     * @param L          Output resolution in bits.
     */
    static bool is_admissible(
        const std::vector<int>& Q_sorted,
        int s, bool quicktaus, int L);

    /**
     * @brief Build the exhaustive-search enumerator for this family.
     *
     * The Params bag carries the structural inputs (`k`, `nb_terms`,
     * `quicktaus`) plus an optional `s` that, if present, fixes the
     * decimation step.
     *
     * @param resolved  Params with k, nb_terms, quicktaus pinned (s optional).
     * @param L         Output resolution in bits.
     * @return          A heap-allocated enumerator.
     * @throws std::invalid_argument  With a `"needs_<reason>"` message
     *         when a required axis is missing.
     */
    static std::unique_ptr<GenEnumerator> make_enumerator(
        const Params& resolved, int L);

    /** @brief Family display name — returns the canonical "TauswortheGen" string. */
    std::string name() const override;
    /** @brief Human-readable parameter summary (polynomial, s, quicktaus flag). */
    std::string display_str() const override;
    /** @brief Seed the state from the leading `k` bits of `init_bv`. */
    void init(const BitVect& init_bv) override;
    /** @brief Apply `s` state-advance steps via the quicktaus or general path. */
    void next() override;
    /** @brief Deep copy this generator (state included). */
    std::unique_ptr<Generator> copy() const override;
    /** @brief Characteristic polynomial (the trinomial / polynomial Q). */
    BitVect char_poly() const override;
    /** @brief Pack the canonical state into `out_words` (raw uint64_t form). */
    void get_transition_state(uint64_t* out_words, int out_nwords) const override;

    /** @brief Decimation step (number of advances per `next()`). */
    int s() const { return s_; }
    /** @brief True iff the quicktaus fast path is in use. */
    bool quicktaus() const { return quicktaus_; }
    /** @brief Sorted exponent list of the connection polynomial. */
    const std::vector<int>& Q() const { return Q_; }

private:
    std::vector<int> Q_;
    int NbCoeff_;
    int s_;
    bool quicktaus_;
    int gen_kms_;

    void next_quick();
    void next_general();
};

}  // namespace regpoly::core
