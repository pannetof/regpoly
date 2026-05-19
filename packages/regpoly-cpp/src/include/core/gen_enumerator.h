// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "params.h"

#include <NTL/ZZ.h>

/**
 * @file gen_enumerator.h
 * @brief Abstract interface for the exhaustive-search generator-space enumerator.
 *
 * Defines `GenEnumerator` (the abstract interface the worker and
 * API talk to), the registry helpers `make_gen_enumerator` /
 * `family_is_enumerable`, and reusable combinatorial primitives
 * (`unrank_combination`, `mixed_radix_decode`, `binomial_zz`,
 * `parse_zz`, `zz_to_dec`) so each concrete enumerator does not
 * reimplement the istream/ostream dance.
 *
 * Counts are exposed as decimal strings because k-dependent
 * combinatorial totals routinely exceed `2^63`. Python bindings
 * lift them to arbitrary-precision ints.
 *
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief Abstract interface for the exhaustive-search generator-space enumerator.
 *
 * Concrete enumerators live next to the generator class they
 * enumerate (for example `TausworthePolyEnumerator` in
 * `gen_tausworthe.cpp`). The worker and API only ever talk to the
 * base class.
 *
 * Counts and indices are exposed as decimal strings because
 * k-dependent combinatorial totals routinely exceed `2^63`; Python
 * bindings lift them to arbitrary-precision ints.
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   auto en = make_gen_enumerator("TauswortheGen", resolved, 32);
 *   if (en) {
 *       auto total = en->size_dec();        // decimal string
 *       auto axes  = en->axes();
 *       Params p   = en->at("0");           // first combination
 *   }
 * @endcode
 *
 * @ingroup core
 */
class GenEnumerator {
public:
    /**
     * @brief Per-axis metadata for the UI summary; outer-most axis first.
     * @ingroup core
     */
    struct Axis {
        std::string name;        ///< Human-readable axis name.
        std::string size_dec;    ///< Decimal axis size (may exceed `2^63`).
        std::string describe;    ///< Human-readable hint for the UI.
    };

    virtual ~GenEnumerator() = default;

    /**
     * @brief Total number of admissible combinations as a decimal string.
     * @return  Decimal representation; may exceed `2^63`.
     */
    virtual std::string size_dec() const = 0;

    /**
     * @brief Concrete parameters for the `idx`-th combination.
     *
     * Must be admissible by construction. Deterministic and stable
     * across processes — same inputs produce the same `Params`.
     *
     * @param idx_dec  Index as a decimal string (may exceed `2^63`).
     * @return         Concrete `Params` for that combination.
     */
    virtual Params at(const std::string& idx_dec) const = 0;

    /**
     * @brief Per-axis metadata for the UI summary; outer-most axis first.
     * @return  Vector of `Axis` records.
     */
    virtual std::vector<Axis> axes() const = 0;
};

// ─── Registry (implemented in factory.cpp) ────────────────────────────────

/**
 * @brief Construct a `GenEnumerator` for a given family if one is registered.
 *
 * Routes through `GeneratorRegistry::find(family)->make_enumerator`.
 * Returns `nullptr` when the family has no exhaustive enumerator;
 * throws when the specific `resolved` inputs are insufficient
 * (e.g. `TauswortheGen` without `nb_terms`).
 *
 * @param family    Family name (canonical or alias).
 * @param resolved  Resolved `Params` after YAML / spec processing.
 * @param L         Output word width in bits.
 * @return          Heap-allocated enumerator, or `nullptr` if the
 *                  family is non-enumerable.
 * @throws std::invalid_argument  With a `needs_*` reason if `resolved`
 *                                is missing a required input.
 */
std::unique_ptr<GenEnumerator> make_gen_enumerator(
    const std::string& family, const Params& resolved, int L);

/**
 * @brief Probe whether a family has a registered exhaustive enumerator.
 *
 * @param family  Family name (canonical or alias).
 * @return        True iff the family is enumerable.
 *
 * @see :py:func:`regpoly.introspection.family_is_enumerable`
 */
bool family_is_enumerable(const std::string& family);

// ─── Reusable primitives for future enumerators ───────────────────────────

/**
 * @brief Lexicographic unranking of a `k`-subset of `{0, 1, ..., n-1}`.
 *
 * Returns the `idx`-th `k`-subset in lexicographic order via the
 * combinatorial number system. `idx` is accepted as an `NTL::ZZ`
 * so callers can pass large offsets without truncation.
 *
 * @param n    Universe size.
 * @param k    Subset size.
 * @param idx  Rank (0-based) within `C(n,k)` admissible subsets.
 * @return     Strictly ascending vector of `k` indices into `[0, n)`.
 * @throws std::invalid_argument  If `k < 0` or `k > n`.
 * @throws std::out_of_range      If `idx < 0` or `idx >= C(n,k)`.
 */
std::vector<int> unrank_combination(int n, int k, const NTL::ZZ& idx);

/**
 * @brief Mixed-radix decode: given axis sizes and a flat `idx`, return per-axis indices.
 *
 * The first axis varies slowest (`sizes[0]` is the outer digit).
 *
 * @param sizes  Per-axis sizes (outer-most first).
 * @param idx    Flat index in `[0, product(sizes))`.
 * @return       Per-axis indices, same length as `sizes`.
 * @throws std::out_of_range  If `idx < 0` or `idx >= product(sizes)`.
 */
std::vector<NTL::ZZ> mixed_radix_decode(
    const std::vector<NTL::ZZ>& sizes, const NTL::ZZ& idx);

/**
 * @brief Binomial coefficient `C(n, k)` as an `NTL::ZZ`.
 *
 * Returns 0 when `k < 0` or `k > n`. Used by enumerators that mix
 * combinatorial subsets with mixed-radix axes.
 *
 * @param n  Universe size.
 * @param k  Subset size.
 * @return   `C(n, k)` as `NTL::ZZ`; `0` when out of range.
 */
NTL::ZZ binomial_zz(int n, int k);

/**
 * @brief Parse a decimal string into an `NTL::ZZ`.
 *
 * Kept tiny because every concrete enumerator otherwise reimplements
 * the istream/ostream dance.
 *
 * @param dec  Decimal-text representation of an integer.
 * @return     Parsed `NTL::ZZ` value.
 * @throws std::invalid_argument  If `dec` is not a decimal integer.
 */
NTL::ZZ      parse_zz(const std::string& dec);

/**
 * @brief Format an `NTL::ZZ` as a decimal string.
 *
 * @param z  Value to format.
 * @return   Decimal-text representation.
 */
std::string  zz_to_dec(const NTL::ZZ& z);

}  // namespace regpoly::core
