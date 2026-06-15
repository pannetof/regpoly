// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.

#pragma once
#include "f2_linear_source.h"
#include "i_testable.h"
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * @file digital_net.h
 * @brief Abstract base class for F_2 digital (t,m,s)-net generators.
 * @ingroup core
 *
 * A `DigitalNet` reinterprets the `Recurrence` interface as a fixed list
 * of `s_max` generating matrices `C_1, ..., C_{s_max}` over F_2, each of
 * size m x m. The point set of the net is
 *
 *   { (C_1 · i, C_2 · i, ..., C_{s_max} · i)  :  i in F_2^m }
 *
 * where i is the m-bit "input index" — what `state_` carries on this
 * subclass. The Recurrence-shaped contract maps to net operations as:
 *
 *   - `init(bv)`      sets the index to `bv` and resets the coordinate
 *                     counter `j_` to 1.
 *   - `next()`        advances the coordinate counter from `j_` to
 *                     `j_ + 1`. Throws `std::out_of_range` when called
 *                     with `j_ == s_max_` (no `C_{s_max + 1}` exists).
 *   - `get_output()`  returns `generating_matrix(j_) · state_` as an
 *                     m-bit `BitVect`. Throws when `j_` is outside
 *                     `[1, s_max_]`.
 *
 * With this contract, the existing equidistribution machinery
 * (`GaussMatrix::prepare`, `test_me_matricial`,
 * `dimension_equid`) computes the digital net's per-resolution
 * equidistribution dimension without modification, provided
 * `s_max >= kg`, where `kg` is the combined state size the kernel
 * passes as `indice_max`. For a J=1 single DigitalNet, kg = m, so
 * `s_max >= m` is enforced in the constructor; for combined
 * generators with J > 1, the ComboEnumerator/CombinedF2LinearSource builder
 * validates `s_max >= kg_combined` against every contained net.
 *
 * `init(bv == 0)` is legal for a digital net — index 0 is the first
 * point of the point set and is a meaningful seed. This relaxes the
 * non-zero contract documented on `Recurrence::init`. The
 * equidistribution kernel never actually seeds with zero (it walks
 * basis vectors), so the relaxation does not affect existing code.
 *
 * The t-value test (see `analyses/tvalue_runner.h`) reuses the same
 * `GaussMatrix::prepare`-built matrix and walks compositions to
 * compute the standard (t,m,s)-net t-value.
 *
 * Tempering chains and `CombinedF2LinearSource` wrapping are allowed:
 * tempering with an invertible matrix T yields a net `C_j' = T · C_j`
 * with the same t-value, and combining nets via CombinedF2LinearSource
 * produces a new digital net with stacked generating matrices —
 * a standard QMC construction (XOR-product of nets).
 */

namespace regpoly::core {

/**
 * @brief Abstract base for F_2 digital-net generators.
 *
 * Concrete subclasses (`SobolNet`, `NiederreiterF2Net`, ...) supply
 * the generating matrices via `generating_matrix(j)`. Subclasses
 * still implement the `Recurrence` pure-virtuals `name()`,
 * `display_str()`, and `copy()` themselves; `DigitalNet` itself
 * implements `init()`, `next()`, `get_output()`, and
 * `compute_default_test_method()`.
 *
 * @ingroup core
 */
class DigitalNet : public F2LinearSource {
public:
    /**
     * @brief Exception type for `next()` past `s_max`.
     *
     * Thrown when callers ask for a coordinate beyond the
     * construction-time `s_max` of the net. Derived from
     * `std::out_of_range` so callers can catch the base type.
     */
    struct DimensionExceededError : public std::out_of_range {
        using std::out_of_range::out_of_range;
    };

    /**
     * @brief Construct a digital net of dimension `m` with `s_max` coordinates.
     *
     * @param m      Input-index width in bits. Sets `k() = m` and
     *               `L() = m` on the base `Recurrence`.
     * @param s_max  Number of stored generating matrices.
     *               Must satisfy `s_max >= m` (the kernel passes
     *               `indice_max = kg = m` for a J=1 net; smaller
     *               `s_max` would cause `next()` to throw mid-prepare).
     *
     * @throws std::invalid_argument  When `m <= 0`, `s_max < m`,
     *                                or `m` exceeds the v1 cap of 64.
     */
    DigitalNet(int m, int s_max);

    // ── Recurrence interface, partially implemented ────────────────────

    /**
     * @brief Seed the input index and reset the coordinate counter.
     *
     * @param init_bv  Input index, width `k() == m`. Zero is legal
     *                 (the net's first point); the constructor's
     *                 non-zero contract is relaxed for digital nets.
     */
    void init(const BitVect& init_bv) override;

    /**
     * @brief Advance the coordinate counter from `j_` to `j_ + 1`.
     *
     * @throws DimensionExceededError  When called with `j_ == s_max_`.
     *                                 Callers must respect the
     *                                 `s_max >= kg` invariant.
     */
    void next() override;

    /**
     * @brief Return `generating_matrix(j_) * state_` as an m-bit `BitVect`.
     *
     * Overrides `F2LinearSource::raw_output()`. Digital nets carry no
     * tempering chain (chain is empty by default), so the
     * `F2LinearSource::get_output()` wrapper returns these bits
     * unchanged.
     *
     * @throws DimensionExceededError  When `j_` is outside `[1, s_max_]`.
     */
    BitVect raw_output() const override;

    // ── Net-specific accessors ────────────────────────────────────────

    /** @brief Maximum coordinate index (length of stored matrix list). */
    int s_max() const { return s_max_; }
    /** @brief Input-index width in bits (equals `k()`). */
    int m() const { return k_; }
    /** @brief Current coordinate counter (1..s_max valid; 0 if uninit). */
    int current_j() const { return j_; }

    /**
     * @brief Read-only view of the current input index `state_`.
     *
     * Phase 3: DigitalNet re-anchored to `F2LinearSource` (was a
     * `Recurrence` subclass pre-Phase-3). The `state_` accessor was
     * inherited from `Recurrence` then; now it's exposed here directly.
     */
    const BitVect& state() const { return state_; }

    /**
     * @brief Recommended test method for digital nets.
     *
     * - `"equidistribution"` → `"matricial"` (the matricial kernel
     *   computes the net's per-resolution equidistribution profile).
     * - `"tvalue"`           → `"schmid"`     (primal Schmid kernel).
     * - Anything else        → `std::nullopt`.
     *
     * Phase 3: direct override of `ITestable::default_test_method`.
     * No memoisation — the lookup is a constant table.
     */
    std::optional<std::string> default_test_method(
        const std::string& test_type) const override;

protected:
    /**
     * @brief Return row-wise access to the j-th generating matrix.
     *
     * Implementations return `m` rows (each a `BitVect` of width `m`).
     * Row r is the r-th row of `C_j` (0-based row index, 0-based
     * column index within the row).
     *
     * @param j  Coordinate index, 1-based; must satisfy `1 <= j <= s_max`.
     * @return   `const std::vector<BitVect>&` of length `m`.
     */
    virtual const std::vector<BitVect>& generating_matrix(int j) const = 0;

    BitVect state_;  ///< Input index (m bits). Pre-Phase-3 this was
                     ///< inherited from Recurrence::state_; now an
                     ///< explicit DigitalNet member since F2LinearSource
                     ///< carries no state.
    int j_;          ///< Coordinate counter, 1..s_max_ valid; 0 = uninitialised.
    int s_max_;      ///< Construction-time matrix-list length.
};

}  // namespace regpoly::core
