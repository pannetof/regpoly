// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "generator.h"
#include "param_spec.h"
#include <memory>
#include <string>
#include <vector>

/**
 * @file cellular_automata.h
 * @brief Maximal-length null-boundary Elementary Cellular Automaton PRNG.
 *
 * `CellularAutomataGen` is the k-cell linear maximal-length Elementary
 * Cellular Automaton (ECA) of Bhuvaneswari & Bhattacharjee (arXiv:2603.19656,
 * 2026), restricted to the two linear rules over GF(2):
 *
 * @verbatim
 * Rule 90  on cell i: S_i' = S_{i-1} XOR S_{i+1}        (linear, balanced)
 * Rule 150 on cell i: S_i' = S_{i-1} XOR S_i XOR S_{i+1}
 * @endverbatim
 *
 * The state is packed in a k-bit BitVect with cell 0 at the MSB. One
 * CA step collapses to three bitvector ops:
 *
 * @verbatim
 * state <- (state << 1) XOR (state >> 1) XOR (rules & state)
 *          ^ brings cell i+1   ^ brings cell i-1  ^ middle bit only
 *          ^ to cell i         ^ to cell i        ^ at rule-150 cells
 * @endverbatim
 *
 * `next()` applies this `s` times in a row (time-spacing). T^s is
 * never materialised in the hot path; the matrix is only assembled
 * lazily by the base `Generator::transition_matrix()` helper when
 * the equidistribution machinery asks.
 *
 * Catalog: `docs/library/bhuvaneswari-bhattacharjee-2026-cellular-automata.yaml`.
 *
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief k-cell linear maximal-length ECA over rules 90 / 150 (Bhuvaneswari & Bhattacharjee 2026).
 *
 * Structural parameters: `k` (number of cells), `rule150_positions`
 * (0-indexed positions that use rule 150 — all others use rule 90),
 * and `s` (time-spacing, `s >= 1`; the paper uses `s = 1` for
 * diagnostics and `s ∈ [2, 10]` for the proposed combined PRNGs).
 * Registered as `"CellularAutomataGen"` (alias `"CA"`).
 *
 * CA's primitive characteristic polynomial path makes the matricial
 * equidistribution method the right default — even when wrapped in
 * a `CombinedGenerator`, the paper's convention treats the combined
 * matrix as a single full-rank target.
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   Params p;
 *   p.set_int("k", 32);
 *   p.set_int("s", 1);
 *   p.set_int_vec("rule150_positions",
 *       {1, 5, 6, 12, 15, 16, 18, 19, 20, 22, 23, 24, 25, 27, 29, 31});
 *   auto gen = create_generator("CellularAutomataGen", p, 32);
 * @endcode
 *
 * @see :py:class:`regpoly.core.generator.Generator`
 * @ingroup core
 */
class CellularAutomataGen : public Generator {
public:
    /**
     * @brief Construct a CellularAutomataGen with explicit rule layout.
     *
     * Most callers should go through `create_generator("CellularAutomataGen", ...)`
     * rather than this constructor directly.
     *
     * @param k                   Number of cells in the automaton.
     * @param rule150_positions   0-indexed cell positions that use rule 150;
     *                            all other cells use rule 90.
     * @param s                   Time-spacing (number of one-step CA updates
     *                            applied per `next()`; must be >= 1).
     * @param L                   Output resolution in bits.
     */
    CellularAutomataGen(int k,
                        const std::vector<int>& rule150_positions,
                        int s, int L);

    /**
     * @brief Build a CellularAutomataGen from a Params dict (registry factory hook).
     * @param params  Parameter dict with keys k, rule150_positions, s.
     * @param L       Output resolution in bits.
     * @return        A constructed CellularAutomataGen as a polymorphic Generator pointer.
     * @throws std::runtime_error  If a required parameter is missing or invalid.
     */
    static std::unique_ptr<Generator> from_params(const Params& params, int L);

    /**
     * @brief Parameter specs declared by this family.
     * @return  Vector of ParamSpec records consumed by the factory and the
     *          Python introspection helpers.
     */
    static std::vector<ParamSpec> param_specs();

    /** @brief Family display name — returns the canonical "CellularAutomataGen" string. */
    std::string name() const override;
    /** @brief Human-readable parameter summary (cells, rule layout, time-spacing). */
    std::string display_str() const override;
    /** @brief Seed the state from the leading `k` bits of `init_bv`. */
    void init(const BitVect& init_bv) override;
    /** @brief Apply `s` one-step CA updates (compose the rule-90/150 transition `s` times). */
    void next() override;
    /** @brief Deep copy this generator (state included). */
    std::unique_ptr<Generator> copy() const override;

    /** @brief Time-spacing (number of one-step updates applied per `next()`). */
    int s() const { return s_; }
    /** @brief Per-cell rule-150 mask (bit i = 1 iff cell i uses rule 150). */
    const BitVect& rules() const { return rules_; }
    /** @brief 0-indexed positions where rule 150 is active. */
    const std::vector<int>& rule150_positions() const { return rule150_positions_; }

protected:
    /**
     * @brief Default test-method override — picks `matricial` for CA.
     *
     * CA's primitive char poly path makes the matricial method the
     * right default — even when wrapped in a CombinedGenerator, the
     * paper's convention treats the combined matrix as a single
     * full-rank target (rank >= t·l check on (t+1)·l rows).
     */
    std::optional<std::string>
    compute_default_test_method(const std::string& test_type) const override;

private:
    int s_;
    std::vector<int> rule150_positions_;  // 0-indexed; kept for copy/display
    BitVect rules_;                        // k bits, bit i = 1 iff rule 150

    void next_one_step_();
};

}  // namespace regpoly::core
