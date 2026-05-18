// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "generator.h"
#include "param_spec.h"
#include <memory>
#include <string>
#include <vector>

// CellularAutomataGen — k-cell linear maximal-length Elementary Cellular
// Automaton (ECA) over the null boundary, using only rules 90 and 150.
//
// Rule 90  on cell i: S_i' = S_{i-1} XOR S_{i+1}        (linear, balanced)
// Rule 150 on cell i: S_i' = S_{i-1} XOR S_i XOR S_{i+1}
//
// The state is packed in a k-bit BitVect with cell 0 at the MSB.
// One CA step collapses to three bitvector ops:
//
//     state ← (state << 1)   XOR   (state >> 1)   XOR   (rules & state)
//             ^ brings cell i+1     ^ brings cell i-1   ^ middle bit only
//             ^ to cell i           ^ to cell i         ^ at rule-150 cells
//
// next() applies this s_ times in a row (time-spacing). T^s is never
// materialized in the hot path; the matrix is only assembled lazily by
// the base Generator::transition_matrix() helper when equidistribution
// machinery asks.
//
// Reference: Bhuvaneswari A & Bhattacharjee K. "Cellular Automata based
// Resource Efficient Maximally Equidistributed Pseudo-Random Number
// Generators", arXiv:2603.19656, 2026.
class CellularAutomataGen : public Generator {
public:
    // rule150_positions: 0-indexed cell positions that use rule 150;
    //                    all other cells use rule 90.
    // s: time-spacing (s >= 1; paper uses 1 for diagnostics and 2..10
    //    for the proposed combined PRNGs).
    CellularAutomataGen(int k,
                        const std::vector<int>& rule150_positions,
                        int s, int L);

    static std::unique_ptr<Generator> from_params(const Params& params, int L);
    static std::vector<ParamSpec> param_specs();

    std::string name() const override;
    std::string display_str() const override;
    void init(const BitVect& init_bv) override;
    void next() override;
    std::unique_ptr<Generator> copy() const override;

    int s() const { return s_; }
    const BitVect& rules() const { return rules_; }
    const std::vector<int>& rule150_positions() const { return rule150_positions_; }

protected:
    // CA's primitive char poly path makes the matricial method the right
    // default — even when wrapped in a CombinedGenerator, the paper's
    // convention treats the combined matrix as a single full-rank target
    // (rank ≥ t·ℓ check on (t+1)·ℓ rows).
    std::optional<std::string>
    compute_default_test_method(const std::string& test_type) const override;

private:
    int s_;
    std::vector<int> rule150_positions_;  // 0-indexed; kept for copy/display
    BitVect rules_;                        // k bits, bit i = 1 iff rule 150

    void next_one_step_();
};
