// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.

#include "tvalue_runner.h"

#include "gauss.h"
#include "transformation.h"

#include <climits>
#include <stdexcept>
#include <vector>

namespace regpoly::core {

namespace {

struct UnpackedGenerator {
    std::vector<Generator*> gens;
    std::vector<std::vector<Transformation*>> trans;
    std::vector<int> gen_k;
};

// Mirrors equidistribution_runner.cpp's unpack_for_kernel.
UnpackedGenerator unpack_for_kernel(const Generator& gen) {
    UnpackedGenerator u;
    u.gens = gen.components();
    u.trans = gen.tempering_chains();
    u.gen_k.reserve(u.gens.size());
    for (auto* c : u.gens) u.gen_k.push_back(c->k());
    return u;
}

// Eliminate `len` columns starting at `col_start` against `M`, with
// `current_rank` as the next free pivot row. Updates `M` in place and
// returns the new rank. Pivot search and row XOR use the same
// primitives as the equidistribution kernel.
int eliminate_columns(GaussMatrix& M, int col_start, int len, int current_rank) {
    const int kg = M.nrows();
    for (int c = col_start; c < col_start + len; ++c) {
        if (current_rank >= kg) return current_rank;
        const int i = M.find_pivot(c, current_rank);
        if (i < 0) continue;
        if (i != current_rank) M.swap_rows(current_rank, i);
        M.eliminate_column(current_rank, c, current_rank + 1);
        ++current_rank;
    }
    return current_rank;
}

// Prefix-sharing branch-and-bound DFS over compositions of `remaining`
// into `s - idx` non-negative parts each in [0, m]. The Pirsic–Schmid
// (2001) algorithm for the t-value of a digital net: blocks 0..idx-1
// have already been eliminated into `working` (with rank
// `current_rank`); for each value of l_idx we clone `working`,
// eliminate block idx's first l_idx columns, and recurse.
//
// Two cuts (the speedup vs the naive enumeration):
//
//   1. Prefix reuse — block idx's `l_idx` elimination is shared across
//      every composition extending the (l_0, ..., l_{idx-1}) prefix.
//      Saves the dominant `O(d kg)` repeat-work each leaf used to do.
//
//   2. Rank-deficit branch-and-bound — each added column can raise the
//      rank by at most 1, so if `current_rank + remaining < d_target`
//      no descendant can hit `d_target`. Prune.
//
// Sibling siblings at the same depth ARE enumerated in increasing
// `l_idx` order, so the column set added between consecutive siblings
// differs by exactly one column — the Gray-code-style adjacency that
// the literature names. We do not exploit it for in-place rewind
// (which would require an undo-able Gauss tableau); we re-clone
// `working` per sibling because that keeps the implementation simple
// and the clone cost is dominated by the saved prefix work.
bool any_composition_full_rank_dfs(
    const GaussMatrix& working, int current_rank,
    int L, int m, int s, int d_target,
    int idx, int remaining)
{
    // Rank-deficit prune: best achievable rank from this state is
    // current_rank + min(remaining, free_rows_left). If even the best
    // case undershoots d_target, no child can reach it.
    const int free_rows = working.nrows() - current_rank;
    const int max_gain = std::min(remaining, free_rows);
    if (current_rank + max_gain < d_target) return false;

    if (idx == s) {
        return remaining == 0 && current_rank == d_target;
    }

    const int positions_left = s - 1 - idx;
    const int min_l = std::max(0, remaining - positions_left * m);
    const int max_l = std::min(remaining, m);

    for (int l = min_l; l <= max_l; ++l) {
        GaussMatrix M = working.copy();
        const int new_rank = eliminate_columns(M, idx * L, l, current_rank);
        if (any_composition_full_rank_dfs(
                M, new_rank, L, m, s, d_target, idx + 1, remaining - l)) {
            return true;
        }
    }
    return false;
}

}  // namespace

TValueResult run_tvalue_profile_schmid(
    const Generator& gen,
    int kg, int m, int s_max,
    const std::vector<int>& delta, int max_t_sum)
{
    if (s_max < 2) {
        // s = 1: trivially t(1) = 0 for any non-degenerate net.
        TValueResult res;
        res.tvals = std::vector<int>(std::max(s_max + 1, 2), 0);
        res.se = 0;
        res.verified = true;
        return res;
    }
    if (static_cast<int>(delta.size()) < s_max + 1) {
        throw std::invalid_argument(
            "run_tvalue_profile_schmid: delta size " + std::to_string(delta.size())
            + " is less than s_max + 1 (" + std::to_string(s_max + 1) + ")");
    }

    auto u = unpack_for_kernel(gen);
    // M has `kg` rows × `s_max * m` columns. Block j (1-based j; 0-based
    // index j-1) occupies columns [(j-1) * m, j * m).
    GaussMatrix M = GaussMatrix::prepare(u.gens, u.gen_k, u.trans,
                                          kg, /*indice_max=*/s_max, /*L=*/m);

    TValueResult res;
    res.tvals.assign(s_max + 1, 0);
    res.se = 0;
    res.verified = true;

    for (int s = 2; s <= s_max; ++s) {
        // Walk d from m down to 0 (m is the per-block column width; this
        // is also the input-dim bound for DigitalNets where kg = m).
        int d_max = 0;
        for (int d = m; d > 0; --d) {
            if (any_composition_full_rank_dfs(
                    M, /*current_rank=*/0, /*L=*/m, /*m=*/m,
                    s, /*d_target=*/d, /*idx=*/0, /*remaining=*/d)) {
                d_max = d;
                break;
            }
        }
        res.tvals[s] = m - d_max;
        res.se += res.tvals[s];

        if (res.tvals[s] > delta[s] || res.se > max_t_sum) {
            res.verified = false;
            break;
        }
    }
    return res;
}

TValueResult run_tvalue_profile_dual(
    const Generator& /*gen*/,
    int /*kg*/, int /*m*/, int /*s_max*/,
    const std::vector<int>& /*delta*/, int /*max_t_sum*/)
{
    throw std::runtime_error(
        "run_tvalue_profile_dual: Niederreiter–Pirsic dual method is "
        "registered but not yet implemented. Use 'schmid' for v1.");
}

}  // namespace regpoly::core
