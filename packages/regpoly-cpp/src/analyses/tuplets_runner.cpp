// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "tuplets_runner.h"

#include "gauss.h"
#include "transformation.h"

#include <algorithm>
#include <climits>
#include <vector>

namespace {

// Mirror Python's float(sys.maxsize) sentinel exactly: on 64-bit Python
// sys.maxsize == LLONG_MAX, and casting to float gives the same value
// as the C++ static_cast below.
const double PY_MAXSIZE_AS_DOUBLE =
    static_cast<double>(LLONG_MAX);

const double PY_NEG_MAXSIZE_AS_DOUBLE =
    -static_cast<double>(LLONG_MAX);

struct UnpackedGenerator {
    std::vector<Generator*> gens;
    std::vector<std::vector<Transformation*>> trans;
    std::vector<int> gen_k;
};

// Polymorphic unpack via Generator::components() / tempering_chains().
// Default (primitive) returns {*this} / {{}}; CombinedGenerator overrides.
UnpackedGenerator unpack_for_kernel(const Generator& gen) {
    UnpackedGenerator u;
    u.gens = gen.components();
    u.trans = gen.tempering_chains();
    u.gen_k.reserve(u.gens.size());
    for (auto* c : u.gens) u.gen_k.push_back(c->k());
    return u;
}

// Iterate over C(range_max - 1, count - 1) — i.e. all (count-1)-subsets
// of {1, 2, ..., range_max-1}, in lexicographic order. Matches Python's
// itertools.combinations(range(1, range_max), count - 1).
//
// The first index in the produced tuple is always 0 (caller prepends it),
// so this helper just yields the (count-1)-suffix.
//
// `cb` is invoked with the current subset (a vector of size count-1)
// each iteration. Returning false from `cb` aborts the enumeration.
template <typename Cb>
void for_each_combination(int range_max, int count, Cb cb) {
    if (count <= 0) {
        cb(std::vector<int>{});
        return;
    }
    if (range_max <= count) return;  // empty enumeration

    std::vector<int> picks(count);
    for (int i = 0; i < count; ++i) picks[i] = 1 + i;

    while (true) {
        if (!cb(picks)) return;

        int i = count - 1;
        while (i >= 0 && picks[i] == range_max - count + i) --i;
        if (i < 0) return;
        ++picks[i];
        for (int j = i + 1; j < count; ++j) picks[j] = picks[j - 1] + 1;
    }
}

}  // namespace

TupletsRunResult run_tuplets(
    const Generator& gen,
    int kg, int L,
    int tupd,
    const std::vector<int>& tuph,
    double threshold,
    int testtype)
{
    auto u = unpack_for_kernel(gen);

    // The C/Python algorithm picks indice_max = max(tuph[1..tupd]).
    int indice_max = 0;
    for (int i = 1; i <= tupd; ++i)
        indice_max = std::max(indice_max, tuph[i]);

    // One working matrix; mutated in place across resolution_equid
    // calls (matches the C reference and the Python code's behaviour
    // with a shared tMat object).
    GaussMatrix working = GaussMatrix::prepare(
        u.gens, u.gen_k, u.trans, kg, indice_max, L);

    TupletsRunResult out;
    out.tupd = tupd;
    out.tuph = tuph;
    out.gap = std::vector<double>(tuph[1] + 1, 0.0);
    out.DELTA = std::vector<double>(tupd + 1, PY_NEG_MAXSIZE_AS_DOUBLE);
    out.pourcentage = std::vector<double>(tupd + 1, 0.0);
    out.firstpart_max = PY_NEG_MAXSIZE_AS_DOUBLE;
    out.firstpart_sum = 0.0;
    out.secondpart_max = PY_NEG_MAXSIZE_AS_DOUBLE;
    out.secondpart_sum = 0.0;

    // ── First part: successive dimensions dim = 1..tuph[1] ─────────────
    int maxindsucc = std::min(tuph[1], kg);

    for (int dim = 1; dim <= maxindsucc; ++dim) {
        std::vector<int> indices(dim);
        for (int i = 0; i < dim; ++i) indices[i] = i;
        int l_t = working.resolution_equid(kg, dim, L, indices);
        out.gap[dim] = static_cast<double>(std::min(L, kg / dim) - l_t);

        out.firstpart_sum += out.gap[dim];
        if (out.gap[dim] > out.firstpart_max)
            out.firstpart_max = out.gap[dim];

        if (testtype == TUPLETS_TYPE_MAX && out.firstpart_max > threshold) {
            out.firstpart_max = PY_MAXSIZE_AS_DOUBLE;
            break;
        }
        if (testtype == TUPLETS_TYPE_SUM && out.firstpart_sum > threshold) {
            out.firstpart_sum = PY_MAXSIZE_AS_DOUBLE;
            break;
        }
    }

    // ── Second part: non-successive dimensions dim = 2..d ──────────────
    bool first_ok =
        (testtype == TUPLETS_TYPE_MAX && out.firstpart_max <= threshold) ||
        (testtype == TUPLETS_TYPE_SUM && out.firstpart_sum <= threshold);

    if (first_ok && tupd > 1) {
        for (int dim_idx = 2; dim_idx <= tupd; ++dim_idx) {
            out.DELTA[dim_idx] = PY_NEG_MAXSIZE_AS_DOUBLE;
            double bound = static_cast<double>(std::min(L, kg / dim_idx));
            int nbposs = 0;
            int nbeczero = 0;
            bool stop = false;

            for_each_combination(
                tuph[dim_idx], dim_idx - 1,
                [&](const std::vector<int>& rest) -> bool {
                    std::vector<int> indices;
                    indices.reserve(rest.size() + 1);
                    indices.push_back(0);
                    for (int x : rest) indices.push_back(x);

                    int l_t = working.resolution_equid(
                        kg, dim_idx, L, indices);
                    double gap_val = bound - static_cast<double>(l_t);

                    if (gap_val > out.DELTA[dim_idx])
                        out.DELTA[dim_idx] = gap_val;
                    if (gap_val == 0.0) ++nbeczero;
                    else out.secondpart_sum += gap_val;
                    ++nbposs;

                    if (testtype == TUPLETS_TYPE_MAX &&
                        out.DELTA[dim_idx] > threshold) {
                        out.secondpart_max = out.DELTA[dim_idx] =
                            PY_MAXSIZE_AS_DOUBLE;
                        stop = true;
                        return false;
                    }
                    if (testtype == TUPLETS_TYPE_SUM &&
                        out.secondpart_sum > threshold) {
                        out.secondpart_sum = out.DELTA[dim_idx] =
                            PY_MAXSIZE_AS_DOUBLE;
                        stop = true;
                        return false;
                    }
                    return true;
                });

            if (nbposs > 0)
                out.pourcentage[dim_idx] =
                    static_cast<double>(nbeczero) /
                    static_cast<double>(nbposs);

            if (out.secondpart_max < out.DELTA[dim_idx])
                out.secondpart_max = out.DELTA[dim_idx];

            if (stop) break;
        }
    } else {
        out.secondpart_max = PY_MAXSIZE_AS_DOUBLE;
        out.secondpart_sum = PY_MAXSIZE_AS_DOUBLE;
    }

    return out;
}
