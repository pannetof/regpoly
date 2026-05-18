// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "combination.h"

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

// Phase 2.4b: equidistribution / collision-free / tuplets search loop
// in C++. Replaces the per-combo loop body of
// regpoly.search.seek.Seek.run(). Iterates the supplied Combination,
// for each combination + retry runs the configured tests in order,
// short-circuits when an equidistribution or tuplets test fails, and
// emits callbacks so Python can persist / display the result.
//
// nbtries > 1 is supported via an optional on_prep callback fired
// before every iteration; Python re-randomizes tempering parameters
// in that callback if needed. Test invocation is purely C++.
//
// Internally this loop builds polymorphic Test objects (see test.h)
// from the SeekTestSpec configs below — adding a new method or test
// type does not require editing this loop.

// Legacy enum kept for backward compatibility with Python callers that
// already speak in SeekTestKind. Equidistribution variants are now
// disambiguated by `SeekTestSpec::method_name` instead; the
// Equidistribution* enumerators are kept as aliases for SeekTestSpec
// to remain API-compatible.
namespace regpoly::core {

enum class SeekTestKind : uint8_t {
    Equidistribution                  = 0,
    CollisionFree                     = 6,
    Tuplets                           = 7,

    // ── Deprecated aliases (translate into Equidistribution + method_name) ──
    EquidistributionMatricial         = 0,
    EquidistributionLattice           = 1,
    EquidistributionHarase            = 2,
    EquidistributionNotPrimitive      = 3,
    EquidistributionSimdNotPrimitive  = 4,
    EquidistributionNothing           = 5,
};

struct SeekTestSpec {
    // Test type. The new canonical form uses {Equidistribution,
    // CollisionFree, Tuplets}; legacy values like
    // EquidistributionMatricial map to Equidistribution + method_name.
    SeekTestKind kind = SeekTestKind::Equidistribution;

    // For kind == Equidistribution: which method to use. Looked up via
    // MethodRegistry. When `kind` is one of the legacy enumerators
    // (EquidistributionMatricial / Lattice / Harase / NotPrimitive /
    // SimdNotPrimitive / Nothing) and method_name is empty, the loop
    // derives the method name from the enumerator.
    std::string method_name;

    // Equidistribution: the test object's own L (max resolution to
    // test); delta size eq_L_max_test+1; mse cap.
    int eq_L_max_test = 0;
    std::vector<int> eq_delta;
    int eq_mse = 0;

    // Tuplets: shape parameters (matches run_tuplets signature).
    int tup_d = 0;
    std::vector<int> tup_h;
    double tup_threshold = 0.0;
    int tup_testtype = 0;
};

// Per-iteration result handed to on_iter (only the relevant subset is
// populated based on which tests ran).
struct SeekIterResult {
    bool selected = false;     // passed all short-circuiting tests

    // Equidistribution snapshot (populated when an equidist test ran).
    bool me_ran = false;
    bool me_verified = false;
    bool me_is_me = false;     // se == 0 (full equidistribution)
    int me_se = 0;
    int me_test_L = 0;         // the test's L (echoed for caller)
    std::vector<int> me_ecart; // size me_test_L + 1; ecart[0] unused

    // Tuplets snapshot.
    bool tup_ran = false;
    bool tup_verified = false;
    bool tup_is_ok = false;
    int tup_firstpart_max = 0;
    int tup_firstpart_sum = 0;
    int tup_secondpart_max = 0;
    int tup_secondpart_sum = 0;

    // Collision-free snapshot.
    bool cf_ran = false;
    bool cf_verified = false;
    int cf_secf = 0;
    std::vector<int> cf_ecart_cf;  // size kg + 1
};

struct SeekProgress {
    int64_t nbgen;             // total iterations executed so far
    int64_t nb_select;         // selections so far
    int64_t nb_me;             // selections that hit full equidistribution
    double elapsed_seconds;
};

struct SeekResult {
    int64_t nbgen;
    int64_t nb_select;
    int64_t nb_me;
    double elapsed_seconds;
};

using SeekOnPrepFn =
    std::function<void(Combination&, bool /*is_retry*/)>;
using SeekOnIterFn =
    std::function<void(Combination&, const SeekIterResult&)>;
using SeekOnProgressFn =
    std::function<void(const SeekProgress&)>;

// Run the seek loop. `comb` must be reset()-able (Seek's caller
// already did this when binding YAML to comb). Each iteration:
//   1. on_prep(comb, is_retry) — optional, may re-randomize.
//   2. Run each SeekTestSpec in order on the live combination state.
//      EquidistributionTests can short-circuit (fail = stop running
//      remaining tests for this iteration). CollisionFree never
//      short-circuits; it always runs if listed, after the equidist
//      test's L is known.
//   3. If selected (all short-circuits passed and an equidist test
//      ran), call on_iter(comb, iter_result).
//   4. Periodically, call on_progress.
//   5. Advance via comb.next() when nbtries retries are exhausted.
SeekResult run_seek_search(
    Combination& comb,
    const std::vector<SeekTestSpec>& tests,
    int nbtries,
    int progress_interval,
    const SeekOnPrepFn& on_prep,
    const SeekOnIterFn& on_iter,
    const SeekOnProgressFn& on_progress);

}  // namespace regpoly::core
