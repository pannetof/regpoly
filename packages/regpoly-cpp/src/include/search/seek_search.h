// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "combination.h"

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

/**
 * @file seek_search.h
 * @brief Per-combo equidistribution / collision-free / tuplets search loop.
 * @ingroup core
 *
 * Phase 2.4b: equidistribution / collision-free / tuplets search loop
 * in C++. Replaces the per-combo loop body of
 * `regpoly.search.seek.Seek.run()`. Iterates the supplied
 * `Combination`, for each combination + retry runs the configured
 * tests in order, short-circuits when an equidistribution or tuplets
 * test fails, and emits callbacks so Python can persist / display
 * the result.
 *
 * `nbtries > 1` is supported via an optional `on_prep` callback fired
 * before every iteration; Python re-randomises tempering parameters
 * in that callback if needed. Test invocation is purely C++.
 *
 * Internally this loop builds polymorphic `Test` objects (see
 * `test.h`) from the `SeekTestSpec` configs below — adding a new
 * method or test type does not require editing this loop. The Python
 * `regpoly.search.seek.Seek` class wraps `run_seek_search` and is
 * responsible for YAML parsing and result presentation; the inner
 * search loop lives here.
 *
 * @see :py:class:`regpoly.search.seek.Seek`
 */

namespace regpoly::core {

/**
 * @brief Test-type discriminator for `SeekTestSpec`.
 *
 * The new canonical form uses the three high-level kinds
 * (`Equidistribution`, `CollisionFree`, `Tuplets`); the
 * `Equidistribution*` enumerators are deprecated aliases that
 * translate to `Equidistribution` plus a `method_name` selected by
 * the legacy enumerator name.
 */
enum class SeekTestKind : uint8_t {
    Equidistribution                  = 0,  ///< Equidistribution (method picked via `method_name`).
    CollisionFree                     = 6,  ///< Collision-free rank test (never short-circuits).
    Tuplets                           = 7,  ///< Tuplets uniformity test.

    // ── Deprecated aliases (translate into Equidistribution + method_name) ──
    EquidistributionMatricial         = 0,  ///< Deprecated: prefer `Equidistribution` + `method_name = "matricial"`.
    EquidistributionLattice           = 1,  ///< Deprecated: prefer `Equidistribution` + `method_name = "lattice"`.
    EquidistributionHarase            = 2,  ///< Deprecated: prefer `Equidistribution` + `method_name = "harase"`.
    EquidistributionNotPrimitive      = 3,  ///< Deprecated: prefer `Equidistribution` + `method_name = "notprimitive"`.
    EquidistributionSimdNotPrimitive  = 4,  ///< Deprecated: prefer `Equidistribution` + `method_name = "simd_notprimitive"`.
    EquidistributionNothing           = 5,  ///< Deprecated: prefer `Equidistribution` + `method_name = "nothing"`.
};

/**
 * @brief Configuration for one search test (built from YAML).
 *
 * The loop runs every `SeekTestSpec` in order against the live
 * combination state. Equidistribution and tuplets tests
 * short-circuit on failure; collision-free always runs (it never
 * short-circuits) and uses `eq_L_max_test` from the preceding
 * equidistribution test where applicable.
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   SeekTestSpec spec;
 *   spec.kind = SeekTestKind::Equidistribution;
 *   spec.method_name = "matricial";
 *   spec.eq_L_max_test = 32;
 *   spec.eq_delta.assign(33, 0);   // accept only se == 0
 *   spec.eq_mse = 0;
 *   std::vector<SeekTestSpec> tests{spec};
 *   // pass `tests` to run_seek_search.
 * @endcode
 *
 * @ingroup core
 */
struct SeekTestSpec {
    /**
     * @brief Test type. New canonical form: one of `Equidistribution`,
     *        `CollisionFree`, `Tuplets`. Legacy `Equidistribution*`
     *        enumerators are mapped onto `Equidistribution` plus an
     *        inferred `method_name`.
     */
    SeekTestKind kind = SeekTestKind::Equidistribution;

    /**
     * @brief For `kind == Equidistribution`: which method to use.
     *
     * Looked up via `MethodRegistry`. When `kind` is one of the legacy
     * enumerators (`EquidistributionMatricial` / `Lattice` / `Harase` /
     * `NotPrimitive` / `SimdNotPrimitive` / `Nothing`) and
     * `method_name` is empty, the loop derives the method name from
     * the enumerator.
     */
    std::string method_name;

    int eq_L_max_test = 0;     ///< Maximum resolution to test (equidistribution).
    std::vector<int> eq_delta; ///< Per-resolution gap budget (size `eq_L_max_test + 1`).
    int eq_mse = 0;            ///< Upper bound on the sum of equidistribution gaps.

    int tup_d = 0;                 ///< Tuplets dimension.
    std::vector<int> tup_h;        ///< Tuplets shape parameters (1-indexed; `tup_h[0]` unused).
    double tup_threshold = 0.0;    ///< Tuplets rejection threshold.
    int tup_testtype = 0;          ///< Tuplets test type (0 = SUM, 1 = MAX).
};

/**
 * @brief Per-iteration result handed to the seek loop's `on_iter` callback.
 *
 * Only the subset of fields relevant to the tests that actually ran
 * is populated; the rest stay at their defaults.
 *
 * @ingroup core
 */
struct SeekIterResult {
    bool selected = false;     ///< True iff every short-circuiting test passed.

    // Equidistribution snapshot (populated when an equidist test ran).
    bool me_ran = false;       ///< True iff an equidistribution test ran.
    bool me_verified = false;  ///< True iff the method certified its result.
    bool me_is_me = false;     ///< True iff `se == 0` (full equidistribution).
    int me_se = 0;             ///< Sum of equidistribution gaps.
    int me_test_L = 0;         ///< Echo of `eq_L_max_test` used.
    std::vector<int> me_ecart; ///< Per-resolution gap; size `me_test_L + 1`; `ecart[0]` unused.

    // Tuplets snapshot.
    bool tup_ran = false;             ///< True iff the tuplets test ran.
    bool tup_verified = false;        ///< True iff the tuplets kernel certified its result.
    bool tup_is_ok = false;           ///< True iff the tuplets test passed.
    int tup_firstpart_max = 0;        ///< Max gap on the first part.
    int tup_firstpart_sum = 0;        ///< Sum of gaps on the first part.
    int tup_secondpart_max = 0;       ///< Max gap on the second part.
    int tup_secondpart_sum = 0;       ///< Sum of gaps on the second part.

    // Collision-free snapshot.
    bool cf_ran = false;              ///< True iff the collision-free test ran.
    bool cf_verified = false;         ///< True iff the collision-free kernel certified its result.
    int cf_secf = 0;                  ///< Sum of collision-free gaps.
    std::vector<int> cf_ecart_cf;     ///< Per-dimension gap; size `kg + 1`.
};

/**
 * @brief Periodic progress snapshot for the seek loop.
 *
 * @ingroup core
 */
struct SeekProgress {
    int64_t nbgen;              ///< Total iterations executed so far.
    int64_t nb_select;          ///< Selections so far.
    int64_t nb_me;              ///< Selections that hit full equidistribution (`se == 0`).
    double elapsed_seconds;     ///< Wall-clock time since the driver started.
};

/**
 * @brief Loop-level summary returned by `run_seek_search`.
 *
 * @ingroup core
 */
struct SeekResult {
    int64_t nbgen;              ///< Total iterations executed.
    int64_t nb_select;          ///< Total selections.
    int64_t nb_me;              ///< Selections with `se == 0`.
    double elapsed_seconds;     ///< Wall-clock time.
};

/// Optional pre-iteration callback. Fires before every iteration; useful for re-randomising tempering.
using SeekOnPrepFn =
    std::function<void(Combination&, bool /*is_retry*/)>;
/// Per-iteration callback. Fires once after every iteration's test stack completes.
using SeekOnIterFn =
    std::function<void(Combination&, const SeekIterResult&)>;
/// Periodic progress callback.
using SeekOnProgressFn =
    std::function<void(const SeekProgress&)>;

/**
 * @brief Run the seek-search loop over a configured `Combination`.
 *
 * Each iteration:
 *  1. Fires `on_prep(comb, is_retry)` (optional, may re-randomise).
 *  2. Runs each `SeekTestSpec` in order against the live combination
 *     state. Equidistribution tests can short-circuit (fail = stop
 *     running remaining tests for this iteration). Collision-free
 *     never short-circuits; it always runs if listed, after the
 *     equidistribution test's `L` is known.
 *  3. If selected (all short-circuits passed and an equidistribution
 *     test ran), fires `on_iter(comb, iter_result)`.
 *  4. Periodically fires `on_progress`.
 *  5. Advances via `comb.next()` when `nbtries` retries are exhausted.
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   namespace yc = regpoly::yaml_config;
 *   auto cfg   = yc::load_seek_config("shared/yaml/equidist/mt19937.yaml");
 *   auto built = yc::build_search(cfg);
 *   run_seek_search(*built.combination, cfg.tests, cfg.nbtries,
 *                   100, nullptr,
 *                   [](Combination&, const SeekIterResult& r) {
 *                       // Persist r.
 *                   },
 *                   nullptr);
 * @endcode
 *
 * @param comb              Configured iterator (must be `reset()`-able).
 * @param tests             Test specs run in order on every iteration.
 * @param nbtries           Number of tries per combo (1 = no inner retry loop).
 * @param progress_interval Iterations between `on_progress` callbacks.
 * @param on_prep           Optional pre-iteration callback (may be null).
 * @param on_iter           Per-iteration callback for selected combos.
 * @param on_progress       Periodic progress callback.
 * @return                  Loop-level totals.
 *
 * @see :py:meth:`regpoly.search.seek.Seek.run`
 */
SeekResult run_seek_search(
    Combination& comb,
    const std::vector<SeekTestSpec>& tests,
    int nbtries,
    int progress_interval,
    const SeekOnPrepFn& on_prep,
    const SeekOnIterFn& on_iter,
    const SeekOnProgressFn& on_progress);

}  // namespace regpoly::core
