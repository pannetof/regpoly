// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "seek_search.h"

#include "combined.h"
#include "equidistribution_method.h"
#include "me_helpers.h"               // build_combined_from_combination
#include "test.h"

#include <chrono>
#include <climits>
#include <memory>
#include <stdexcept>
#include <vector>

using namespace regpoly::core;


namespace regpoly::core {

namespace {

// Translate a SeekTestSpec into a polymorphic Test instance. This is
// the ONLY translation site between the legacy spec API and the
// polymorphic dispatch — adding a new method or test type means
// editing the registries (factory.cpp's register_all_methods et al.),
// NOT this function.
std::unique_ptr<Test> spec_to_test(const SeekTestSpec& spec) {
    // Legacy enumerator → method_name mapping. Used when the caller
    // populated `kind` with a deprecated alias and left `method_name`
    // empty.
    auto resolve_method_name = [&]() -> std::string {
        if (!spec.method_name.empty()) return spec.method_name;
        switch (spec.kind) {
            case SeekTestKind::EquidistributionMatricial:        return "matricial";
            case SeekTestKind::EquidistributionLattice:          return "lattice";
            case SeekTestKind::EquidistributionHarase:           return "harase";
            case SeekTestKind::EquidistributionNotPrimitive:     return "notprimitive";
            case SeekTestKind::EquidistributionSimdNotPrimitive: return "simd_notprimitive";
            case SeekTestKind::EquidistributionNothing:          return "nothing";
            // Equidistribution (canonical) requires explicit method_name.
            default: return "matricial";  // historical default
        }
    };

    // Test-type discrimination collapses to three cases.
    if (spec.kind == SeekTestKind::CollisionFree) {
        return std::unique_ptr<Test>(new CollisionFreeTestRunner());
    }
    if (spec.kind == SeekTestKind::Tuplets) {
        return std::unique_ptr<Test>(new TupletsTestRunner(
            spec.tup_d, spec.tup_h,
            spec.tup_threshold, spec.tup_testtype));
    }

    // Everything else is an equidistribution variant.
    auto method = MethodRegistry::create(resolve_method_name());
    return std::unique_ptr<Test>(new EquidistributionTestRunner(
        std::move(method),
        spec.eq_L_max_test,
        spec.eq_delta,
        spec.eq_mse));
}

}  // namespace

SeekResult run_seek_search(
    Combination& comb,
    const std::vector<SeekTestSpec>& tests,
    int nbtries,
    int progress_interval,
    const SeekOnPrepFn& on_prep,
    const SeekOnIterFn& on_iter,
    const SeekOnProgressFn& on_progress)
{
    if (nbtries < 1) nbtries = 1;
    if (progress_interval < 1) progress_interval = 1;

    // Build polymorphic Test instances once per call. Each Test
    // owns its EquidistributionMethod (if applicable) and is reused
    // across iterations.
    std::vector<std::unique_ptr<Test>> runners;
    runners.reserve(tests.size());
    for (const auto& spec : tests) runners.push_back(spec_to_test(spec));

    int64_t nbgen = 0;
    int64_t nb_select = 0;
    int64_t nb_me = 0;
    int no_try = 1;
    const auto t_start = std::chrono::steady_clock::now();

    while (true) {
        ++nbgen;

        const bool is_retry = (no_try > 1);
        if (on_prep) on_prep(comb, is_retry);

        // Build a fresh CombinedGenerator from the current combination
        // state. Cheap — copies the active components + tempering chain.
        auto combined_owned = build_combined_from_combination(comb);
        Generator& combined = *combined_owned;

        SeekIterResult iter;
        bool passed = true;

        for (const auto& runner : runners) {
            runner->run(iter, combined, comb.k_g(), comb.L());
            if (!runner->passed(iter)) {
                passed = false;
                break;
            }
        }

        iter.selected = passed && iter.me_ran;
        if (iter.selected) {
            ++nb_select;
            if (iter.me_is_me) ++nb_me;
            if (on_iter) on_iter(comb, iter);
        }

        if (on_progress && (nbgen % progress_interval == 0)) {
            const auto now = std::chrono::steady_clock::now();
            on_progress(SeekProgress{
                nbgen, nb_select, nb_me,
                std::chrono::duration<double>(now - t_start).count()
            });
        }

        // Retry / advance.
        if (no_try < nbtries) {
            ++no_try;
        } else {
            if (!comb.next()) break;
            no_try = 1;
        }
    }

    const auto t_end = std::chrono::steady_clock::now();
    return SeekResult{
        nbgen, nb_select, nb_me,
        std::chrono::duration<double>(t_end - t_start).count()
    };
}

}  // namespace regpoly::core
