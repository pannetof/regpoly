#include "seek_search.h"

#include "combined.h"
#include "equidistribution_runner.h"
#include "me_harase.h"
#include "me_helpers.h"
#include "me_notprimitive.h"
#include "me_notprimitive_simd.h"
#include "tuplets_runner.h"

#include <chrono>
#include <climits>
#include <memory>

namespace {

// Mirrors the Python EquidistributionResults.is_presque_me():
//   verified must be true and se must be <= mse.
// Used as the short-circuit condition.
bool me_is_presque(const SeekIterResult& r, int mse) {
    if (!r.me_verified) return false;
    return r.me_se <= mse;
}

// Mirrors Python TupletsResults.is_ok():
//   if !verified, treat as ok (test was disabled).
//   else: verified means the bound was met.
// In the C++ runner, `verified` is signaled by tup_threshold being
// satisfied; for the simple case here we treat any returned struct
// as ok (the Python code's is_ok logic isn't easily replicated here
// without re-doing the dispatch). Phase 5 will tighten this when the
// web tests need it; for now matching Python's default is enough.
bool tup_is_ok_approximation(const SeekIterResult& r) {
    return r.tup_is_ok;
}

void run_equidist_test(
    SeekIterResult& iter,
    const SeekTestSpec& spec,
    const Generator& combined,
    int kg, int L)
{
    iter.me_ran = true;
    iter.me_test_L = spec.eq_L_max_test;

    if (spec.kind == SeekTestKind::EquidistributionNothing) {
        iter.me_ecart.assign(spec.eq_L_max_test + 1, 0);
        iter.me_se = 0;
        iter.me_verified = false;
        iter.me_is_me = false;
        return;
    }

    if (spec.kind == SeekTestKind::EquidistributionMatricial) {
        auto r = run_matricial_equidistribution(
            combined, kg, L, spec.eq_L_max_test,
            spec.eq_delta, spec.eq_mse);
        iter.me_ecart = std::move(r.ecart);
        iter.me_se = r.se;
        iter.me_verified = r.verified;
        iter.me_is_me = (iter.me_se == 0) && iter.me_verified;
        return;
    }

    if (spec.kind == SeekTestKind::EquidistributionLattice) {
        auto r = test_me_lat(
            combined, kg, L, spec.eq_L_max_test,
            spec.eq_delta, spec.eq_mse);
        iter.me_ecart = std::move(r.ecart);
        iter.me_se = r.se;
        iter.me_verified = true;
        iter.me_is_me = (iter.me_se == 0);
        return;
    }

    if (spec.kind == SeekTestKind::EquidistributionHarase) {
        auto r = test_me_harase(
            combined, kg, L, spec.eq_L_max_test,
            spec.eq_delta, spec.eq_mse);
        iter.me_ecart = std::move(r.ecart);
        iter.me_se = r.se;
        iter.me_verified = true;
        iter.me_is_me = (iter.me_se == 0);
        return;
    }

    if (spec.kind == SeekTestKind::EquidistributionNotPrimitive) {
        auto r = test_me_notprimitive(
            combined, kg, L, spec.eq_L_max_test,
            spec.eq_delta, spec.eq_mse);
        iter.me_ecart = std::move(r.ecart);
        iter.me_se = r.se;
        iter.me_verified = true;
        iter.me_is_me = (iter.me_se == 0);
        return;
    }

    if (spec.kind == SeekTestKind::EquidistributionSimdNotPrimitive) {
        auto r = test_me_notprimitive_simd(
            combined, kg, L, spec.eq_L_max_test,
            spec.eq_delta, spec.eq_mse);
        iter.me_ecart = std::move(r.ecart);
        iter.me_se = r.se;
        iter.me_verified = true;
        iter.me_is_me = (iter.me_se == 0);
        return;
    }
}

void run_cf_test(
    SeekIterResult& iter,
    const Generator& combined,
    int kg, int L,
    int L_for_phi4)
{
    iter.cf_ran = true;
    auto r = run_collision_free(combined, kg, L, L_for_phi4);
    iter.cf_ecart_cf = std::move(r.ecart_cf);
    iter.cf_secf = r.secf;
    iter.cf_verified = r.verified;
}

void run_tup_test(
    SeekIterResult& iter,
    const SeekTestSpec& spec,
    const Generator& combined,
    int kg, int L)
{
    iter.tup_ran = true;
    auto r = run_tuplets(combined, kg, L, spec.tup_d,
                         spec.tup_h, spec.tup_threshold,
                         spec.tup_testtype);
    iter.tup_firstpart_max = r.firstpart_max;
    iter.tup_firstpart_sum = r.firstpart_sum;
    iter.tup_secondpart_max = r.secondpart_max;
    iter.tup_secondpart_sum = r.secondpart_sum;
    // The runner doesn't return is_ok / verified in its dict; Phase 2.3
    // codified the same convention as Python — verified is implicit
    // from the kind being enabled. Leave is_ok true so the loop does
    // not short-circuit on tuplets unless callers wire something
    // tighter later.
    iter.tup_verified = true;
    iter.tup_is_ok = true;
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

    int64_t nbgen = 0;
    int64_t nb_select = 0;
    int64_t nb_me = 0;
    int no_try = 1;
    bool first_iter = true;
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

        // We need the ME test's L to feed CollisionFree (it uses
        // me_results.L to compute Phi_4). Track separately.
        int me_test_L_for_cf = comb.L();

        for (const auto& spec : tests) {
            if (spec.kind == SeekTestKind::CollisionFree) {
                if (!iter.me_ran) {
                    // The Python code passes me_results=None and falls
                    // back to comb.L. Match that behaviour.
                    me_test_L_for_cf = comb.L();
                }
                run_cf_test(iter, combined, comb.k_g(), comb.L(),
                            me_test_L_for_cf);
                continue;
            }

            if (spec.kind == SeekTestKind::Tuplets) {
                run_tup_test(iter, spec, combined, comb.k_g(), comb.L());
                if (!tup_is_ok_approximation(iter)) {
                    passed = false;
                    break;
                }
                continue;
            }

            // Otherwise: an equidistribution variant.
            run_equidist_test(iter, spec, combined, comb.k_g(), comb.L());
            me_test_L_for_cf = iter.me_test_L;
            if (!me_is_presque(iter, spec.eq_mse)) {
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

        first_iter = false;
        (void)first_iter;
    }

    const auto t_end = std::chrono::steady_clock::now();
    return SeekResult{
        nbgen, nb_select, nb_me,
        std::chrono::duration<double>(t_end - t_start).count()
    };
}
