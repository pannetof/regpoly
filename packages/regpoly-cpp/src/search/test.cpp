// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "test.h"

#include "equidistribution_runner.h"  // run_collision_free
#include "seek_search.h"              // SeekIterResult
#include "tuplets_runner.h"

#include <utility>

using namespace regpoly::core;


// ── EquidistributionTestRunner ──────────────────────────────────────────

namespace regpoly::core {

EquidistributionTestRunner::EquidistributionTestRunner(
    std::unique_ptr<EquidistributionMethod> method,
    int eq_L_max_test,
    std::vector<int> delta,
    int mse)
    : method_(std::move(method)),
      eq_L_max_test_(eq_L_max_test),
      delta_(std::move(delta)),
      mse_(mse)
{}

void EquidistributionTestRunner::run(
    SeekIterResult& iter,
    const Generator& gen,
    int kg, int L) const
{
    iter.me_ran = true;
    iter.me_test_L = eq_L_max_test_;

    auto r = method_->run(gen, kg, L, eq_L_max_test_, delta_, mse_);

    iter.me_ecart    = std::move(r.ecart);
    iter.me_se       = r.se;
    iter.me_verified = r.verified;
    iter.me_is_me    = (iter.me_se == 0) && r.verified;
}

bool EquidistributionTestRunner::passed(const SeekIterResult& iter) const {
    // Mirrors Python EquidistributionResults.is_presque_me():
    //   verified must be true and se must be <= mse.
    if (!iter.me_verified) return false;
    return iter.me_se <= mse_;
}

// ── CollisionFreeTestRunner ─────────────────────────────────────────────

void CollisionFreeTestRunner::run(
    SeekIterResult& iter,
    const Generator& gen,
    int kg, int L) const
{
    iter.cf_ran = true;
    // Phi_4 takes the equidistribution test's L if one ran, else
    // the generator's L. seek_search.cpp tracks me_test_L_for_cf
    // directly; here we read iter.me_test_L (set by an earlier
    // EquidistributionTestRunner::run) or fall back to L.
    const int L_for_phi4 = iter.me_ran ? iter.me_test_L : L;
    auto r = run_collision_free(gen, kg, L, L_for_phi4);
    iter.cf_ecart_cf = std::move(r.ecart_cf);
    iter.cf_secf     = r.secf;
    iter.cf_verified = r.verified;
}

// ── TupletsTestRunner ───────────────────────────────────────────────────

TupletsTestRunner::TupletsTestRunner(
    int d,
    std::vector<int> h,
    double threshold,
    int testtype)
    : d_(d), h_(std::move(h)),
      threshold_(threshold), testtype_(testtype)
{}

void TupletsTestRunner::run(
    SeekIterResult& iter,
    const Generator& gen,
    int kg, int L) const
{
    iter.tup_ran = true;
    auto r = run_tuplets(gen, kg, L, d_, h_, threshold_, testtype_);
    iter.tup_firstpart_max  = r.firstpart_max;
    iter.tup_firstpart_sum  = r.firstpart_sum;
    iter.tup_secondpart_max = r.secondpart_max;
    iter.tup_secondpart_sum = r.secondpart_sum;
    // The runner doesn't expose verified / is_ok directly; preserve
    // the prior seek_search convention: verified = true, is_ok = true.
    iter.tup_verified = true;
    iter.tup_is_ok    = true;
}

bool TupletsTestRunner::passed(const SeekIterResult& iter) const {
    return iter.tup_is_ok;
}

}  // namespace regpoly::core
