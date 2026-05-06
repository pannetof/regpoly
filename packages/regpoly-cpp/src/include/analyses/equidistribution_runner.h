// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "generator.h"
#include <vector>

// Phase 2.3: equidistribution test orchestration in C++.
//
// run_matricial_equidistribution(gen, kg, L, Lmax, delta, mse) packages the
// outer loop that the Python EquidistributionTest.run() previously
// owned for the matricial method:
//
//   1. prepare_mat to build a GaussMatrix.
//   2. Walk l = 1..Lmax, calling dimension_equid(mat.copy(), kg, l, L).
//      Track ecart[l] = floor(kg/l) - dim.  Honour the "verif" state
//      machine (re-test missed resolutions when a non-zero gap appears
//      at l).  Bail out as soon as ecart[l] > delta[l] or se > mse.
//   3. Fill remaining unmeasured ecart slots with 0 (verified) or
//      INT_MAX (unverified beyond maxl) following the Python convention.
//
// The function consumes a single Generator& (typically a CombinedGenerator
// per Phase 1) and returns ecart, se, and verified.

struct MatricialEquidResult {
    std::vector<int> ecart;  // size Lmax+1; ecart[0] unused
    int se;
    bool verified;
};

MatricialEquidResult run_matricial_equidistribution(
    const Generator& gen,
    int kg, int L, int Lmax,
    const std::vector<int>& delta, int mse);

// Phase 2.3: collision-free orchestration in C++.
//
// run_collision_free(gen, kg, L, Lmax) walks t = kg..2 (descending),
// computing rank deficits via rang_cf for every t in Phi_4. Returns
// ecart_cf indexed 0..kg and the sum of gaps.

struct CollisionFreeResult {
    std::vector<int> ecart_cf;  // size kg+1
    int secf;
    bool verified;
};

// `L_for_phi4` is the L value used to compute Phi_4. In the existing
// Python orchestration this is `me_results.L` if available, otherwise
// `gen.L()`. The caller passes the value it would have computed.
CollisionFreeResult run_collision_free(
    const Generator& gen,
    int kg, int L, int L_for_phi4);
