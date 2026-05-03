#pragma once
#include "bitvect.h"
#include "generator.h"
#include "transformation.h"
#include "me_helpers.h"   // MeLatResult
#include <vector>

// SIMD-aware non-primitive equidistribution test.
//
// Pipeline:
//   1–3. Recover χ_f, factor it, select largest Mersenne primitive φ
//        (shared with `test_me_notprimitive` via `chi_recovery.h`).
//   4.   Compute V-projected seeds via apply_polynomial(ψ, ·).
//   Fast: if `view_proto.simd_lane_count() == 1` AND p == kg, dispatch
//         to `test_me_lat`.  Otherwise fall through.
//   5.   For each (sm, wm) ∈ {0..elementNo-1} × {1..elementNo},
//        - advance generator clones by sm next() calls  (= start_mode);
//        - run a SIMD-PIS reduction (port of Saito's
//          AlgorithmSIMDEquidistribution.hpp) with the V-projected
//          seed as the PIS "rand" vector, where each PIS step reads
//          elementNo consecutive next() outputs and assembles them
//          as a 128-bit super-word with weight_mode blending against
//          a cached `previous` super-word, then extracts v MSBs of
//          each lane interleaved into the basis "next" BitVect;
//        - rescale e = min_count·elementNo − (elementNo − wm).
//        Aggregate per-sm MAX over wm, then MIN over sm.
//   6.   Emit ecart[v] = floor(p/v) − k(v).
//
// Generic across generators via the abstract Generator::simd_lane_count()
// hint.  For `simd_lane_count() == 1`, the algorithm collapses to
// behaviour identical to plain `test_me_notprimitive`.

MeLatResult test_me_notprimitive_simd(
    const std::vector<Generator*>& gens,
    const std::vector<std::vector<Transformation*>>& trans,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse);
