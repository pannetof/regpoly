#pragma once

#include "generator.h"
#include <vector>

// Phase 2.3: Tuplets test orchestration in C++.
//
// Computes Delta(t_1, ..., t_d) for a combined generator's current state
// space. Mirrors regpoly.analyses.tuplets_test.TupletsTest.run().
//
// Conventions match the Python algorithm exactly:
//   - tuph is 1-indexed: tuph[0] is unused, tuph[1..d] holds the t values.
//   - testtype: 1 = MAX (matches _MAX_TYPE), 0 = SUM (_SUM_TYPE).
//   - DOUBLE_INF stand-in: std::numeric_limits<double>::max() is used in
//     places where the Python code uses sys.maxsize for double values
//     (firstpart_max, secondpart_max, DELTA). The Python sentinel is the
//     int sys.maxsize cast to float, ~9.22e18; the C++ counterpart uses
//     the same numeric value via DBL of LLONG_MAX so result comparisons
//     across languages stay byte-identical for the typical test
//     parameter regimes.

constexpr int TUPLETS_TYPE_SUM = 0;
constexpr int TUPLETS_TYPE_MAX = 1;

struct TupletsRunResult {
    int tupd;
    std::vector<int> tuph;          // 1-indexed (size tupd+1)
    std::vector<double> gap;         // 1-indexed (size tuph[1]+1)
    std::vector<double> DELTA;       // 1-indexed (size tupd+1)
    std::vector<double> pourcentage; // 1-indexed (size tupd+1)
    double firstpart_max;
    double firstpart_sum;
    double secondpart_max;
    double secondpart_sum;
};

TupletsRunResult run_tuplets(
    const Generator& gen,
    int kg, int L,
    int tupd,
    const std::vector<int>& tuph,
    double threshold,
    int testtype);
