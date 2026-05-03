#pragma once

#include "temper_optimizer.h"
#include "transformation.h"

#include <cstdint>
#include <string>
#include <vector>

// Phase 2.4d: tempering bitmask optimization driver in C++.
//
// Replaces the recursive optimize(v) + run_minimize() body of
// packages/regpoly/src/regpoly/search/tempering_optimizer.py with a
// pure-C++ implementation that operates on:
//
//   * a TemperOptCache (already-existing C++ stateful kernel that
//     knows how to compute one resolution's gap incrementally), and
//   * a list of ParamLocator entries — one per optimizable bitmask
//     parameter (currently `b` for tempMK and `b` for laggedTempering).
//
// Safe-mask construction stays in Python (it is a small, structural
// computation that depends on `mu`/width and changes infrequently).
// The driver receives the precomputed safe_masks and just consumes them.

struct TemperParamLocator {
    Transformation* trans;       // not owned; lifetime is the caller's
    std::string param_name;      // e.g. "b"
    int width;                   // bit width (used by the perturbation step)
    int64_t current_value;       // shadow of trans's current value;
                                 // updated in place by the driver and
                                 // readable by the caller after the
                                 // driver returns.
};

struct TemperingOptimizerConfig {
    int max_essais;              // budget for one run_once() call
    std::vector<int> delta;      // size L+1; delta[v] is the gap budget
                                 // at resolution v. Use INT_MAX for
                                 // "minimize" mode.
    int mse;                     // upper bound on the cumulative SE
    int n_restarts;              // > 1 enables run_minimize() iterative
                                 // tightening; 1 = single run_once
    uint64_t random_seed;        // 0 = OS entropy
};

struct TemperingOptResult {
    int se;                      // sum of gaps for v = 1..L
    std::vector<int> gaps;       // size L+1; gaps[0] unused
    double elapsed_seconds;
    int essais;                  // number of perturbation iterations
                                 // consumed (run_minimize: number of
                                 // run_once calls)
};

// Run a single recursive optimization pass.
//
// safe_masks is shaped [L+1][P] where safe_masks[v][p] is the bitmask
// of bits that are "safe" to perturb in parameter p at resolution v.
// safe_masks[0] is unused.
//
// On entry, params[p].current_value must reflect the live value of
// trans->update(...) so the perturbation step can XOR against it. On
// return, the locators carry the best-found values and the
// underlying Transformation objects have those values applied.
TemperingOptResult run_tempering_optimizer_once(
    const TemperingOptimizerConfig& cfg,
    TemperOptCache& cache,
    std::vector<TemperParamLocator>& params,
    const std::vector<std::vector<uint64_t>>& safe_masks);

// Iterative delta-tightening loop (n_restarts > 1).
//
// Re-randomizes every locator's current_value, runs run_once with the
// previous best as the new (delta, mse) cap, accepts on improvement,
// counts consecutive non-improvements as failures, stops when
// failures == n_restarts.
TemperingOptResult run_tempering_optimizer_minimize(
    const TemperingOptimizerConfig& cfg,
    TemperOptCache& cache,
    std::vector<TemperParamLocator>& params,
    const std::vector<std::vector<uint64_t>>& safe_masks);
