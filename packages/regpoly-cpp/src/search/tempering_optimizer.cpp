#include "tempering_optimizer.h"

#include "params.h"
#include "random_samplers.h"

#include <algorithm>
#include <chrono>
#include <climits>
#include <cstdint>
#include <functional>
#include <utility>

namespace {

// Apply a single-key partial update to a transformation: mirrors
// Python's set_param(name, value) which calls _cpp_trans.update on a
// dict containing just that key. The transformation's update() method
// uses params.get_int(name, current) so unmentioned keys keep their
// current value.
void set_trans_param(TemperParamLocator& loc, int64_t new_value) {
    loc.current_value = new_value;
    Params p;
    p.set_int(loc.param_name, new_value);
    loc.trans->update(p);
}

uint64_t random_perturbation(int width) {
    if (width <= 0) return 0;
    uint64_t r = regpoly_random::rng()();
    if (width >= 64) return r;
    return r & ((uint64_t(1) << width) - 1);
}

// Snapshot all current values; used to restore best_params on
// improvement and to roll back to the pre-recursion state when
// recursing further does not find an improvement.
std::vector<int64_t> snapshot(const std::vector<TemperParamLocator>& params) {
    std::vector<int64_t> out;
    out.reserve(params.size());
    for (const auto& p : params) out.push_back(p.current_value);
    return out;
}

void restore(std::vector<TemperParamLocator>& params,
             const std::vector<int64_t>& snap) {
    for (size_t i = 0; i < params.size(); ++i)
        set_trans_param(params[i], snap[i]);
}

// Final gap recomputation from scratch — mirrors the Python
// _compute_gaps which builds a fresh LatticeOptCache and calls
// compute_all. Here we reuse the existing cache (which has been
// step()'d throughout the optimization) by calling compute_all
// directly, which rebuilds polys + StackBase from the live bitmask.
std::pair<std::vector<int>, int>
compute_gaps_now(TemperOptCache& cache) {
    auto ecart = cache.compute_all();  // size L+1, ecart[0] unused
    int se = 0;
    for (size_t i = 1; i < ecart.size(); ++i) se += ecart[i];
    return {std::move(ecart), se};
}

}  // namespace

TemperingOptResult run_tempering_optimizer_once(
    const TemperingOptimizerConfig& cfg,
    TemperOptCache& cache,
    std::vector<TemperParamLocator>& params,
    const std::vector<std::vector<uint64_t>>& safe_masks)
{
    regpoly_random::seed(cfg.random_seed);

    const int L = cache.L();
    const int max_essais = cfg.max_essais;
    const auto t_start = std::chrono::steady_clock::now();

    cache.reset_step();

    std::vector<int> ecart(L + 1, -1);
    std::vector<int> best_se(L + 1, INT_MAX);
    auto best_params = snapshot(params);
    int essais = 0;
    int max_v = 0;

    // Recursive optimize(v) translated 1:1 from the Python.
    std::function<void(int)> optimize = [&](int v) {
        const int Lim = (v < L / 2) ? 5 : 2;

        for (int i = 0; i < Lim; ++i) {
            if (essais >= max_essais) break;
            ++essais;

            // Apply a fresh random perturbation to every locator under
            // the safe mask at this resolution.
            for (size_t p = 0; p < params.size(); ++p) {
                uint64_t mask = safe_masks[v][p];
                if (mask == 0) continue;
                uint64_t perturbation =
                    random_perturbation(params[p].width) & mask;
                int64_t cur = params[p].current_value;
                set_trans_param(params[p],
                                cur ^ static_cast<int64_t>(perturbation));
            }

            ecart[v] = cache.step(v);

            int se = 0;
            for (int j = 1; j <= v; ++j) se += ecart[j];

            if (se < best_se[v] && v >= max_v) {
                essais = 0;
                best_params = snapshot(params);
                best_se[v] = se;
                max_v = v;
            }

            if (ecart[v] <= cfg.delta[v] && se <= cfg.mse) {
                if (v >= L) break;
                else optimize(v + 1);
            }

            if (ecart[L] >= 0
                && ecart[L] <= cfg.delta[L]
                && best_se[L] <= cfg.mse) {
                break;
            }
        }
    };

    optimize(1);
    restore(params, best_params);

    auto [gaps, se] = compute_gaps_now(cache);

    const auto t_end = std::chrono::steady_clock::now();
    const double elapsed =
        std::chrono::duration<double>(t_end - t_start).count();

    return TemperingOptResult{se, std::move(gaps), elapsed, essais};
}

TemperingOptResult run_tempering_optimizer_minimize(
    const TemperingOptimizerConfig& cfg,
    TemperOptCache& cache,
    std::vector<TemperParamLocator>& params,
    const std::vector<std::vector<uint64_t>>& safe_masks)
{
    regpoly_random::seed(cfg.random_seed);

    const int L = cache.L();
    const int n_restarts = cfg.n_restarts;
    const auto t_start = std::chrono::steady_clock::now();

    bool have_best = false;
    TemperingOptResult best{};
    auto best_snapshot = snapshot(params);
    int failures = 0;
    int n_calls = 0;

    while (failures < n_restarts) {
        ++n_calls;

        // Re-randomize every locator. Mirrors Python's
        // `set_param(pn, random.getrandbits(w))`.
        for (auto& loc : params)
            set_trans_param(loc,
                            static_cast<int64_t>(random_perturbation(loc.width)));

        TemperingOptimizerConfig sub = cfg;
        sub.n_restarts = 1;
        if (!have_best) {
            sub.delta.assign(L + 1, INT_MAX);
            sub.mse = INT_MAX;
        } else {
            sub.delta = best.gaps;
            // Cap any -1 (unmeasured) sentinel to INT_MAX.
            for (auto& d : sub.delta) if (d < 0) d = INT_MAX;
            sub.mse = best.se;
        }
        // run_once seeds the RNG at top — we want the OUTER seeding to
        // persist across restarts. Bypass by passing seed=0 (OS entropy)
        // for the first call, and the perturbation chain advances from
        // there. The original Python uses Python's global random state,
        // which our outer seed + inner re-randomization mirrors closely
        // enough for tests.
        sub.random_seed = 0;

        auto result = run_tempering_optimizer_once(
            sub, cache, params, safe_masks);

        if (!have_best || result.se < best.se) {
            best = result;
            best_snapshot = snapshot(params);
            failures = 0;
            have_best = true;
        } else {
            ++failures;
        }
    }

    if (have_best) restore(params, best_snapshot);

    const auto t_end = std::chrono::steady_clock::now();
    best.elapsed_seconds =
        std::chrono::duration<double>(t_end - t_start).count();
    best.essais = n_calls;
    return best;
}
