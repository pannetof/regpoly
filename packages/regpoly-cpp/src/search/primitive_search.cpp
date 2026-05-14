// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "primitive_search.h"

#include "factory.h"
#include "generator_registry.h"
#include "param_spec.h"
#include "params.h"
#include "primitivity.h"
#include "random_samplers.h"
#include "well.h"

#include <chrono>
#include <memory>
#include <stdexcept>

namespace {

// Copy every key from `src` into `dst`, overwriting existing values.
void merge_into(Params& dst, const Params& src) {
    for (const auto& kv : src.ints())        dst.set_int(kv.first, kv.second);
    for (const auto& kv : src.bools())       dst.set_bool(kv.first, kv.second);
    for (const auto& kv : src.strings())     dst.set_string(kv.first, kv.second);
    for (const auto& kv : src.int_vecs())    dst.set_int_vec(kv.first, kv.second);
    for (const auto& kv : src.uint_vecs())   dst.set_uint_vec(kv.first, kv.second);
    for (const auto& kv : src.struct_maps()) dst.set_struct_map(kv.first, kv.second);
}

}  // namespace

int64_t run_primitive_search(
    const PrimitiveSearchConfig& cfg,
    const OnHitFn& on_hit,
    const OnProgressFn& on_progress)
{
    regpoly_random::seed(cfg.random_seed);

    auto specs = get_gen_param_specs(cfg.family);
    const int progress_every = cfg.progress_interval > 0
                                   ? cfg.progress_interval
                                   : 100;

    // Cost-cap setup: only WELL families honour max_cost. Resolve
    // family aliases (WELLRNG, Carry2Gen → WELLGen) so callers using
    // legacy names trigger the same path.
    const std::string canonical_family = [&]() {
        const auto* info = GeneratorRegistry::find(cfg.family);
        return info ? info->canonical_name : cfg.family;
    }();
    const bool cost_capped =
        (cfg.max_cost > 0 && canonical_family == "WELLGen");
    if (cost_capped && cfg.fixed_params.has_struct_map("matrices")) {
        throw std::invalid_argument(
            "PrimitiveSearch: WELL `matrices` pin and `max_cost` are "
            "mutually exclusive — pinning `matrices` fixes the cost; "
            "max_cost only applies when the search varies matrices.");
    }
    const int well_w = cost_capped
        ? static_cast<int>(cfg.structural_params.get_int("w", 32))
        : 32;

    const auto t_start = std::chrono::steady_clock::now();
    int64_t tries = 0;

    while (true) {
        if (cfg.max_tries > 0 && tries >= cfg.max_tries) break;

        const auto now = std::chrono::steady_clock::now();
        const double elapsed =
            std::chrono::duration<double>(now - t_start).count();
        if (cfg.max_seconds > 0.0 && elapsed >= cfg.max_seconds) break;

        ++tries;

        Params p;
        // Cost-bounded WELL: pre-populate the per-iteration scratch
        // Params with a freshly-sampled cost-bounded `matrices` BEFORE
        // build_iteration_params runs, so the random-sampler dispatch
        // skips matrices (already present) instead of complaining that
        // it has no rand_type. The pre-loop check has already
        // guaranteed the user did not pin matrices in fixed_params.
        if (cost_capped) {
            p.set_struct_map(
                "matrices",
                WELLGen::random_matrices(well_w, cfg.max_cost,
                                         regpoly_random::rng()));
        }
        // Merge structural + fixed + sampled-where-needed.
        merge_into(p, cfg.structural_params);
        merge_into(p, cfg.fixed_params);
        bool sample_failed = false;
        for (const auto& spec : specs) {
            if (p.has(spec.name)) continue;
            if (spec.has_default) continue;
            if (spec.structural)
                throw std::invalid_argument(
                    "PrimitiveSearch: structural parameter '" + spec.name +
                    "' must be provided in structural_params");
            const std::string& rt = spec.rand_type;
            if (rt.empty() || rt == "none")
                throw std::invalid_argument(
                    "PrimitiveSearch: parameter '" + spec.name +
                    "' has no default and no random sampler — fixed_params "
                    "must pin it");
            try {
                if (!regpoly_random::sample_param_into(spec, p, cfg.L))
                    throw std::invalid_argument(
                        "PrimitiveSearch: no sampler registered for rand_type '"
                        + rt + "' (parameter '" + spec.name + "')");
            } catch (const std::invalid_argument&) {
                // A sampler can fail when its rand_args reference another
                // just-sampled param whose value makes the range
                // degenerate (e.g. q's "1,t-1" when t==1). Abandon this
                // candidate and start over rather than propagating the
                // exception out of the search loop.
                sample_failed = true;
                break;
            }
        }
        if (sample_failed) continue;

        std::unique_ptr<Generator> gen;
        try {
            gen = create_generator(cfg.family, p, cfg.L);
        } catch (const std::exception&) {
            // Bad random parameter combo (e.g. inadmissible Tausworthe
            // poly + s pairing); skip and continue. Mirrors
            // PrimitiveSearch.run()'s try/except around Generator.create.
            continue;
        }

        if (is_full_period(*gen)) {
            TestedGenerator hit{
                cfg.family, gen->k(), cfg.L, std::move(p), tries};
            on_hit(hit);
        }

        if (tries % progress_every == 0) {
            const auto t_now = std::chrono::steady_clock::now();
            const double e =
                std::chrono::duration<double>(t_now - t_start).count();
            on_progress(SearchProgress{tries, e});
        }
    }

    // Final progress event so the caller always sees the closing tally.
    const auto t_end = std::chrono::steady_clock::now();
    const double elapsed_total =
        std::chrono::duration<double>(t_end - t_start).count();
    on_progress(SearchProgress{tries, elapsed_total});

    return tries;
}
