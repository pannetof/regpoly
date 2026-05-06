// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "primitive_search.h"

#include "factory.h"
#include "param_spec.h"
#include "params.h"
#include "primitivity.h"
#include "random_samplers.h"

#include <chrono>
#include <memory>
#include <stdexcept>

namespace {

// Copy every key from `src` into `dst`, overwriting existing values.
void merge_into(Params& dst, const Params& src) {
    for (const auto& kv : src.ints())      dst.set_int(kv.first, kv.second);
    for (const auto& kv : src.bools())     dst.set_bool(kv.first, kv.second);
    for (const auto& kv : src.strings())   dst.set_string(kv.first, kv.second);
    for (const auto& kv : src.int_vecs())  dst.set_int_vec(kv.first, kv.second);
    for (const auto& kv : src.uint_vecs()) dst.set_uint_vec(kv.first, kv.second);
}

// Build the per-iteration full Params bag: structural + fixed +
// (randomized values for every spec that is not already present, has
// no default, and is not structural). Throws if a required spec has
// no random sampler and is not provided.
Params build_iteration_params(
    const PrimitiveSearchConfig& cfg,
    const std::vector<ParamSpec>& specs)
{
    Params p;
    merge_into(p, cfg.structural_params);
    merge_into(p, cfg.fixed_params);

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

        if (!regpoly_random::sample_param_into(spec, p, cfg.L))
            throw std::invalid_argument(
                "PrimitiveSearch: no sampler registered for rand_type '"
                + rt + "' (parameter '" + spec.name + "')");
    }

    return p;
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

    const auto t_start = std::chrono::steady_clock::now();
    int64_t tries = 0;

    while (true) {
        if (cfg.max_tries > 0 && tries >= cfg.max_tries) break;

        const auto now = std::chrono::steady_clock::now();
        const double elapsed =
            std::chrono::duration<double>(now - t_start).count();
        if (cfg.max_seconds > 0.0 && elapsed >= cfg.max_seconds) break;

        ++tries;

        Params p = build_iteration_params(cfg, specs);

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
