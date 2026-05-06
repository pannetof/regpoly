// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

// Single-Generator& adapter overloads for the lattice equidistribution
// kernels. Each adapter unpacks the input via the virtual
// Generator::components() / Generator::tempering_chains() (whose default
// returns {*this} / {{}} for primitives, overridden by CombinedGenerator
// to return its components and per-component tempering) and forwards
// to the existing vector overload.

#include "generator.h"
#include "transformation.h"
#include "me_helpers.h"
#include "me_harase.h"
#include "me_notprimitive.h"
#include "me_notprimitive_simd.h"
#include "temper_optimizer.h"

#include <vector>

namespace {

struct Unpacked {
    std::vector<Generator*> gens;
    std::vector<std::vector<Transformation*>> trans;
};

Unpacked unpack(const Generator& gen) {
    return {gen.components(), gen.tempering_chains()};
}

}  // namespace

MeLatResult test_me_lat(
    const Generator& gen,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse)
{
    auto u = unpack(gen);
    return test_me_lat(u.gens, u.trans, kg, L, maxL, delta, mse);
}

MeLatResult test_me_harase(
    const Generator& gen,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse)
{
    auto u = unpack(gen);
    return test_me_harase(u.gens, u.trans, kg, L, maxL, delta, mse);
}

int compute_kv(const Generator& gen, int kg, int v) {
    auto u = unpack(gen);
    return compute_kv(u.gens, u.trans, kg, v);
}

MeLatResult test_me_notprimitive(
    const Generator& gen,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse)
{
    auto u = unpack(gen);
    return test_me_notprimitive(u.gens, u.trans, kg, L, maxL, delta, mse);
}

MeLatResult test_me_notprimitive_simd(
    const Generator& gen,
    int kg, int L, int maxL,
    const std::vector<int>& delta, int mse)
{
    auto u = unpack(gen);
    return test_me_notprimitive_simd(u.gens, u.trans, kg, L, maxL, delta, mse);
}

PISCache::PISCache(const Generator& gen, int kg, int L)
    : PISCache(unpack(gen).gens, unpack(gen).trans, kg, L) {}

TemperOptCache::TemperOptCache(const Generator& gen, int kg, int L)
    : TemperOptCache(unpack(gen).gens, unpack(gen).trans, kg, L) {}
