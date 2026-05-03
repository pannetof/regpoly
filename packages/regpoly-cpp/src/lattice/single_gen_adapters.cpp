// Phase 1 — single-Generator& adapter overloads for the lattice
// equidistribution kernels.
//
// Each adapter detects whether the input is a CombinedGenerator. If so,
// it unpacks the components and per-component tempering chains and
// forwards to the existing vector overload, which already implements
// the full algorithm. If the input is a primitive (non-Combined)
// Generator, it is treated as a 1-component combination with no
// tempering.
//
// This preserves exact behaviour of the existing pytest cross-checks
// while exposing the new public API (single Generator&) that Phase 2+
// migrates to consume natively.

#include "combined.h"
#include "generator.h"
#include "transformation.h"
#include "me_helpers.h"
#include "me_harase.h"
#include "me_notprimitive.h"
#include "me_notprimitive_simd.h"
#include "temper_optimizer.h"

#include <vector>

namespace {

// Build (gens, trans) pair from a Generator&.
//   - CombinedGenerator: unpack components and tempering chains.
//   - Primitive Generator: 1-component, no tempering.
struct Unpacked {
    std::vector<Generator*> gens;
    std::vector<std::vector<Transformation*>> trans;
};

Unpacked unpack(const Generator& gen) {
    if (auto* cg = dynamic_cast<const CombinedGenerator*>(&gen)) {
        return {cg->raw_component_pointers(), cg->raw_tempering_pointers()};
    }
    return {
        std::vector<Generator*>{const_cast<Generator*>(&gen)},
        std::vector<std::vector<Transformation*>>{{}},
    };
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
