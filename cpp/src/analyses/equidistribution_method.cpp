// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "equidistribution_method.h"

#include "combined_f2_linear_source.h"
#include "equidistribution_runner.h"  // test_me_matricial
#include "me_helpers.h"               // test_me_lat
#include "me_harase.h"                // test_me_harase
#include "me_notprimitive.h"          // test_me_notprimitive
#include "me_notprimitive_simd.h"     // test_me_notprimitive_simd

#include <deque>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

using namespace regpoly::core;


// ── Concrete method classes ─────────────────────────────────────────────

namespace regpoly::core {

namespace {

// Lattice-family kernels (`test_me_lat`, `test_me_harase`,
// `test_me_notprimitive`, `test_me_notprimitive_simd`) take
// `const CombinedF2LinearSource&`. The polymorphic
// `EquidistributionMethod::run` interface takes `const ITestable&`
// so the registry can dispatch any source uniformly; this helper
// downcasts at the entry of every lattice-family method.
const CombinedF2LinearSource& as_combined(
    const ITestable& gen, const char* method_name)
{
    auto* cs = dynamic_cast<const CombinedF2LinearSource*>(&gen);
    if (!cs) {
        throw std::invalid_argument(
            std::string("EquidistributionMethod[") + method_name +
            "]: gen must be a CombinedF2LinearSource (got " +
            gen.name() + "). Wrap a single source with "
            "regpoly.make_combined(gen, Lmax=...) before invoking "
            "lattice-family equidistribution methods.");
    }
    return *cs;
}

class MatricialMethod : public EquidistributionMethod {
public:
    std::string name() const override { return "matricial"; }
    EquidistributionMethodResult run(
        const ITestable& gen, int kg, int L, int maxL,
        const std::vector<int>& delta, int mse) const override
    {
        auto r = test_me_matricial(
            gen, kg, L, maxL, delta, mse);
        return {std::move(r.ecart), r.se, r.verified};
    }
};

// Helper: post-processing of an EquidistributionResult from a
// lattice-family kernel. Now reads `verified` from the unified result
// type (Phase 2 step 7, Q6) — lattice kernels set it true by construction.
EquidistributionMethodResult wrap_lat(EquidistributionResult r) {
    return {std::move(r.ecart), r.se, r.verified};
}

class LatticeMethod : public EquidistributionMethod {
public:
    std::string name() const override { return "lattice"; }
    EquidistributionMethodResult run(
        const ITestable& gen, int kg, int L, int maxL,
        const std::vector<int>& delta, int mse) const override
    {
        return wrap_lat(test_me_lat(
            as_combined(gen, "lattice"), kg, L, maxL, delta, mse));
    }
};

class HaraseMethod : public EquidistributionMethod {
public:
    std::string name() const override { return "harase"; }
    EquidistributionMethodResult run(
        const ITestable& gen, int kg, int L, int maxL,
        const std::vector<int>& delta, int mse) const override
    {
        return wrap_lat(test_me_harase(
            as_combined(gen, "harase"), kg, L, maxL, delta, mse));
    }
};

class NotPrimitiveMethod : public EquidistributionMethod {
public:
    std::string name() const override { return "notprimitive"; }
    EquidistributionMethodResult run(
        const ITestable& gen, int kg, int L, int maxL,
        const std::vector<int>& delta, int mse) const override
    {
        return wrap_lat(test_me_notprimitive(
            as_combined(gen, "notprimitive"), kg, L, maxL, delta, mse));
    }
};

class SimdNotPrimitiveMethod : public EquidistributionMethod {
public:
    std::string name() const override { return "simd_notprimitive"; }
    EquidistributionMethodResult run(
        const ITestable& gen, int kg, int L, int maxL,
        const std::vector<int>& delta, int mse) const override
    {
        return wrap_lat(test_me_notprimitive_simd(
            as_combined(gen, "simd_notprimitive"), kg, L, maxL, delta, mse));
    }
};

// "nothing" — the disabled-test sentinel. Returns zeros and verified=false.
class NothingMethod : public EquidistributionMethod {
public:
    std::string name() const override { return "nothing"; }
    EquidistributionMethodResult run(
        const ITestable&, int, int, int maxL,
        const std::vector<int>&, int) const override
    {
        return {std::vector<int>(maxL + 1, 0), /*se=*/0, /*verified=*/false};
    }
};

// ── Registry storage ───────────────────────────────────────────────────

struct Slot {
    std::string name;
    EquidistributionMethodRegistry::FactoryFn factory;
};

std::unordered_map<std::string, Slot*>& by_name() {
    static std::unordered_map<std::string, Slot*> m;
    return m;
}

std::deque<Slot>& storage() {
    static std::deque<Slot> v;
    return v;
}

std::vector<std::string>& order() {
    static std::vector<std::string> v;
    return v;
}

// One-shot installer for the built-in methods. Mirrors the pattern in
// factory.cpp/register_all_generators — function-local static lambda
// runs exactly once.
void register_builtin_methods() {
    static const int once = []{
        EquidistributionMethodRegistry::reg("matricial",
            []{ return std::unique_ptr<EquidistributionMethod>(new MatricialMethod); });
        EquidistributionMethodRegistry::reg("lattice",
            []{ return std::unique_ptr<EquidistributionMethod>(new LatticeMethod); });
        EquidistributionMethodRegistry::reg("harase",
            []{ return std::unique_ptr<EquidistributionMethod>(new HaraseMethod); });
        EquidistributionMethodRegistry::reg("notprimitive",
            []{ return std::unique_ptr<EquidistributionMethod>(new NotPrimitiveMethod); });
        EquidistributionMethodRegistry::reg("simd_notprimitive",
            []{ return std::unique_ptr<EquidistributionMethod>(new SimdNotPrimitiveMethod); });
        EquidistributionMethodRegistry::reg("nothing",
            []{ return std::unique_ptr<EquidistributionMethod>(new NothingMethod); });
        return 0;
    }();
    (void)once;
}

}  // namespace

// ── Public registry API ─────────────────────────────────────────────────

int EquidistributionMethodRegistry::reg(const std::string& name, FactoryFn factory) {
    auto& m = by_name();
    if (m.count(name)) return 0;
    auto& slot = storage().emplace_back(Slot{name, std::move(factory)});
    m.emplace(name, &slot);
    order().push_back(name);
    return 0;
}

std::unique_ptr<EquidistributionMethod>
EquidistributionMethodRegistry::create(const std::string& name) {
    register_builtin_methods();
    auto& m = by_name();
    auto it = m.find(name);
    if (it == m.end()) {
        throw std::invalid_argument("Unknown equidistribution method: " + name);
    }
    return it->second->factory();
}

bool EquidistributionMethodRegistry::has(const std::string& name) {
    register_builtin_methods();
    return by_name().count(name) != 0;
}

std::vector<std::string> EquidistributionMethodRegistry::names() {
    register_builtin_methods();
    return order();
}

}  // namespace regpoly::core
