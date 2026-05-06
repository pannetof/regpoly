// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "equidistribution_method.h"

#include "equidistribution_runner.h"  // run_matricial_equidistribution
#include "me_helpers.h"               // test_me_lat
#include "me_harase.h"                // test_me_harase
#include "me_notprimitive.h"          // test_me_notprimitive
#include "me_notprimitive_simd.h"     // test_me_notprimitive_simd

#include <deque>
#include <stdexcept>
#include <unordered_map>
#include <utility>

// ── Concrete method classes ─────────────────────────────────────────────

namespace {

class MatricialMethod : public EquidistributionMethod {
public:
    std::string name() const override { return "matricial"; }
    EquidistributionMethodResult run(
        const Generator& gen, int kg, int L, int maxL,
        const std::vector<int>& delta, int mse) const override
    {
        auto r = run_matricial_equidistribution(
            gen, kg, L, maxL, delta, mse);
        return {std::move(r.ecart), r.se, r.verified};
    }
};

// Helper: every lattice-family method shares the post-processing of an
// MeLatResult — verified is always true on return and se is the kernel
// result.
EquidistributionMethodResult wrap_lat(MeLatResult r) {
    return {std::move(r.ecart), r.se, /*verified=*/true};
}

class LatticeMethod : public EquidistributionMethod {
public:
    std::string name() const override { return "lattice"; }
    EquidistributionMethodResult run(
        const Generator& gen, int kg, int L, int maxL,
        const std::vector<int>& delta, int mse) const override
    {
        return wrap_lat(test_me_lat(gen, kg, L, maxL, delta, mse));
    }
};

class HaraseMethod : public EquidistributionMethod {
public:
    std::string name() const override { return "harase"; }
    EquidistributionMethodResult run(
        const Generator& gen, int kg, int L, int maxL,
        const std::vector<int>& delta, int mse) const override
    {
        return wrap_lat(test_me_harase(gen, kg, L, maxL, delta, mse));
    }
};

class NotPrimitiveMethod : public EquidistributionMethod {
public:
    std::string name() const override { return "notprimitive"; }
    EquidistributionMethodResult run(
        const Generator& gen, int kg, int L, int maxL,
        const std::vector<int>& delta, int mse) const override
    {
        return wrap_lat(test_me_notprimitive(gen, kg, L, maxL, delta, mse));
    }
};

class SimdNotPrimitiveMethod : public EquidistributionMethod {
public:
    std::string name() const override { return "simd_notprimitive"; }
    EquidistributionMethodResult run(
        const Generator& gen, int kg, int L, int maxL,
        const std::vector<int>& delta, int mse) const override
    {
        return wrap_lat(test_me_notprimitive_simd(gen, kg, L, maxL, delta, mse));
    }
};

// "nothing" — the disabled-test sentinel. Returns zeros and verified=false.
class NothingMethod : public EquidistributionMethod {
public:
    std::string name() const override { return "nothing"; }
    EquidistributionMethodResult run(
        const Generator&, int, int, int maxL,
        const std::vector<int>&, int) const override
    {
        return {std::vector<int>(maxL + 1, 0), /*se=*/0, /*verified=*/false};
    }
};

// ── Registry storage ───────────────────────────────────────────────────

struct Slot {
    std::string name;
    MethodRegistry::FactoryFn factory;
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
        MethodRegistry::reg("matricial",
            []{ return std::unique_ptr<EquidistributionMethod>(new MatricialMethod); });
        MethodRegistry::reg("lattice",
            []{ return std::unique_ptr<EquidistributionMethod>(new LatticeMethod); });
        MethodRegistry::reg("harase",
            []{ return std::unique_ptr<EquidistributionMethod>(new HaraseMethod); });
        MethodRegistry::reg("notprimitive",
            []{ return std::unique_ptr<EquidistributionMethod>(new NotPrimitiveMethod); });
        MethodRegistry::reg("simd_notprimitive",
            []{ return std::unique_ptr<EquidistributionMethod>(new SimdNotPrimitiveMethod); });
        MethodRegistry::reg("nothing",
            []{ return std::unique_ptr<EquidistributionMethod>(new NothingMethod); });
        return 0;
    }();
    (void)once;
}

}  // namespace

// ── Public registry API ─────────────────────────────────────────────────

int MethodRegistry::reg(const std::string& name, FactoryFn factory) {
    auto& m = by_name();
    if (m.count(name)) return 0;
    auto& slot = storage().emplace_back(Slot{name, std::move(factory)});
    m.emplace(name, &slot);
    order().push_back(name);
    return 0;
}

std::unique_ptr<EquidistributionMethod>
MethodRegistry::create(const std::string& name) {
    register_builtin_methods();
    auto& m = by_name();
    auto it = m.find(name);
    if (it == m.end()) {
        throw std::invalid_argument("Unknown equidistribution method: " + name);
    }
    return it->second->factory();
}

bool MethodRegistry::has(const std::string& name) {
    register_builtin_methods();
    return by_name().count(name) != 0;
}

std::vector<std::string> MethodRegistry::names() {
    register_builtin_methods();
    return order();
}
