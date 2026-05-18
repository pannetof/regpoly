// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "generator.h"

#include <functional>
#include <memory>
#include <string>
#include <vector>

// Polymorphic equidistribution method.
//
// Replaces the SeekTestKind / if-chain dispatch that lived in
// seek_search.cpp:run_equidist_test. Each concrete method (matricial,
// lattice, harase, notprimitive, simd_notprimitive, nothing) is a
// stateless object that wraps the corresponding free-function kernel
// and reports its result through a uniform value type.
//
// Adding a new method `BazMethod` means:
//   1. Write a `BazMethod : public EquidistributionMethod` subclass
//      that calls the kernel free function in run().
//   2. Register it via `MethodRegistry::reg("baz", ...)` in factory.cpp.
// No existing dispatch site needs editing.

namespace regpoly::core {

struct EquidistributionMethodResult {
    std::vector<int> ecart;  // size maxL+1; ecart[0] unused
    int  se = 0;
    bool verified = false;
    // For some methods (lattice, harase, notprimitive, simd_notprimitive)
    // verified is unconditionally true on return; "matricial" sets it
    // from the kernel's own state machine; "nothing" sets it false.
    // Callers MUST consult this flag before treating se as authoritative.
};

class EquidistributionMethod {
public:
    virtual ~EquidistributionMethod() = default;
    virtual std::string name() const = 0;
    virtual EquidistributionMethodResult run(
        const Generator& gen,
        int kg, int L, int maxL,
        const std::vector<int>& delta,
        int mse) const = 0;
};


// Registry — string -> factory function. The factory builds a fresh
// instance per lookup (cheap; methods are stateless).
class MethodRegistry {
public:
    using FactoryFn = std::function<std::unique_ptr<EquidistributionMethod>()>;

    static int reg(const std::string& name, FactoryFn factory);

    // Build a fresh instance, or throw std::invalid_argument if the
    // name is not registered.
    static std::unique_ptr<EquidistributionMethod> create(const std::string& name);

    // Non-throwing query: does the registry know `name`?
    static bool has(const std::string& name);

    // All registered method names, in registration order.
    static std::vector<std::string> names();
};

}  // namespace regpoly::core
