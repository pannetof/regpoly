// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "generator.h"

#include <functional>
#include <memory>
#include <string>
#include <vector>

/**
 * @file equidistribution_method.h
 * @brief Polymorphic equidistribution-method interface + registry.
 * @ingroup core
 *
 * Polymorphic equidistribution method. Replaces the `SeekTestKind` /
 * if-chain dispatch that lived in
 * `seek_search.cpp::run_equidist_test`. Each concrete method
 * (`matricial`, `lattice`, `harase`, `notprimitive`,
 * `simd_notprimitive`, `nothing`) is a stateless object that wraps
 * the corresponding free-function kernel and reports its result
 * through a uniform value type.
 *
 * Adding a new method `BazMethod` means:
 *  1. Write a `BazMethod : public EquidistributionMethod` subclass
 *     that calls the kernel free function in `run()`.
 *  2. Register it via `MethodRegistry::reg("baz", ...)` in
 *     `factory.cpp`.
 *
 * No existing dispatch site needs editing.
 */

namespace regpoly::core {

/**
 * @brief Uniform result type returned by every `EquidistributionMethod::run`.
 *
 * Carries the per-resolution gap vector `ecart`, the cumulative SE,
 * and the `verified` flag callers must consult before treating `se`
 * as authoritative. `verified` is unconditionally true for the
 * `lattice` / `harase` / `notprimitive` / `simd_notprimitive`
 * methods; `matricial` derives it from the kernel's own state
 * machine, and `nothing` always returns false.
 *
 * @ingroup core
 */
struct EquidistributionMethodResult {
    std::vector<int> ecart;  ///< Per-resolution gap; size `maxL + 1`; `ecart[0]` unused.
    int  se = 0;             ///< Cumulative equidistribution gap.
    bool verified = false;   ///< True iff the kernel certified its result.
};

/**
 * @brief Abstract interface implemented by every equidistribution method.
 *
 * Concrete subclasses wrap one of the algorithmic kernels exposed
 * under the `lattice/` headers. Each method is stateless and
 * constructed fresh per dispatch via `MethodRegistry`.
 *
 * @ingroup core
 */
class EquidistributionMethod {
public:
    virtual ~EquidistributionMethod() = default;

    /**
     * @brief Stable method identifier used by the registry.
     * @return  Method name (`"matricial"`, `"lattice"`, etc.).
     */
    virtual std::string name() const = 0;

    /**
     * @brief Run the wrapped kernel on `gen`.
     *
     * @param gen   Generator (typically a `CombinedGenerator`) to evaluate.
     * @param kg    Combined state size.
     * @param L     Output word width.
     * @param maxL  Maximum resolution to test.
     * @param delta Per-resolution gap budget (size `maxL + 1`).
     * @param mse   Upper bound on the cumulative gap.
     * @return      Uniform `EquidistributionMethodResult`.
     */
    virtual EquidistributionMethodResult run(
        const Generator& gen,
        int kg, int L, int maxL,
        const std::vector<int>& delta,
        int mse) const = 0;
};


/**
 * @brief Name-string registry for equidistribution methods.
 *
 * Maps method names (`"matricial"`, `"lattice"`, `"harase"`,
 * `"notprimitive"`, `"simd_notprimitive"`, `"nothing"`) to factory
 * closures. Registration happens at static-init time in
 * `factory.cpp`; lookup builds a fresh instance per call (cheap —
 * methods are stateless). Mirrors the pattern used by
 * `GeneratorRegistry` in `core/`.
 *
 * @ingroup core
 */
class MethodRegistry {
public:
    /// Factory closure type stored under each registered name.
    using FactoryFn = std::function<std::unique_ptr<EquidistributionMethod>()>;

    /**
     * @brief Register a factory under `name`.
     *
     * Returns a dummy int so the call can be used as the initialiser
     * of a static variable (standard "register-on-static-init" idiom).
     *
     * @param name     Method name.
     * @param factory  Factory closure producing fresh instances.
     * @return         An ignorable int (always `0`).
     */
    static int reg(const std::string& name, FactoryFn factory);

    /**
     * @brief Build a fresh instance for `name`.
     * @param name  Registered method name.
     * @return      Owning pointer to a new instance.
     * @throws std::invalid_argument  If `name` is not registered.
     */
    static std::unique_ptr<EquidistributionMethod> create(const std::string& name);

    /**
     * @brief Non-throwing lookup.
     * @param name  Method name to test.
     * @return      True iff the registry knows `name`.
     */
    static bool has(const std::string& name);

    /**
     * @brief Enumerate every registered method name.
     * @return  Names in registration order.
     */
    static std::vector<std::string> names();
};

}  // namespace regpoly::core
