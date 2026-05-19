// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "param_spec.h"
#include "params.h"
#include "transformation.h"

#include <deque>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

/**
 * @file transformation_registry.h
 * @brief Registry of `Transformation` subclasses, mirror of `GeneratorRegistry`.
 *
 * Smaller mirror of `GeneratorRegistry`: stores one `Info` per
 * transformation type (canonical name + `from_params` + `param_specs`),
 * with no per-class pybind11 binder — transformations are exposed
 * via a generic `create_transformation()` / `get_trans_param_specs()`
 * surface, so no `py::class_<Bar, Transformation>` is needed.
 *
 * @ingroup core
 */

namespace regpoly::core {

/**
 * @brief Registry of `Transformation` subclasses, mirror of `GeneratorRegistry`.
 *
 * Adding a new transformation type means:
 *
 *  1. Write `bar.h` / `bar.cpp` with the `Transformation` subclass.
 *  2. Add ONE line in `factory.cpp`'s `register_all_transformations`
 *     block:
 *     `TR::reg("bar", &BarTrans::from_params, &BarTrans::param_specs);`
 *  3. (No pybind11 binding needed — `Transformation`s are exposed
 *     via a generic `create_transformation()` / `get_trans_param_specs()`
 *     surface, so no per-class `py::class_<Bar, Transformation>`
 *     registration is required.)
 *
 * The factory entry points `create_transformation()` and
 * `get_trans_param_specs()` are pure registry lookups.
 *
 * @ingroup core
 */
class TransformationRegistry {
public:
    /** @brief Factory callable: build a `Transformation` from a `Params` bag. */
    using FromParamsFn = std::function<
        std::unique_ptr<Transformation>(const Params& params)>;
    /** @brief Callable returning the subclass's declarative parameter specs. */
    using ParamSpecsFn = std::function<std::vector<ParamSpec>()>;

    /**
     * @brief Per-type registry record.
     * @ingroup core
     */
    struct Info {
        std::string canonical_name;  ///< Registered (canonical) type name.
        FromParamsFn from_params;    ///< Subclass `from_params` factory.
        ParamSpecsFn param_specs;    ///< Subclass `param_specs` callable.
    };

    /**
     * @brief Register a concrete transformation type.
     *
     * Idempotent under the same canonical name — calling twice with
     * identical factories is a defensive no-op (guards against duplicate
     * translation-unit linkage in rare build configurations).
     *
     * @param canonical_name  Type name as it appears in YAML / SQL.
     * @param from_params     Factory callable taking a `Params` bag.
     * @param param_specs     Callable returning the declarative specs.
     * @return                Sentinel int (always 0) for use in a static-init expression.
     */
    static int reg(const std::string& canonical_name,
                   FromParamsFn from_params,
                   ParamSpecsFn param_specs);

    /**
     * @brief Resolve a transformation type name to its registry record (throwing).
     *
     * @param name  Canonical type name.
     * @return      Const reference to the registered `Info`.
     * @throws std::invalid_argument  If `name` is not registered.
     */
    static const Info& lookup(const std::string& name);

    /**
     * @brief Resolve a transformation type name to its registry record (non-throwing).
     *
     * @param name  Canonical type name.
     * @return      Pointer to the registered `Info`, or `nullptr` if unknown.
     */
    static const Info* find(const std::string& name);
};

}  // namespace regpoly::core
