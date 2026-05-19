// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "generator.h"
#include "transformation.h"
#include "params.h"
#include "param_spec.h"
#include "gen_enumerator.h"
#include <memory>
#include <string>
#include <vector>

/**
 * @file factory.h
 * @brief Public entry points into the registry-backed runtime factory.
 *
 * Every generator and transformation lookup goes through these four
 * free functions. They are pure thin shims over `GeneratorRegistry`
 * and `TransformationRegistry`; no per-family if-chains live here.
 * The pybind11 binding-registration hook lives at the bottom and is
 * called from `_regpoly_cpp` module init.
 *
 * @ingroup core
 */

// Forward declaration for pybind11 module type
namespace pybind11 { class module_; }

namespace regpoly::core {

/**
 * @brief Construct a concrete `Generator` of the given family from `params`.
 *
 * Looks `family` up in `GeneratorRegistry`, then dispatches to the
 * subclass's `from_params` factory. Aliases (legacy YAML / SQL
 * spellings registered via `REGISTER_GENERATOR_ALIAS`) are resolved
 * to the canonical entry.
 *
 * @param family  Canonical family name or registered alias.
 * @param params  Typed parameter bag (see `Params`).
 * @param L       Output word width in bits.
 * @return        Heap-allocated concrete generator owned by the caller.
 * @throws std::invalid_argument  If `family` is not a registered name.
 *
 * @code{.cpp}
 *   using namespace regpoly::core;
 *   Params p;
 *   p.set_int("w", 32); p.set_int("r", 624); p.set_int("m", 397);
 *   p.set_int("p", 31); p.set_int("a", 0x9908b0df);
 *   auto gen = create_generator("MTGen", p, 32);
 *   gen->init(BitVect(gen->k()));
 *   gen->next();
 * @endcode
 *
 * @see :py:func:`regpoly.core.generator.Generator.create`
 */
std::unique_ptr<Generator> create_generator(
    const std::string& family, const Params& params, int L);

/**
 * @brief Construct a concrete `Transformation` of the given type from `params`.
 *
 * Mirror of `create_generator` for the tempering chain. No SIMD /
 * per-class binding plumbing — `TransformationRegistry` carries a
 * single `from_params` callable per type.
 *
 * @param type    Registered transformation type name (e.g. `"permut"`,
 *                `"tempMK"`, `"tempMK2"`, `"laggedTempering"`).
 * @param params  Typed parameter bag.
 * @return        Heap-allocated transformation owned by the caller.
 * @throws std::invalid_argument  If `type` is not a registered name.
 *
 * @see :py:func:`regpoly.core.transformation.Transformation.create`
 */
std::unique_ptr<Transformation> create_transformation(
    const std::string& type, const Params& params);

/**
 * @brief Return the declarative parameter specs of a generator family.
 *
 * Wraps `GeneratorRegistry::lookup(family).param_specs()`. Used by
 * the YAML loader to validate parameter names / types and by the
 * random sampler to dispatch on `rand_type`.
 *
 * @param family  Canonical family name or registered alias.
 * @return        Declarative spec vector (one entry per parameter).
 * @throws std::invalid_argument  If `family` is not a registered name.
 *
 * @see :py:func:`regpoly.introspection.get_gen_param_specs`
 */
std::vector<ParamSpec> get_gen_param_specs(const std::string& family);

/**
 * @brief Return the declarative parameter specs of a transformation type.
 *
 * Wraps `TransformationRegistry::lookup(type).param_specs()`.
 *
 * @param type  Registered transformation type name.
 * @return      Declarative spec vector (one entry per parameter).
 * @throws std::invalid_argument  If `type` is not a registered name.
 *
 * @see :py:func:`regpoly.introspection.get_trans_param_specs`
 */
std::vector<ParamSpec> get_trans_param_specs(const std::string& type);

/**
 * @brief Pybind11 binding hook: walks the registry and emits one `py::class_` per family.
 *
 * Internal plumbing; called once from `_regpoly_cpp`'s `PYBIND11_MODULE`
 * block. Each generator subclass that registered a non-null bind
 * function gets a Python-side class exposing its specific surface;
 * the others remain accessible only through `create_generator`.
 *
 * @param m  Target pybind11 module.
 */
void register_generator_types(pybind11::module_& m);

}  // namespace regpoly::core
