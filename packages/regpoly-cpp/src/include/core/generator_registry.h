// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once

#include "gen_enumerator.h"
#include "generator.h"
#include "param_spec.h"
#include "params.h"

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

/**
 * @file generator_registry.h
 * @brief Self-registering registry of every concrete `Generator` subclass.
 *
 * Defines `GeneratorRegistry` plus three macros
 * (`REGISTER_GENERATOR`, `REGISTER_GENERATOR_ENUMERABLE`,
 * `REGISTER_GENERATOR_ALIAS`) that subclasses use to register
 * themselves at translation-unit init time. The factory entry
 * points in `factory.h` are pure registry lookups — adding a new
 * family touches no central file beyond writing the subclass and
 * its single `REGISTER_*` line.
 *
 * @ingroup core
 */

// Forward declaration for pybind11 module (avoid pulling pybind11 into
// regpoly_core's public surface).
namespace pybind11 { class module_; }

namespace regpoly::core {

/**
 * @brief Self-registering registry of every concrete `Generator` subclass.
 *
 * Each subclass (e.g. `FooGen`) registers itself once per build by
 * writing a single line at the bottom of `foogen.cpp`:
 *
 * @code{.cpp}
 *   REGISTER_GENERATOR(FooGen, "FooGen", &bind_foogen);
 * @endcode
 *
 * where `bind_foogen` is a small function that instantiates a
 * `py::class_<FooGen, Generator>(m, "FooGen")` (only needed when the
 * pybind11 extension is built; otherwise pass `nullptr`).
 *
 * Aliases (legacy YAML / SQL spellings) are registered alongside:
 *
 * @code{.cpp}
 *   REGISTER_GENERATOR_ALIAS("FooGen", "OldFooName");
 * @endcode
 *
 * The factory entry points (`create_generator`,
 * `get_gen_param_specs`, `register_generator_types`,
 * `family_is_enumerable`, `make_gen_enumerator`) look up the
 * registry instead of branching on string literals; adding a new
 * family touches no central file.
 *
 * @ingroup core
 */
class GeneratorRegistry {
public:
    /** @brief Factory callable: build a `Generator` from a `Params` bag and width `L`. */
    using FromParamsFn = std::function<
        std::unique_ptr<Generator>(const Params& params, int L)>;
    /** @brief Callable returning the subclass's declarative parameter specs. */
    using ParamSpecsFn = std::function<std::vector<ParamSpec>()>;
    /** @brief Pybind11 binder: emits the per-class `py::class_` if Python is built. */
    using BindFn = std::function<void(pybind11::module_&)>;
    /** @brief Optional enumerator factory; empty when the family is non-enumerable. */
    using EnumeratorFn = std::function<
        std::unique_ptr<GenEnumerator>(const Params& resolved, int L)>;

    /**
     * @brief Per-family registry record.
     * @ingroup core
     */
    struct Info {
        std::string canonical_name;    ///< Registered (canonical) family name.
        FromParamsFn from_params;      ///< Subclass `from_params` factory.
        ParamSpecsFn param_specs;      ///< Subclass `param_specs` callable.
        BindFn bind;                   ///< Pybind11 binder (may be empty).
        EnumeratorFn make_enumerator;  ///< Enumerator factory (empty → non-enumerable).
    };

    /**
     * @brief Register a concrete generator family.
     *
     * Idempotent under the same canonical name — calling twice with
     * the same name and identical factories is a no-op (defensive
     * against translation-unit duplicate-link rare cases).
     *
     * @param canonical_name   Registered family name.
     * @param from_params      Subclass `from_params` factory.
     * @param param_specs      Subclass `param_specs` callable.
     * @param bind             Optional pybind11 binder (empty allowed).
     * @param make_enumerator  Optional enumerator factory (empty allowed).
     * @return                 Sentinel int (always 0) for use in a static-init expression.
     */
    static int reg(const std::string& canonical_name,
                   FromParamsFn from_params,
                   ParamSpecsFn param_specs,
                   BindFn bind = {},
                   EnumeratorFn make_enumerator = {});

    /**
     * @brief Map an alternate spelling onto a canonical entry.
     *
     * Caller is responsible for the canonical entry being
     * registered first (this enforces the convention that
     * `REGISTER_GENERATOR_ALIAS` lines appear after the matching
     * `REGISTER_GENERATOR` line in `factory.cpp`).
     *
     * @param alias           Alternate spelling.
     * @param canonical_name  Existing canonical family name.
     * @return                Sentinel int (always 0).
     * @throws std::logic_error  If `canonical_name` is not yet registered.
     */
    static int reg_alias(const std::string& alias,
                         const std::string& canonical_name);

    /**
     * @brief Resolve a family name (canonical or alias) to its registry record.
     *
     * @param name  Family name (canonical or alias).
     * @return      Const reference to the registered `Info`.
     * @throws std::invalid_argument  If `name` is not registered.
     */
    static const Info& lookup(const std::string& name);

    /**
     * @brief Resolve a family name (canonical or alias) to its registry record (non-throwing).
     *
     * @param name  Family name (canonical or alias).
     * @return      Pointer to the registered `Info`, or `nullptr` if unknown.
     */
    static const Info* find(const std::string& name);

    /**
     * @brief All canonical entries in registration order. Aliases excluded.
     * @return  Const reference to the registration-order vector.
     */
    static const std::vector<Info*>& canonical_entries();

    /**
     * @brief Vector of `(alias, canonical_name)` pairs in registration order.
     *
     * Used by `factory_bindings.cpp` to mirror legacy spellings as
     * Python module attributes.
     *
     * @return  Const reference to the alias-pair vector.
     */
    static const std::vector<std::pair<std::string, std::string>>& aliases();
};

/**
 * @brief Self-register a `Generator` subclass at translation-unit init time.
 *
 * Expands to a static-initialiser line at file scope. Place once at
 * the bottom of each `*gen.cpp`.
 *
 * @param Cls         Subclass type name.
 * @param Name        Canonical factory name (matches family strings in YAML).
 * @param BindFnPtr   Optional pybind11 binder; pass `nullptr` if the
 *                    Python extension is not being built.
 */
#define REGISTER_GENERATOR(Cls, Name, BindFnPtr)                  \
    namespace {                                                   \
    [[maybe_unused]] static const int _regpoly_gen_reg_##Cls =    \
        GeneratorRegistry::reg(                                   \
            (Name),                                               \
            &Cls::from_params,                                    \
            &Cls::param_specs,                                    \
            (BindFnPtr));                                         \
    }

/**
 * @brief Self-register a `Generator` subclass together with an enumerator factory.
 *
 * Same as `REGISTER_GENERATOR` plus a fourth argument that supplies
 * the exhaustive `GenEnumerator` factory. `TauswortheGen` is the
 * canonical example.
 *
 * @param Cls         Subclass type name.
 * @param Name        Canonical factory name.
 * @param BindFnPtr   Optional pybind11 binder; pass `nullptr` to skip.
 * @param EnumFnPtr   `make_enumerator`-style factory pointer.
 */
#define REGISTER_GENERATOR_ENUMERABLE(Cls, Name, BindFnPtr, EnumFnPtr) \
    namespace {                                                       \
    [[maybe_unused]] static const int _regpoly_gen_reg_##Cls =        \
        GeneratorRegistry::reg(                                       \
            (Name),                                                   \
            &Cls::from_params,                                        \
            &Cls::param_specs,                                        \
            (BindFnPtr),                                              \
            (EnumFnPtr));                                             \
    }

/**
 * @brief Map an alternate (legacy) spelling onto a canonical generator name.
 *
 * Place near the corresponding `REGISTER_GENERATOR` line. The
 * canonical name must already be registered when this macro
 * expansion runs.
 *
 * @param CanonicalName  Existing canonical family name.
 * @param Alias          Alternate spelling to register.
 */
#define REGISTER_GENERATOR_ALIAS(CanonicalName, Alias)                \
    namespace {                                                       \
    [[maybe_unused]] static const int                                 \
        _regpoly_gen_alias_##__LINE__ =                               \
            GeneratorRegistry::reg_alias((Alias), (CanonicalName));   \
    }

}  // namespace regpoly::core
