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

// Forward declaration for pybind11 module (avoid pulling pybind11 into
// regpoly_core's public surface).
namespace pybind11 { class module_; }

namespace regpoly::core {

// Self-registering registry of Generator subclass types.
//
// Each subclass (FooGen) registers itself once per build by writing a
// single line at the bottom of foogen.cpp:
//
//     REGISTER_GENERATOR(FooGen, "FooGen", &bind_foogen);
//
// where `bind_foogen` is a small function that instantiates a
// `py::class_<FooGen, Generator>(m, "FooGen")` (only needed when the
// pybind11 extension is built; otherwise pass nullptr).
//
// Aliases (legacy YAML / SQL spellings) are registered alongside:
//
//     REGISTER_GENERATOR_ALIAS("FooGen", "OldFooName");
//
// The factory entry points (create_generator, get_gen_param_specs,
// register_generator_types, family_is_enumerable, make_gen_enumerator)
// look up the registry instead of branching on string literals; adding
// a new family touches no central file.
class GeneratorRegistry {
public:
    using FromParamsFn = std::function<
        std::unique_ptr<Generator>(const Params& params, int L)>;
    using ParamSpecsFn = std::function<std::vector<ParamSpec>()>;
    using BindFn = std::function<void(pybind11::module_&)>;
    using EnumeratorFn = std::function<
        std::unique_ptr<GenEnumerator>(const Params& resolved, int L)>;

    struct Info {
        std::string canonical_name;
        FromParamsFn from_params;
        ParamSpecsFn param_specs;
        BindFn bind;                  // may be empty
        EnumeratorFn make_enumerator; // may be empty (non-enumerable)
    };

    // Register a concrete generator type. Idempotent under the same
    // canonical name — calling twice with the same name and identical
    // factories is a no-op (defensive against translation-unit
    // duplicate-link rare cases). Returns a sentinel int for use in
    // a static-initialiser expression.
    static int reg(const std::string& canonical_name,
                   FromParamsFn from_params,
                   ParamSpecsFn param_specs,
                   BindFn bind = {},
                   EnumeratorFn make_enumerator = {});

    // Map an alternate spelling onto a canonical entry. Caller is
    // responsible for the canonical entry being registered first.
    static int reg_alias(const std::string& alias,
                         const std::string& canonical_name);

    // Lookup; throws std::invalid_argument if name is unknown.
    static const Info& lookup(const std::string& name);

    // Returns nullptr if unknown (non-throwing variant).
    static const Info* find(const std::string& name);

    // All canonical entries in registration order. Aliases excluded.
    static const std::vector<Info*>& canonical_entries();

    // Vector of (alias, canonical_name) pairs, in registration order.
    // Used by factory_bindings.cpp to mirror legacy spellings as
    // Python module attributes.
    static const std::vector<std::pair<std::string, std::string>>& aliases();
};

// Register a Generator subclass with the central registry. Place at
// file scope (e.g. at the bottom of the subclass's .cpp file). The
// macro expands to a static-storage int initialiser whose evaluation
// runs the registration during dynamic init of the translation unit.
#define REGISTER_GENERATOR(Cls, Name, BindFnPtr)                  \
    namespace {                                                   \
    [[maybe_unused]] static const int _regpoly_gen_reg_##Cls =    \
        GeneratorRegistry::reg(                                   \
            (Name),                                               \
            &Cls::from_params,                                    \
            &Cls::param_specs,                                    \
            (BindFnPtr));                                         \
    }

// Same with an enumerator function (TauswortheGen is the only current
// example).
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

// Map an alternate (legacy) spelling onto a canonical generator name.
// Place near the corresponding REGISTER_GENERATOR.
#define REGISTER_GENERATOR_ALIAS(CanonicalName, Alias)                \
    namespace {                                                       \
    [[maybe_unused]] static const int                                 \
        _regpoly_gen_alias_##__LINE__ =                               \
            GeneratorRegistry::reg_alias((Alias), (CanonicalName));   \
    }

}  // namespace regpoly::core
