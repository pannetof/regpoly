// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "generator_registry.h"

#include <deque>
#include <stdexcept>

namespace {

// Lazily-initialised singletons. Function-local statics avoid the
// static-initialisation-order fiasco: REGISTER_GENERATOR macros run at
// dynamic-init time and call reg(), which calls these accessors — the
// tables are constructed on first use.
//
// Storage is a deque, NOT a vector: deque guarantees pointer stability
// across emplace_back. by_name() and canonical_order() store raw
// pointers into storage() and would dangle if a vector reallocated.
std::unordered_map<std::string, GeneratorRegistry::Info*>& by_name() {
    static std::unordered_map<std::string, GeneratorRegistry::Info*> m;
    return m;
}

std::deque<GeneratorRegistry::Info>& storage() {
    static std::deque<GeneratorRegistry::Info> v;
    return v;
}

std::vector<GeneratorRegistry::Info*>& canonical_order() {
    static std::vector<GeneratorRegistry::Info*> v;
    return v;
}

std::vector<std::pair<std::string, std::string>>& alias_pairs() {
    static std::vector<std::pair<std::string, std::string>> v;
    return v;
}

}  // namespace

int GeneratorRegistry::reg(const std::string& canonical_name,
                           FromParamsFn from_params,
                           ParamSpecsFn param_specs,
                           BindFn bind,
                           EnumeratorFn make_enumerator) {
    auto& m = by_name();
    auto it = m.find(canonical_name);
    if (it != m.end()) {
        // Already registered. Defensively treat as no-op rather than
        // throwing — repeat registration can happen if the same
        // translation unit is linked from multiple targets.
        return 0;
    }
    auto& slot = storage().emplace_back(Info{
        canonical_name,
        std::move(from_params),
        std::move(param_specs),
        std::move(bind),
        std::move(make_enumerator),
    });
    m.emplace(canonical_name, &slot);
    canonical_order().push_back(&slot);
    return 0;
}

int GeneratorRegistry::reg_alias(const std::string& alias,
                                 const std::string& canonical_name) {
    auto& m = by_name();
    auto it = m.find(canonical_name);
    if (it == m.end()) {
        // Canonical not yet registered; defer would be cleaner but
        // every current call site declares the canonical first. Throw
        // to surface the misuse loudly.
        throw std::logic_error(
            "GeneratorRegistry::reg_alias: canonical name not registered: "
            + canonical_name);
    }
    m.emplace(alias, it->second);
    alias_pairs().emplace_back(alias, canonical_name);
    return 0;
}

const GeneratorRegistry::Info& GeneratorRegistry::lookup(const std::string& name) {
    auto* p = find(name);
    if (!p) throw std::invalid_argument("Unknown generator family: " + name);
    return *p;
}

const GeneratorRegistry::Info* GeneratorRegistry::find(const std::string& name) {
    auto& m = by_name();
    auto it = m.find(name);
    return it == m.end() ? nullptr : it->second;
}

const std::vector<GeneratorRegistry::Info*>& GeneratorRegistry::canonical_entries() {
    return canonical_order();
}

const std::vector<std::pair<std::string, std::string>>& GeneratorRegistry::aliases() {
    return alias_pairs();
}
