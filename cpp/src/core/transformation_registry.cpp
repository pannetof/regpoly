// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "transformation_registry.h"

#include <stdexcept>

using namespace regpoly::core;


namespace regpoly::core {

namespace {

// Function-local statics to avoid static-initialisation order issues.
// Storage is a deque so pointer stability is guaranteed across reg().
std::unordered_map<std::string, TransformationRegistry::Info*>& by_name() {
    static std::unordered_map<std::string, TransformationRegistry::Info*> m;
    return m;
}

std::deque<TransformationRegistry::Info>& storage() {
    static std::deque<TransformationRegistry::Info> v;
    return v;
}

}  // namespace

int TransformationRegistry::reg(const std::string& canonical_name,
                                FromParamsFn from_params,
                                ParamSpecsFn param_specs) {
    auto& m = by_name();
    if (m.find(canonical_name) != m.end()) return 0;
    auto& slot = storage().emplace_back(Info{
        canonical_name,
        std::move(from_params),
        std::move(param_specs),
    });
    m.emplace(canonical_name, &slot);
    return 0;
}

const TransformationRegistry::Info& TransformationRegistry::lookup(const std::string& name) {
    auto* p = find(name);
    if (!p) throw std::invalid_argument("Unknown transformation type: " + name);
    return *p;
}

const TransformationRegistry::Info* TransformationRegistry::find(const std::string& name) {
    auto& m = by_name();
    auto it = m.find(name);
    return it == m.end() ? nullptr : it->second;
}

}  // namespace regpoly::core
