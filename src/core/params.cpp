// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "params.h"

using namespace regpoly::core;


namespace regpoly::core {

void ParamBag::set_int(const std::string& key, int64_t val) { ints_[key] = val; }
void ParamBag::set_bool(const std::string& key, bool val) { bools_[key] = val; }
void ParamBag::set_string(const std::string& key, const std::string& val) { strings_[key] = val; }
void ParamBag::set_int_vec(const std::string& key, const std::vector<int>& val) { int_vecs_[key] = val; }
void ParamBag::set_uint_vec(const std::string& key, const std::vector<uint64_t>& val) { uint_vecs_[key] = val; }
void ParamBag::set_struct_map(const std::string& key, StructMap val) { struct_maps_[key] = std::move(val); }

int64_t ParamBag::get_int(const std::string& key, int64_t def) const {
    auto it = ints_.find(key);
    return it != ints_.end() ? it->second : def;
}

bool ParamBag::get_bool(const std::string& key, bool def) const {
    auto it = bools_.find(key);
    return it != bools_.end() ? it->second : def;
}

std::string ParamBag::get_string(const std::string& key, const std::string& def) const {
    auto it = strings_.find(key);
    return it != strings_.end() ? it->second : def;
}

std::vector<int> ParamBag::get_int_vec(const std::string& key) const {
    auto it = int_vecs_.find(key);
    if (it != int_vecs_.end()) return it->second;
    // Fall back to uint_vecs_ — the Python binding stores list params as
    // int_vec or uint_vec based on a per-value heuristic that doesn't know
    // the target type; allow callers to read either side regardless.
    auto uit = uint_vecs_.find(key);
    if (uit != uint_vecs_.end()) {
        std::vector<int> out;
        out.reserve(uit->second.size());
        for (uint64_t v : uit->second) out.push_back(static_cast<int>(v));
        return out;
    }
    return std::vector<int>{};
}

std::vector<uint64_t> ParamBag::get_uint_vec(const std::string& key) const {
    auto it = uint_vecs_.find(key);
    if (it != uint_vecs_.end()) return it->second;
    // Mirrored fallback: see get_int_vec().
    auto iit = int_vecs_.find(key);
    if (iit != int_vecs_.end()) {
        std::vector<uint64_t> out;
        out.reserve(iit->second.size());
        for (int v : iit->second) out.push_back(static_cast<uint64_t>(v));
        return out;
    }
    return std::vector<uint64_t>{};
}

const StructMap& ParamBag::get_struct_map(const std::string& key) const {
    static const StructMap empty;
    auto it = struct_maps_.find(key);
    return it != struct_maps_.end() ? it->second : empty;
}

bool ParamBag::has(const std::string& key) const {
    return ints_.count(key) || bools_.count(key) || strings_.count(key)
        || int_vecs_.count(key) || uint_vecs_.count(key)
        || struct_maps_.count(key);
}

bool ParamBag::has_struct_map(const std::string& key) const {
    return struct_maps_.count(key) > 0;
}

}  // namespace regpoly::core
