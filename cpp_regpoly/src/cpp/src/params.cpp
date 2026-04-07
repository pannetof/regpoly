#include "params.h"

void Params::set_int(const std::string& key, int64_t val) { ints_[key] = val; }
void Params::set_bool(const std::string& key, bool val) { bools_[key] = val; }
void Params::set_string(const std::string& key, const std::string& val) { strings_[key] = val; }
void Params::set_int_vec(const std::string& key, const std::vector<int>& val) { int_vecs_[key] = val; }
void Params::set_uint_vec(const std::string& key, const std::vector<uint64_t>& val) { uint_vecs_[key] = val; }

int64_t Params::get_int(const std::string& key, int64_t def) const {
    auto it = ints_.find(key);
    return it != ints_.end() ? it->second : def;
}

bool Params::get_bool(const std::string& key, bool def) const {
    auto it = bools_.find(key);
    return it != bools_.end() ? it->second : def;
}

std::string Params::get_string(const std::string& key, const std::string& def) const {
    auto it = strings_.find(key);
    return it != strings_.end() ? it->second : def;
}

std::vector<int> Params::get_int_vec(const std::string& key) const {
    auto it = int_vecs_.find(key);
    return it != int_vecs_.end() ? it->second : std::vector<int>{};
}

std::vector<uint64_t> Params::get_uint_vec(const std::string& key) const {
    auto it = uint_vecs_.find(key);
    return it != uint_vecs_.end() ? it->second : std::vector<uint64_t>{};
}

bool Params::has(const std::string& key) const {
    return ints_.count(key) || bools_.count(key) || strings_.count(key)
        || int_vecs_.count(key) || uint_vecs_.count(key);
}
