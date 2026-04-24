#pragma once
#include <string>
#include <unordered_map>
#include <vector>
#include <cstdint>

class Params {
public:
    void set_int(const std::string& key, int64_t val);
    void set_bool(const std::string& key, bool val);
    void set_string(const std::string& key, const std::string& val);
    void set_int_vec(const std::string& key, const std::vector<int>& val);
    void set_uint_vec(const std::string& key, const std::vector<uint64_t>& val);

    int64_t get_int(const std::string& key, int64_t def = 0) const;
    bool get_bool(const std::string& key, bool def = false) const;
    std::string get_string(const std::string& key, const std::string& def = "") const;
    std::vector<int> get_int_vec(const std::string& key) const;
    std::vector<uint64_t> get_uint_vec(const std::string& key) const;
    bool has(const std::string& key) const;

    // Read-only views for serialisers (e.g. params_to_dict in bindings).
    const std::unordered_map<std::string, int64_t>& ints() const { return ints_; }
    const std::unordered_map<std::string, bool>& bools() const { return bools_; }
    const std::unordered_map<std::string, std::string>& strings() const { return strings_; }
    const std::unordered_map<std::string, std::vector<int>>& int_vecs() const { return int_vecs_; }
    const std::unordered_map<std::string, std::vector<uint64_t>>& uint_vecs() const { return uint_vecs_; }

private:
    std::unordered_map<std::string, int64_t> ints_;
    std::unordered_map<std::string, bool> bools_;
    std::unordered_map<std::string, std::string> strings_;
    std::unordered_map<std::string, std::vector<int>> int_vecs_;
    std::unordered_map<std::string, std::vector<uint64_t>> uint_vecs_;
};
