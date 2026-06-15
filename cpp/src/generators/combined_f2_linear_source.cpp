// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Francois Panneton, Ph.D.

#include "combined_f2_linear_source.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>

namespace regpoly::core {

namespace {

int sum_k(const std::vector<std::unique_ptr<F2LinearSource>>& comps) {
    int total = 0;
    for (const auto& c : comps) total += c->k();
    return total;
}

int min_L(const std::vector<std::unique_ptr<F2LinearSource>>& comps, int Lmax) {
    if (comps.empty()) return Lmax;
    int L = comps[0]->L();
    for (size_t i = 1; i < comps.size(); ++i)
        L = std::min(L, comps[i]->L());
    return std::min(L, Lmax);
}

}  // namespace

CombinedF2LinearSource::CombinedF2LinearSource(
    std::vector<std::unique_ptr<F2LinearSource>> components,
    int Lmax)
    : components_(std::move(components))
{
    if (components_.empty())
        throw std::invalid_argument(
            "CombinedF2LinearSource: at least one component is required");
    k_ = sum_k(components_);
    L_ = min_L(components_, Lmax);

    prefix_k_.reserve(components_.size() + 1);
    prefix_k_.push_back(0);
    for (const auto& c : components_)
        prefix_k_.push_back(prefix_k_.back() + c->k());
}

CombinedF2LinearSource::CombinedF2LinearSource(
    std::vector<std::unique_ptr<F2LinearSource>> components,
    std::vector<TemperingChain> tempering_chains,
    int Lmax)
    : components_(std::move(components))
{
    if (components_.empty())
        throw std::invalid_argument(
            "CombinedF2LinearSource: at least one component is required");
    if (tempering_chains.size() != components_.size()) {
        throw std::invalid_argument(
            "CombinedF2LinearSource: tempering chain count must equal "
            "component count");
    }
    for (size_t j = 0; j < components_.size(); ++j) {
        components_[j]->set_tempering(std::move(tempering_chains[j]));
    }
    k_ = sum_k(components_);
    L_ = min_L(components_, Lmax);
    prefix_k_.reserve(components_.size() + 1);
    prefix_k_.push_back(0);
    for (const auto& c : components_)
        prefix_k_.push_back(prefix_k_.back() + c->k());
}

std::string CombinedF2LinearSource::name() const {
    std::ostringstream os;
    os << "CombinedF2LinearSource[J=" << components_.size() << "]";
    return os.str();
}

std::string CombinedF2LinearSource::display_str() const {
    std::ostringstream os;
    os << "CombinedF2LinearSource(J=" << components_.size()
       << ", k=" << k_ << ", L=" << L_ << ")";
    for (size_t j = 0; j < components_.size(); ++j) {
        os << "\n  [" << j << "] " << components_[j]->display_str();
    }
    return os.str();
}

void CombinedF2LinearSource::init(const BitVect& init_bv) {
    int total = k_;
    BitVect padded(total);
    int copy_bits = std::min(total, init_bv.nbits());
    for (int i = 0; i < copy_bits; ++i)
        if (init_bv.get_bit(i))
            padded.set_bit(i, 1);

    for (size_t j = 0; j < components_.size(); ++j) {
        int k_j = components_[j]->k();
        BitVect slice(k_j);
        int off = prefix_k_[j];
        for (int i = 0; i < k_j; ++i)
            if (padded.get_bit(off + i))
                slice.set_bit(i, 1);
        components_[j]->init(slice);
    }
}

void CombinedF2LinearSource::next() {
    for (auto& c : components_) c->next();
}

BitVect CombinedF2LinearSource::get_output() const {
    BitVect out(L_);
    for (size_t j = 0; j < components_.size(); ++j) {
        BitVect tempered = components_[j]->get_output();
        int n = std::min(L_, tempered.nbits());
        for (int i = 0; i < n; ++i)
            if (tempered.get_bit(i))
                out.set_bit(i, out.get_bit(i) ^ 1);
    }
    return out;
}

std::unique_ptr<ITestable> CombinedF2LinearSource::copy() const {
    std::vector<std::unique_ptr<F2LinearSource>> clones;
    clones.reserve(components_.size());
    for (const auto& c : components_)
        clones.push_back(c->clone_source());
    return std::make_unique<CombinedF2LinearSource>(std::move(clones), L_);
}

std::vector<const F2LinearSource*> CombinedF2LinearSource::sources() const {
    std::vector<const F2LinearSource*> out;
    out.reserve(components_.size());
    for (const auto& c : components_)
        out.push_back(c.get());
    return out;
}

std::vector<F2LinearSource*> CombinedF2LinearSource::components() const {
    std::vector<F2LinearSource*> out;
    out.reserve(components_.size());
    for (const auto& c : components_)
        out.push_back(c.get());
    return out;
}

std::vector<Recurrence*>
CombinedF2LinearSource::recurrence_components(const char* caller) const {
    std::vector<Recurrence*> out;
    out.reserve(components_.size());
    for (size_t j = 0; j < components_.size(); ++j) {
        Recurrence* r = dynamic_cast<Recurrence*>(components_[j].get());
        if (!r) {
            throw std::invalid_argument(
                std::string(caller) +
                ": component " + std::to_string(j) + " (" +
                components_[j]->name() +
                ") is not a Recurrence; this kernel needs recurrence-"
                "state evolution and cannot operate on a DigitalNet.");
        }
        out.push_back(r);
    }
    return out;
}

std::optional<std::string>
CombinedF2LinearSource::default_test_method(const std::string& test_type) const {
    if (components_.size() == 1)
        return components_[0]->default_test_method(test_type);
    if (test_type == "equidistribution") {
        std::optional<std::string> unanimous;
        for (const auto& c : components_) {
            auto m = c->default_test_method(test_type);
            if (!m.has_value()) {
                unanimous.reset();
                break;
            }
            if (!unanimous.has_value())
                unanimous = m;
            else if (*unanimous != *m) {
                unanimous.reset();
                break;
            }
        }
        if (unanimous.has_value())
            return unanimous;
        return std::string("notprimitive");
    }
    return std::nullopt;
}

}  // namespace regpoly::core
