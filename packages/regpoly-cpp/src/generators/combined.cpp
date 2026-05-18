// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "combined.h"

#include <NTL/GF2X.h>
#include <algorithm>
#include <sstream>
#include <stdexcept>

using namespace regpoly::core;


namespace regpoly::core {

namespace {

NTL::GF2X char_poly_to_gf2x(const Generator& gen) {
    BitVect bv = gen.char_poly();
    int K = gen.k();
    NTL::GF2X p;
    NTL::SetCoeff(p, K);  // implicit leading z^K term
    for (int j = 0; j < K; ++j)
        if (bv.get_bit(j))
            NTL::SetCoeff(p, j);
    return p;
}

}  // namespace

int CombinedGenerator::compute_k(
    const std::vector<std::unique_ptr<Generator>>& comps)
{
    int total = 0;
    for (const auto& c : comps)
        total += c->k();
    return total;
}

int CombinedGenerator::compute_L(
    const std::vector<std::unique_ptr<Generator>>& comps,
    int Lmax)
{
    if (comps.empty())
        return Lmax;
    int L = comps[0]->L();
    for (size_t i = 1; i < comps.size(); ++i)
        L = std::min(L, comps[i]->L());
    return std::min(L, Lmax);
}

CombinedGenerator::CombinedGenerator(
    std::vector<std::unique_ptr<Generator>> components,
    int Lmax)
    : Generator(compute_k(components), compute_L(components, Lmax)),
      components_(std::move(components)),
      tempering_(components_.size())
{
    prefix_k_.reserve(components_.size() + 1);
    prefix_k_.push_back(0);
    for (const auto& c : components_)
        prefix_k_.push_back(prefix_k_.back() + c->k());
    refresh_concatenated_state();
}

CombinedGenerator::CombinedGenerator(
    std::vector<std::unique_ptr<Generator>> components,
    std::vector<ComponentTempering> tempering_chains,
    int Lmax)
    : Generator(compute_k(components), compute_L(components, Lmax)),
      components_(std::move(components)),
      tempering_(std::move(tempering_chains))
{
    if (tempering_.size() != components_.size()) {
        throw std::invalid_argument(
            "CombinedGenerator: tempering chain count must equal "
            "component count");
    }
    prefix_k_.reserve(components_.size() + 1);
    prefix_k_.push_back(0);
    for (const auto& c : components_)
        prefix_k_.push_back(prefix_k_.back() + c->k());
    refresh_concatenated_state();
}

std::string CombinedGenerator::name() const {
    std::ostringstream os;
    os << "CombinedGenerator[J=" << components_.size() << "]";
    return os.str();
}

std::string CombinedGenerator::display_str() const {
    std::ostringstream os;
    os << "CombinedGenerator(J=" << components_.size()
       << ", k=" << k() << ", L=" << L() << ")";
    for (size_t j = 0; j < components_.size(); ++j) {
        os << "\n  [" << j << "] " << components_[j]->display_str();
    }
    return os.str();
}

void CombinedGenerator::init(const BitVect& init_bv) {
    int total = k();
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
    refresh_concatenated_state();
}

void CombinedGenerator::next() {
    for (auto& c : components_)
        c->next();
    refresh_concatenated_state();
}

BitVect CombinedGenerator::get_output() const {
    int Lout = L();
    BitVect out(Lout);
    for (size_t j = 0; j < components_.size(); ++j) {
        BitVect raw = components_[j]->get_output();
        for (auto& t : tempering_[j])
            t->apply(raw);
        int n = std::min(Lout, raw.nbits());
        for (int i = 0; i < n; ++i)
            if (raw.get_bit(i))
                out.set_bit(i, out.get_bit(i) ^ 1);
    }
    return out;
}

std::unique_ptr<Generator> CombinedGenerator::copy() const {
    std::vector<std::unique_ptr<Generator>> comps;
    comps.reserve(components_.size());
    for (const auto& c : components_)
        comps.push_back(c->copy());

    std::vector<ComponentTempering> chains;
    chains.reserve(tempering_.size());
    for (const auto& chain : tempering_) {
        ComponentTempering copy_chain;
        copy_chain.reserve(chain.size());
        for (const auto& t : chain)
            copy_chain.push_back(t->copy());
        chains.push_back(std::move(copy_chain));
    }

    auto cg = std::make_unique<CombinedGenerator>(
        std::move(comps), std::move(chains), L());
    // Prefix is reconstructed by the constructor; state was just
    // refreshed there from the freshly cloned component states.
    return cg;
}

BitVect CombinedGenerator::char_poly() const {
    NTL::GF2X product;
    NTL::SetCoeff(product, 0);  // start with constant polynomial 1
    for (const auto& c : components_)
        product *= char_poly_to_gf2x(*c);

    int K = k();
    BitVect out(K);
    int dp = NTL::deg(product);
    int upper = std::min(K, dp);
    for (int j = 0; j < upper; ++j)
        if (NTL::IsOne(NTL::coeff(product, j)))
            out.set_bit(j, 1);
    // The leading z^K term is implicit per the Generator convention; no
    // bit is stored for it.
    return out;
}

std::vector<Generator*> CombinedGenerator::raw_component_pointers() const {
    std::vector<Generator*> out;
    out.reserve(components_.size());
    for (const auto& c : components_)
        out.push_back(c.get());
    return out;
}

std::vector<std::vector<Transformation*>>
CombinedGenerator::raw_tempering_pointers() const {
    std::vector<std::vector<Transformation*>> out;
    out.reserve(tempering_.size());
    for (const auto& chain : tempering_) {
        std::vector<Transformation*> raw;
        raw.reserve(chain.size());
        for (const auto& t : chain)
            raw.push_back(t.get());
        out.push_back(std::move(raw));
    }
    return out;
}

void CombinedGenerator::refresh_concatenated_state() {
    BitVect concat(k());
    for (size_t j = 0; j < components_.size(); ++j) {
        const BitVect& s = components_[j]->state();
        int k_j = components_[j]->k();
        int off = prefix_k_[j];
        for (int i = 0; i < k_j; ++i)
            if (s.get_bit(i))
                concat.set_bit(off + i, 1);
    }
    set_raw_state(concat);
}

std::optional<std::string> CombinedGenerator::compute_default_test_method(const std::string& test_type) const {
    // J=1 is a degenerate combined that wraps a single primitive: defer
    // to the wrapped component so the wrapper is observationally equal
    // to the primitive (matches the JEquals1MatchesPrimitive invariant).
    // Goes through the component's public default_test_method so the
    // component's own cache is consulted/populated.
    if (components_.size() == 1)
        return components_[0]->default_test_method(test_type);
    if (test_type == "equidistribution") {
        // If every component agrees on a non-empty preferred method,
        // honor it.  This lets families like CellularAutomataGen opt
        // into the matricial method even when wrapped in a combined
        // generator (paper convention treats the combined matrix as a
        // single rank target).  Falls back to "notprimitive" — the safe
        // default for reducible combined characteristic polynomials.
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
