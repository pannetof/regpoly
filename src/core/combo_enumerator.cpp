// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#include "combo_enumerator.h"
#include "combined_f2_linear_source.h"

#include <algorithm>
#include <climits>
#include <sstream>
#include <stdexcept>

using namespace regpoly::core;


// ── F2LinearSourcePool ────────────────────────────────────────────────────

namespace regpoly::core {

F2LinearSourcePool::F2LinearSourcePool()
    : pool_(std::make_shared<SourcePool>()),
      current_index_(0),
      owns_pool_(true) {}

void F2LinearSourcePool::add_source(const F2LinearSource& src) {
    if (!owns_pool_)
        throw std::logic_error(
            "F2LinearSourcePool::add_source: cannot add to a shared pool — "
            "copy_pool_from must be called only AFTER all add_source calls on "
            "the source pool.");
    pool_->push_back(src.clone_source());
}

void F2LinearSourcePool::add_trans(const Transformation& t) {
    trans_.push_back(t.copy());
}

void F2LinearSourcePool::copy_pool_from(const F2LinearSourcePool& other) {
    pool_ = other.pool_;       // shared_ptr aliasing — same vector object
    owns_pool_ = false;
    current_index_ = 0;
}

int F2LinearSourcePool::nb_sources() const {
    return static_cast<int>(pool_->size());
}

int F2LinearSourcePool::nb_trans() const {
    return static_cast<int>(trans_.size());
}

F2LinearSource& F2LinearSourcePool::source_at(int i) const {
    if (i < 0 || i >= static_cast<int>(pool_->size()))
        throw std::out_of_range("F2LinearSourcePool::source_at: index out of range");
    return *(*pool_)[i];
}

F2LinearSource& F2LinearSourcePool::active_source() const {
    return source_at(current_index_);
}

Transformation& F2LinearSourcePool::trans_at(int i) const {
    if (i < 0 || i >= static_cast<int>(trans_.size()))
        throw std::out_of_range("F2LinearSourcePool::trans_at: index out of range");
    return *trans_[i];
}

std::string F2LinearSourcePool::display() const {
    std::ostringstream oss;
    for (size_t i = 0; i < trans_.size(); ++i) {
        if (i > 0) oss << "\n";
        oss << trans_[i]->display_str();
    }
    return oss.str();
}

// ── ComboEnumerator ────────────────────────────────────────────────────────────

ComboEnumerator::ComboEnumerator(int J, int Lmax)
    : J_(J), Lmax_(Lmax), k_g_(0), L_(0),
      indices_(J, -1),
      exhausted_(false)
{
    pools_.reserve(J);
    for (int j = 0; j < J; ++j)
        pools_.push_back(std::make_shared<F2LinearSourcePool>());
}

F2LinearSourcePool& ComboEnumerator::pool(int j) {
    if (j < 0 || j >= J_)
        throw std::out_of_range("ComboEnumerator::pool: index out of range");
    return *pools_[j];
}

const F2LinearSourcePool& ComboEnumerator::pool(int j) const {
    if (j < 0 || j >= J_)
        throw std::out_of_range("ComboEnumerator::pool: index out of range");
    return *pools_[j];
}

F2LinearSource& ComboEnumerator::at(int j) const {
    return pool(j).active_source();
}

void ComboEnumerator::update_stats() {
    k_g_ = 0;
    int min_L = INT_MAX;
    for (int j = 0; j < J_; ++j) {
        F2LinearSource& g = pools_[j]->active_source();
        k_g_ += g.k();
        if (g.L() < min_L) min_L = g.L();
    }
    L_ = (min_L > Lmax_) ? Lmax_ : min_L;
}

int ComboEnumerator::compute_min_index(int j) const {
    int min_idx = 0;
    const auto* my_pool = pools_[j]->pool_id();
    for (int p = 0; p < j; ++p) {
        if (pools_[p]->pool_id() == my_pool)
            min_idx = std::max(min_idx, indices_[p] + 1);
    }
    return min_idx;
}

bool ComboEnumerator::already_used(const F2LinearSource* g, int j) const {
    for (int p = 0; p < j; ++p) {
        const F2LinearSource* other = &pools_[p]->source_at(indices_[p]);
        if (other == g) return true;
    }
    return false;
}

bool ComboEnumerator::place_from(int j) {
    if (j >= J_) return true;
    const int n = pools_[j]->nb_sources();
    if (n == 0) return false;

    const int start = compute_min_index(j);
    for (int i = start; i < n; ++i) {
        F2LinearSource* g = &pools_[j]->source_at(i);
        if (!already_used(g, j)) {
            indices_[j] = i;
            pools_[j]->set_current_index(i);
            if (place_from(j + 1)) return true;
        }
    }
    indices_[j] = -1;
    return false;
}

bool ComboEnumerator::advance_from(int j) {
    if (j < 0) return false;

    const int n = pools_[j]->nb_sources();
    const int start = compute_min_index(j);
    const int try_from = std::max(indices_[j] + 1, start);

    for (int i = try_from; i < n; ++i) {
        F2LinearSource* g = &pools_[j]->source_at(i);
        if (!already_used(g, j)) {
            indices_[j] = i;
            pools_[j]->set_current_index(i);
            if (place_from(j + 1)) return true;
        }
    }
    // Exhausted at slot j — carry left.
    indices_[j] = -1;
    return advance_from(j - 1);
}

bool ComboEnumerator::reset() {
    if (J_ == 0) {
        exhausted_ = true;
        return false;
    }
    for (int j = 0; j < J_; ++j) {
        if (pools_[j]->nb_sources() == 0) {
            exhausted_ = true;
            return false;
        }
    }
    indices_.assign(J_, -1);
    exhausted_ = false;
    if (!place_from(0)) {
        exhausted_ = true;
        return false;
    }
    update_stats();
    return true;
}

bool ComboEnumerator::next() {
    if (exhausted_) return false;
    if (!advance_from(J_ - 1)) {
        exhausted_ = true;
        return false;
    }
    update_stats();
    return true;
}

// ── build_combined_from_enumerator ────────────────────────────────────────

std::unique_ptr<ITestable>
build_combined_from_enumerator(const ComboEnumerator& comb)
{
    // CombinedF2LinearSource holds F2LinearSource components — both
    // Recurrence-typed PRNGs and DigitalNet sources are allowed.
    // Kernels that require Recurrence components (matricial χ-recovery,
    // SIMD path) walk components() and dynamic_cast each; they throw
    // std::invalid_argument if any component is a DigitalNet.
    std::vector<std::unique_ptr<F2LinearSource>> comps;
    std::vector<TemperingChain> chains;
    comps.reserve(comb.J());
    chains.reserve(comb.J());

    for (int j = 0; j < comb.J(); ++j) {
        const F2LinearSourcePool& p = comb.pool(j);
        comps.push_back(p.active_source().clone_source());

        TemperingChain chain;
        for (int ti = 0; ti < p.nb_trans(); ++ti)
            chain.add(p.trans_at(ti).copy());
        chains.push_back(std::move(chain));
    }

    return std::make_unique<CombinedF2LinearSource>(
        std::move(comps), std::move(chains), comb.Lmax());
}

}  // namespace regpoly::core
