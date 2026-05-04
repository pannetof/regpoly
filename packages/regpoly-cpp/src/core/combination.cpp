#include "combination.h"

#include <algorithm>
#include <climits>
#include <sstream>
#include <stdexcept>

// ── Component ──────────────────────────────────────────────────────────────

Component::Component()
    : pool_(std::make_shared<GenPool>()),
      current_gen_(0),
      owns_pool_(true) {}

void Component::add_gen(const Generator& gen) {
    if (!owns_pool_)
        throw std::logic_error(
            "Component::add_gen: cannot add to a shared pool — share_pool_with "
            "must be called only AFTER all add_gen calls on the source "
            "component.");
    pool_->push_back(gen.copy());
}

void Component::add_trans(const Transformation& t) {
    trans_.push_back(t.copy());
}

void Component::share_pool_with(const Component& other) {
    pool_ = other.pool_;       // shared_ptr aliasing — same vector object
    owns_pool_ = false;
    current_gen_ = 0;
}

int Component::nb_gen() const {
    return static_cast<int>(pool_->size());
}

int Component::nb_trans() const {
    return static_cast<int>(trans_.size());
}

Generator& Component::gen_at(int i) const {
    if (i < 0 || i >= static_cast<int>(pool_->size()))
        throw std::out_of_range("Component::gen_at: index out of range");
    return *(*pool_)[i];
}

Generator& Component::active_gen() const {
    return gen_at(current_gen_);
}

Transformation& Component::trans_at(int i) const {
    if (i < 0 || i >= static_cast<int>(trans_.size()))
        throw std::out_of_range("Component::trans_at: index out of range");
    return *trans_[i];
}

std::string Component::display() const {
    std::ostringstream oss;
    for (size_t i = 0; i < trans_.size(); ++i) {
        if (i > 0) oss << "\n";
        oss << trans_[i]->display_str();
    }
    return oss.str();
}

// ── Combination ────────────────────────────────────────────────────────────

Combination::Combination(int J, int Lmax)
    : J_(J), Lmax_(Lmax), k_g_(0), L_(0),
      indices_(J, -1),
      exhausted_(false)
{
    components_.reserve(J);
    for (int j = 0; j < J; ++j)
        components_.push_back(std::make_shared<Component>());
}

Component& Combination::component(int j) {
    if (j < 0 || j >= J_)
        throw std::out_of_range("Combination::component: index out of range");
    return *components_[j];
}

const Component& Combination::component(int j) const {
    if (j < 0 || j >= J_)
        throw std::out_of_range("Combination::component: index out of range");
    return *components_[j];
}

Generator& Combination::at(int j) const {
    return component(j).active_gen();
}

void Combination::update_stats() {
    k_g_ = 0;
    int min_L = INT_MAX;
    for (int j = 0; j < J_; ++j) {
        Generator& g = components_[j]->active_gen();
        k_g_ += g.k();
        if (g.L() < min_L) min_L = g.L();
    }
    L_ = (min_L > Lmax_) ? Lmax_ : min_L;
}

int Combination::compute_min_index(int j) const {
    int min_idx = 0;
    const auto* my_pool = components_[j]->pool_id();
    for (int p = 0; p < j; ++p) {
        if (components_[p]->pool_id() == my_pool)
            min_idx = std::max(min_idx, indices_[p] + 1);
    }
    return min_idx;
}

bool Combination::already_used(const Generator* g, int j) const {
    for (int p = 0; p < j; ++p) {
        const Generator* other = &components_[p]->gen_at(indices_[p]);
        if (other == g) return true;
    }
    return false;
}

bool Combination::place_from(int j) {
    if (j >= J_) return true;
    const int n = components_[j]->nb_gen();
    if (n == 0) return false;

    const int start = compute_min_index(j);
    for (int i = start; i < n; ++i) {
        Generator* g = &components_[j]->gen_at(i);
        if (!already_used(g, j)) {
            indices_[j] = i;
            components_[j]->set_current_gen(i);
            if (place_from(j + 1)) return true;
        }
    }
    indices_[j] = -1;
    return false;
}

bool Combination::advance_from(int j) {
    if (j < 0) return false;

    const int n = components_[j]->nb_gen();
    const int start = compute_min_index(j);
    const int try_from = std::max(indices_[j] + 1, start);

    for (int i = try_from; i < n; ++i) {
        Generator* g = &components_[j]->gen_at(i);
        if (!already_used(g, j)) {
            indices_[j] = i;
            components_[j]->set_current_gen(i);
            if (place_from(j + 1)) return true;
        }
    }
    // Exhausted at slot j — carry left.
    indices_[j] = -1;
    return advance_from(j - 1);
}

bool Combination::reset() {
    if (J_ == 0) {
        exhausted_ = true;
        return false;
    }
    for (int j = 0; j < J_; ++j) {
        if (components_[j]->nb_gen() == 0) {
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

bool Combination::next() {
    if (exhausted_) return false;
    if (!advance_from(J_ - 1)) {
        exhausted_ = true;
        return false;
    }
    update_stats();
    return true;
}
