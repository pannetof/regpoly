#pragma once

#include "combined.h"
#include "generator.h"
#include "transformation.h"

#include <cstdint>
#include <memory>
#include <vector>

// Phase 2.4b-pre: C++ port of regpoly.core.{component,combination}.
//
// Combination is a stateful iterator over the cartesian product of J
// component generator pools, with two constraints (matching the
// existing Python semantics):
//
//   1. Identity-uniqueness: the same Generator object never appears
//      twice in a single combo (this would XOR the component with
//      itself and produce the zero sequence).
//   2. Shared-pool C(n,k) selection: when two component slots reference
//      the same pool object (set up via Component::share_pool_with),
//      the indices in those slots must be strictly increasing — so a
//      pool of size n shared across k slots yields C(n,k) combos
//      rather than the n^k permutations.
//
// On each successful reset()/next(), Combination updates k_g (sum of
// active generator k's) and L (min active L, capped at Lmax).

class Component {
public:
    // Owned pool of generators (each unique_ptr is the sole owner).
    using GenPool = std::vector<std::unique_ptr<Generator>>;
    using TransChain = std::vector<std::unique_ptr<Transformation>>;

    Component();
    Component(const Component&) = delete;
    Component& operator=(const Component&) = delete;
    Component(Component&&) = default;
    Component& operator=(Component&&) = default;

    // Append a deep copy of `gen` to this component's own pool. Throws
    // if this component currently shares its pool with another.
    void add_gen(const Generator& gen);

    // Append a deep copy of `t` to this component's tempering chain.
    void add_trans(const Transformation& t);

    // Make this component reference `other`'s pool (shared_ptr semantics).
    // After this call, add_gen on this component is forbidden.
    void share_pool_with(const Component& other);

    int nb_gen() const;
    int nb_trans() const;
    int current_gen() const { return current_gen_; }
    void set_current_gen(int i) { current_gen_ = i; }

    // Access the generator at index i (does not advance state).
    Generator& gen_at(int i) const;

    // Access the active generator (= gen_at(current_gen_)).
    Generator& active_gen() const;

    // Access the transformation chain.
    const TransChain& trans() const { return trans_; }
    Transformation& trans_at(int i) const;

    // Identity used by Combination's shared-pool detection. Two
    // components share a pool iff their pool_id() pointers compare
    // equal.
    const GenPool* pool_id() const { return pool_.get(); }

    // Diagnostic: concatenate `display_str()` of every transformation,
    // separated by newlines. Mirrors Python's Component.display().
    std::string display() const;

private:
    std::shared_ptr<GenPool> pool_;     // never null after construction
    TransChain trans_;
    int current_gen_;
    bool owns_pool_;
};


class Combination {
public:
    Combination(int J, int Lmax);
    Combination(const Combination&) = delete;
    Combination& operator=(const Combination&) = delete;
    Combination(Combination&&) = default;
    Combination& operator=(Combination&&) = default;

    int J() const { return J_; }
    int Lmax() const { return Lmax_; }
    int k_g() const { return k_g_; }
    int L() const { return L_; }

    Component& component(int j);
    const Component& component(int j) const;

    // Active generator at slot j (after reset()/next()). Equivalent to
    // Python's `comb[j]`.
    Generator& at(int j) const;

    // Reset to the first valid combination. Returns false if no valid
    // combo exists (e.g. some component has no generators).
    bool reset();

    // Advance to the next valid combination. Returns false on exhaustion.
    // After exhaustion, subsequent next() calls keep returning false.
    bool next();

    // True iff the iterator has been exhausted (next() returned false).
    bool exhausted() const { return exhausted_; }

private:
    int J_;
    int Lmax_;
    int k_g_;
    int L_;
    std::vector<std::shared_ptr<Component>> components_;
    std::vector<int> indices_;   // indices_[j] = current index in
                                 // components_[j]'s pool; -1 = unset
    bool exhausted_;

    // Recompute k_g and L from the current active generators.
    void update_stats();

    // Place valid indices at slots j..J-1 starting from a fresh state.
    // On success, indices_[j..J-1] are populated and components_'
    // current_gen are set.
    bool place_from(int j);

    // Try to advance the index at slot j (and reset slots j+1..J-1).
    // On exhaustion at slot j, recursively try slot j-1.
    bool advance_from(int j);

    // Compute the minimum allowed start index at slot j based on
    // shared-pool constraints with slots 0..j-1.
    int compute_min_index(int j) const;

    // Check whether a generator pointer at slot j collides with any
    // already-placed slot < j.
    bool already_used(const Generator* g, int j) const;
};


// Build a CombinedGenerator that snapshots the current active state
// of `comb`. Each active component generator is deep-copied via
// Generator::copy(), and each component's tempering chain is cloned
// via Transformation::copy(). The returned CombinedGenerator is
// independent of `comb` — subsequent next()/reset() on `comb` does
// not affect the returned object.
std::unique_ptr<CombinedGenerator>
build_combined_from_combination(const Combination& comb);
