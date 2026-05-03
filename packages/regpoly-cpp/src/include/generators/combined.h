#pragma once

#include "generator.h"
#include "transformation.h"
#include <memory>
#include <utility>
#include <vector>

// CombinedGenerator — a Generator subclass that owns J component
// generators (each itself a Generator) plus a per-component tempering
// chain. `next()` advances all components in lockstep; the published
// L-bit output is the bitwise XOR of the components' (post-tempering)
// outputs.
//
// The combined characteristic polynomial is the product of the
// components' characteristic polynomials over GF(2). The combined
// transition matrix is block-diagonal across components.
//
// All lattice / equidistribution kernels in Phase 1+ consume a single
// `const Generator&` — primitive or combined, the API is uniform.
//
// Phase 1: tempering chains and the full SIMD-aware story (lane counts,
// output phases) are stored but `simd_lane_count()` / `output_phases()`
// fall through to the base defaults. Phase 2 widens this when the
// search-driver loop is migrated to C++ and SIMD-aware combined
// scenarios become a real concern.

class CombinedGenerator : public Generator {
public:
    using ComponentTempering =
        std::vector<std::unique_ptr<Transformation>>;

    // Components-only constructor (no tempering chains).
    CombinedGenerator(
        std::vector<std::unique_ptr<Generator>> components,
        int Lmax);

    // Components + per-component tempering chains. The chains are
    // applied to each component's output before XOR.
    CombinedGenerator(
        std::vector<std::unique_ptr<Generator>> components,
        std::vector<ComponentTempering> tempering_chains,
        int Lmax);

    std::string name() const override;
    std::string display_str() const override;

    // init(bv): bv must have at least k() bits. The leading k_0 bits go
    // to component 0, next k_1 to component 1, etc. If bv is too short
    // it is zero-padded; if too long the trailing bits are ignored.
    void init(const BitVect& init_bv) override;

    void next() override;

    // get_output(): XOR of components' (tempered) L-bit outputs,
    // truncated to L().
    BitVect get_output() const override;

    std::unique_ptr<Generator> copy() const override;

    // char_poly(): product over GF(2) of components' char polys.
    // Returns k() bits (no leading z^k term, matching Generator
    // convention).
    BitVect char_poly() const override;

    int J() const { return static_cast<int>(components_.size()); }
    const Generator& component(int j) const { return *components_[j]; }

    // Per-component k partition: prefix_k_[j] = sum_{i<j} k_i.
    // prefix_k_[J] = k() (total).
    const std::vector<int>& prefix_k() const { return prefix_k_; }

    // Internal accessors used by the lattice/equidistribution adapter
    // overloads. Phase 1 keeps the kernel implementations on their
    // original `vector<Generator*>` + per-component tempering API; the
    // new `const Generator&` overloads dispatch through these accessors
    // to recover the vector form. Phase 2+ migrates the implementations
    // to consume Generator& natively, at which point these can become
    // private again.
    std::vector<Generator*> raw_component_pointers() const;
    std::vector<std::vector<Transformation*>>
        raw_tempering_pointers() const;

private:
    std::vector<std::unique_ptr<Generator>> components_;
    std::vector<ComponentTempering> tempering_;
    std::vector<int> prefix_k_;

    int compute_k(
        const std::vector<std::unique_ptr<Generator>>& comps);
    int compute_L(
        const std::vector<std::unique_ptr<Generator>>& comps,
        int Lmax);

    void refresh_concatenated_state();
};
