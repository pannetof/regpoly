// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "generator.h"
#include "param_spec.h"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

// Forward declaration; full definition pulled in only by the .cpp.
// This keeps NTL out of every translation unit that includes this header.
namespace regpoly::core {

class GenEnumerator;

class MarsaXorshiftGen : public Generator {
public:
    // Type 1: single word xorshift
    struct Type1Params {
        int a, b, c;
    };

    // Type 2x: two-component
    struct Type2xParams {
        std::vector<int> p;  // 3 values
        std::vector<int> q;  // 3 values
    };

    // Type 3: multi-tap
    struct Tap {
        int position;
        int shift;
    };

    // Type 4: two-component, 2 shifts each
    struct Type4Params {
        std::vector<int> p;  // 2 values
        std::vector<int> q;  // 2 values
    };

    // Type 100 (general): multi-component
    struct MiEntry {
        int position;
        std::vector<int> shifts;
    };

    MarsaXorshiftGen(int type, int w, int r, int m,
                     const Type1Params& t1,
                     const Type2xParams& t2x,
                     const std::vector<Tap>& taps,
                     const Type4Params& t4,
                     const std::vector<MiEntry>& mi,
                     int L);

    static std::unique_ptr<Generator> from_params(const Params& params, int L);
    static std::vector<ParamSpec> param_specs();

    // Build the exhaustive-search enumerator.  Requires `type` and `w`
    // in `resolved`; types 2, 3, 4, 100 additionally require `r`; types
    // 2 and 4 also accept an optional `m` pin; type 100 accepts
    // optional `nb_taps` (default 3) and `shifts_per_tap`
    // (per-tap vector, default `[1, 1, ..., 1]` of length nb_taps).
    // Throws std::invalid_argument with a "needs_<axis>" message
    // when a required input is missing or inconsistent.
    //
    // Axes per type (outermost first; all coordinates are encoded as
    // mixed-radix digits over NTL::ZZ so the enumerator scales to
    // very large totals):
    //
    //   type=1:    pattern (4) × a (w-1) × b (w-1) × c (w-1)
    //   type=2:    m (r-1) × p[0..2] (2w-1) × q[0..2] (2w-1)
    //   type=3:    tap-subset (C(r,3)) × sh[0..2] (2(w-1))
    //   type=4:    m (r-1) × p[0..1] (2w-1) × q[0..1] (2w-1)
    //   type=100:  tap-subset (C(r, nb_taps))
    //               × per-tap shift slots (2(w-1) each, total
    //                 Σ shifts_per_tap[i] slots)
    static std::unique_ptr<GenEnumerator> make_enumerator(
        const Params& resolved, int L);

    std::string name() const override;
    std::string display_str() const override;
    void init(const BitVect& init_bv) override;
    void next() override;
    std::unique_ptr<Generator> copy() const override;

private:
    // All structural and configuration fields are set in the ctor and
    // never mutated afterwards.  Marking them `const` makes that
    // invariant compiler-enforced.
    const int type_;
    const int w_;
    const int r_;
    const int m_;
    const Type1Params t1_;
    const Type2xParams t2x_;
    const std::vector<Tap> taps_;
    const Type4Params t4_;
    const std::vector<MiEntry> mi_;
    const uint64_t wmask_;            // pre-computed (1 << w_) - 1; 0xFFF..FFFF when w_ >= 64

    static uint64_t ShiftR(uint64_t v, int s);
    uint64_t V(int idx) const;
    void SetV(int idx, uint64_t val);

    // BitVect-backed slow path used when w > 64 (the uint64_t kernel
    // cannot represent words wider than 64 bits).  Same recurrence
    // semantics as the fast path.  Defined in marsaxorshift.cpp.
    void next_wide_();
};

}  // namespace regpoly::core
