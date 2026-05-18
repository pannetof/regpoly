// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "generator.h"
#include "param_spec.h"
#include "params.h"

#include <cstdint>
#include <memory>
#include <string>

namespace regpoly::core {

class GenEnumerator;

// Blackman & Vigna (2022) xoroshiro linear engine.
//
// State: r words of w bits each, r >= 2.  Parameters (A, B, C) in [1, w-1].
// (The paper writes "k" for the word count; this codebase reserves k for
// the state size in bits, so we use r — matching the variable name the
// paper uses elsewhere, e.g. in Panneton & L'Ecuyer 2005.)
//
// The update follows the rw × rw matrix from §3.1:
//
//   new state = (old s_1, old s_2, ..., old s_{r-2},
//                rotl(s_0, A) ^ t ^ (t << B),
//                rotl(t, C))
//   where t = s_0 ^ s_{r-1}.
//
// The scrambler (+ / ++ / * / **) lives outside this class — this is the
// raw linear engine, suitable for primitivity / equidistribution analysis.
class XoroshiroGen : public Generator {
public:
    XoroshiroGen(int w, int r, int A, int B, int C, int L);

    static std::unique_ptr<Generator> from_params(const Params& params, int L);
    static std::vector<ParamSpec> param_specs();

    // Build the exhaustive-search enumerator.  Requires `w` and `r` in
    // `resolved` (the structural axes the user pins).  Throws
    // std::invalid_argument with a "needs_<axis>" message when an
    // input is missing or out of range.  Enumeration order (outer
    // first): A × B × C, each in [1, w-1].
    static std::unique_ptr<GenEnumerator> make_enumerator(
        const Params& resolved, int L);

    std::string name() const override;
    std::string display_str() const override;
    void init(const BitVect& init_bv) override;
    void next() override;
    std::unique_ptr<Generator> copy() const override;

    int w() const { return w_; }
    int r() const { return r_; }
    int A() const { return A_; }
    int B() const { return B_; }
    int C() const { return C_; }

private:
    const int w_;
    const int r_;
    const int A_;
    const int B_;
    const int C_;
    const uint64_t wmask_;

    uint64_t rotl_(uint64_t x, int rot) const;
};

}  // namespace regpoly::core
