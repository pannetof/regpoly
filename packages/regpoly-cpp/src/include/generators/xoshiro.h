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

// Blackman & Vigna (2022) xoshiro linear engine.
//
// State: r words of w bits each, r in {4, 8}.  Parameters (A, B) in
// [1, w-1].  (The paper writes "k" for the word count; this codebase
// reserves k for the state size in bits, so we use r.)  The update
// follows the 4w × 4w (resp. 8w × 8w) matrix from §3.2:
//
//   r = 4 (S_4w):
//     t = s[1] << A;  s[2] ^= s[0];  s[3] ^= s[1];
//     s[1] ^= s[2];   s[0] ^= s[3];  s[2] ^= t;
//     s[3] = rotl(s[3], B);
//
//   r = 8 (S_8w):
//     t = s[1] << A;
//     s[2] ^= s[0];   s[5] ^= s[1];  s[1] ^= s[2];
//     s[7] ^= s[3];   s[3] ^= s[4];  s[4] ^= s[5];
//     s[0] ^= s[6];   s[6] ^= s[7];  s[6] ^= t;
//     s[7] = rotl(s[7], B);
//
// Scramblers (+ / ++ / **) are NOT implemented — this is the raw linear
// engine.
class XoshiroGen : public Generator {
public:
    XoshiroGen(int w, int r, int A, int B, int L);

    static std::unique_ptr<Generator> from_params(const Params& params, int L);
    static std::vector<ParamSpec> param_specs();

    // Build the exhaustive-search enumerator.  Requires `w` and `r` in
    // `resolved` (r must be 4 or 8).  Throws std::invalid_argument
    // with a "needs_<axis>" message on bad inputs.  Enumeration order
    // (outer first): A × B, each in [1, w-1].
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

private:
    const int w_;
    const int r_;           // 4 or 8
    const int A_;
    const int B_;
    const uint64_t wmask_;

    uint64_t rotl_(uint64_t x, int rot) const;

    void next_4_();
    void next_8_();
};

}  // namespace regpoly::core
