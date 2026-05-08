// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

#pragma once
#include "generator.h"
#include "param_spec.h"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <functional>

class WELLGen : public Generator {
public:
    // One algorithm slot's transformation. `Mi` selects which paper M-class
    // (Panneton et al. 2006, Table I); the active fields depend on Mi:
    //
    //   M0 = 0       y = 0                            (no args)
    //   M1 = 1       y = x                            (no args)
    //   M2(t) = 2    y = x ≫ t  /  x ≪ −t            (t signed)
    //   M3(t) = 3    y = x ⊕ shift(x, t)              (t signed)
    //   M4(a) = 4    y = (x ≫ 1) ⊕ a if LSB(x)        (a 32-bit)
    //   M5(t,b) = 5  y = x ⊕ (shift(x, t) & b)         (t signed, b 32-bit)
    //   M6(q,t,s,a) = 6  rotate-mask-conditional XOR  (q,t,s ∈ [0,w-1], a 32-bit)
    //
    // M6's d_s mask and (1 << t) test bit are synthesised at runtime from
    // the small-integer fields (matches paper Table I exactly).
    struct MatrixEntry {
        int Mi = 0;        // 0..6
        int q  = 0;        // M6: rotate amount
        int t  = 0;        // M2/M3/M5: signed shift; M6: bit position 0..w-1
        int s  = 0;        // M6: mask selector 0..w-1
        uint32_t a = 0;    // M4, M6: XOR constant
        uint32_t b = 0;    // M5: mask
    };

    WELLGen(int w, int r, int p, int m1, int m2, int m3,
              const std::vector<MatrixEntry>& matrices, int L);

    static std::unique_ptr<Generator> from_params(const Params& params, int L);
    static std::vector<ParamSpec> param_specs();

    std::string name() const override;
    std::string display_str() const override;
    void init(const BitVect& init_bv) override;
    void next() override;
    std::unique_ptr<Generator> copy() const override;
    BitVect get_output() const override;
    void get_transition_state(uint64_t* out_words, int out_nwords) const override;

private:
    int w_;
    int r_;
    int p_;         // output mask bits
    int m1_, m2_, m3_;
    // 8 algorithm slots T0..T7 in the WELL recurrence; each slot's
    // class is one of M0..M6 from Table I of Panneton et al. (2006).
    std::vector<MatrixEntry> matrices_;
    int i_;         // circular buffer pointer
    int state_bits_; // w * r (full state size for circular buffer)
    uint64_t maskp_;   // UPPER mask: p least significant bits
    uint64_t umaskp_;  // LOWER mask: ~maskp
    static constexpr uint64_t M32 = 0xFFFFFFFFULL;

    static uint64_t ShiftR(uint64_t v, int s);
    static uint64_t apply_matrix(const MatrixEntry& m, uint64_t v);
    uint64_t TMAT(int j, uint64_t val) const;
    uint64_t V(int idx) const;
    void SetV(int idx, uint64_t val);

    static int type_cost(int Mi);
    static std::string type_display(const MatrixEntry& m);
};
