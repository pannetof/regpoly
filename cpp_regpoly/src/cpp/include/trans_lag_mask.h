#pragma once
#include "transformation.h"
#include "param_spec.h"
#include <cstdint>

/**
 * Two-reference output transformation (MELG-style).
 *
 * Applies:
 *   y  = output_word ^ (output_word << sigma)     (T1: left-shift XOR)
 *   y ^= state_word_at_lag_L & b                  (T2: lag-mask AND)
 *
 * The lag word is read from the state at bit offset L*w (the L-th
 * w-bit word in the state array).  This works with any generator
 * whose state is an array of w-bit words.
 *
 * Parameters:
 *   w      — word size
 *   sigma  — left shift for T1 (0 < sigma < w)
 *   L      — lag index for T2 (0 < L < num_words)
 *   b      — bit mask for T2
 */
class LaggedTempering : public Transformation {
public:
    LaggedTempering(int w, int sigma, int L, uint64_t b);

    static std::unique_ptr<Transformation> from_params(const Params& params);
    static std::vector<ParamSpec> param_specs();

    std::string name() const override;
    std::string display_str() const override;
    void apply(BitVect& state) const override;
    std::unique_ptr<Transformation> copy() const override;
    void update(const Params& params) override;

private:
    int sigma_;
    int L_;
    uint64_t b_;
};
