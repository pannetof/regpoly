#include "trans_lag_mask.h"
#include <sstream>
#include <iomanip>

LaggedTempering::LaggedTempering(int w, int sigma, int L, uint64_t b)
    : sigma_(sigma), L_(L), b_(b)
{
    w_ = w;
}

std::string LaggedTempering::name() const {
    return "Lagged Tempering";
}

std::string LaggedTempering::display_str() const {
    int hex_digits = (w_ + 3) / 4;
    std::ostringstream oss;
    oss << "Lagged Tempering (sigma=" << sigma_
        << ", L=" << L_ << ", w=" << w_ << ")  b = "
        << std::hex << std::setfill('0') << std::setw(hex_digits) << b_;
    return oss.str();
}

void LaggedTempering::apply(BitVect& state) const {
    int n = state.nbits();

    // --- T1: y = output ^ (output << sigma) ---
    // Extract the first w bits (output word)
    BitVect out_word = state.copy();
    out_word.and_mask(w_);

    BitVect shifted = out_word.copy();
    shifted.lshift(sigma_);
    out_word.xor_with(shifted);
    out_word.and_mask(w_);

    // --- T2: y ^= (state_word_at_lag_L & b) ---
    // Read w bits starting at position L*w in the state
    int lag_offset = L_ * w_;
    if (lag_offset + w_ <= n) {
        // Build mask b as a w-bit BitVect at position 0
        BitVect b_mask(n);
        for (int i = 0; i < w_; i++) {
            if ((b_ >> (w_ - 1 - i)) & 1ULL)
                b_mask.set_bit(i, 1);
        }

        // Read the lag word and place it at position 0
        BitVect lag_word(n);
        for (int i = 0; i < w_; i++) {
            int src = lag_offset + i;
            if (src < n && state.get_bit(src))
                lag_word.set_bit(i, 1);
        }

        // AND with mask b
        for (int i = 0; i < lag_word.nwords(); i++)
            lag_word.data()[i] &= b_mask.data()[i];

        // XOR into output
        out_word.xor_with(lag_word);
        out_word.and_mask(w_);
    }

    // --- Recombine: replace first w bits, keep the rest ---
    BitVect rest = state.copy();
    rest.and_invmask(w_);
    rest.xor_with(out_word);
    state = rest;
}

std::unique_ptr<Transformation> LaggedTempering::copy() const {
    return std::make_unique<LaggedTempering>(w_, sigma_, L_, b_);
}

void LaggedTempering::update(const Params& params) {
    w_ = (int)params.get_int("w", w_);
    sigma_ = (int)params.get_int("sigma", sigma_);
    L_ = (int)params.get_int("L", L_);
    b_ = (uint64_t)params.get_int("b", (int64_t)b_);
}

// ── Factory method ─────────────────────────────────────────────────────

std::unique_ptr<Transformation> LaggedTempering::from_params(const Params& params) {
    int w = (int)params.get_int("w");
    int sigma = (int)params.get_int("sigma");
    int L = (int)params.get_int("L");
    uint64_t b = (uint64_t)params.get_int("b", 0);
    return std::make_unique<LaggedTempering>(w, sigma, L, b);
}

std::vector<ParamSpec> LaggedTempering::param_specs() {
    return {
        {"w",     "int", true,  false, 0, "",        "", false},
        {"sigma", "int", false, false, 0, "range",   "1,w-1", false},
        {"L",     "int", false, false, 0, "",        "", false},
        {"b",     "int", false, false, 0, "bitmask", "w", true},
    };
}
