#pragma once
#include "generateur.h"
#include "param_spec.h"
#include <cstdint>
#include <memory>
#include <string>

class MersenneTwister : public Generateur {
public:
    MersenneTwister(int w, int r, int m, int p, uint64_t a, int L);

    static std::unique_ptr<Generateur> from_params(const Params& params, int L);
    static std::vector<ParamSpec> param_specs();

    std::string name() const override;
    std::string display_str() const override;
    void init(const BitVect& init_bv) override;
    void next() override;
    std::unique_ptr<Generateur> copy() const override;
    BitVect char_poly() const override;
    BitVect get_output() const override;
    void get_transition_state(uint64_t* out_words, int out_nwords) const override;

    int w() const { return w_; }
    int r() const { return r_; }
    int m() const { return m_; }
    int p() const { return p_; }
    uint64_t a() const { return a_; }
    int i_val() const { return i_; }

private:
    int w_, r_, m_, p_;
    uint64_t a_;
    int i_;
    uint64_t uu_, ll_;
    uint64_t maskw_;
    int state_bits_;

    uint64_t V(int idx) const;
    void SetV(int idx, uint64_t val);
    BitVect rotated_state() const;
};
