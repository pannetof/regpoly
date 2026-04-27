#pragma once
#include "f2w_base.h"
#include "param_spec.h"
#include <memory>
#include <string>
#include <algorithm>

class F2wLFSRGen : public F2wBaseGen {
public:
    F2wLFSRGen(int w, int r, int nbcoeff,
               const std::vector<int>& nocoeff,
               const std::vector<uint64_t>& coeff,
               uint64_t modM, bool normal_basis,
               int step_count, int L);

    static std::unique_ptr<Generator> from_params(const Params& params, int L);
    static std::vector<ParamSpec> param_specs();

    std::string type_name() const override { return "LFSR in F_{2^w}"; }
    void init(const BitVect& init_bv) override;
    void next() override;
    std::unique_ptr<Generator> copy() const override;

    void get_transition_state(uint64_t* out_words, int out_nwords) const override;

    int state_bits() const { return state_bits_; }

private:
    int state_bits_;
};
