#pragma once
#include "gen_f2w_base.h"
#include <memory>
#include <string>
#include <algorithm>

class GenF2wLFSR : public GenF2wBase {
public:
    GenF2wLFSR(int w, int r, int nbcoeff,
               const std::vector<int>& nocoeff,
               const std::vector<uint64_t>& coeff,
               uint64_t modM, bool normal_basis,
               int step_count, int L);

    std::string type_name() const override { return "LFSR in F_{2^w}"; }
    void init(const BitVect& init_bv) override;
    void next() override;
    std::unique_ptr<Generateur> copy() const override;

    void get_transition_state(uint64_t* out_words, int out_nwords) const override;

    int state_bits() const { return state_bits_; }

private:
    int state_bits_;
};
