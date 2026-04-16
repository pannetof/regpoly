#pragma once
#include "gen_f2w_base.h"
#include "param_spec.h"
#include <memory>
#include <string>

class GenF2wPolyLCG : public GenF2wBase {
public:
    GenF2wPolyLCG(int w, int r, int nbcoeff,
                  const std::vector<int>& nocoeff,
                  const std::vector<uint64_t>& coeff,
                  uint64_t modM, bool normal_basis,
                  int step_count, int L);

    static std::unique_ptr<Generateur> from_params(const Params& params, int L);
    static std::vector<ParamSpec> param_specs();

    std::string type_name() const override { return "Polynomial LCG in F_{2^w}[z]/P(z)"; }
    void init(const BitVect& init_bv) override;
    void next() override;
    std::unique_ptr<Generateur> copy() const override;
};
