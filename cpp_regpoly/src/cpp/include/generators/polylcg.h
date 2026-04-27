#pragma once
#include "generator.h"
#include "param_spec.h"
#include <memory>
#include <string>

class PolyLCGGen : public Generator {
public:
    PolyLCGGen(int k, const BitVect& poly, int L);

    static std::unique_ptr<Generator> from_params(const Params& params, int L);
    static std::vector<ParamSpec> param_specs();

    std::string name() const override;
    std::string display_str() const override;
    void init(const BitVect& init_bv) override;
    void next() override;
    std::unique_ptr<Generator> copy() const override;
    BitVect char_poly() const override;

    const BitVect& poly() const { return poly_; }

private:
    BitVect poly_;
};
