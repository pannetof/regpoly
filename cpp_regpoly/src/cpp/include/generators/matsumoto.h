#pragma once
#include "generator.h"
#include "param_spec.h"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

class MatsumotoGen : public Generator {
public:
    MatsumotoGen(int type, int n, int m,
                 const std::vector<int>& paramsint,
                 const std::vector<uint32_t>& paramsunsigned,
                 int L);

    static std::unique_ptr<Generator> from_params(const Params& params, int L);
    static std::vector<ParamSpec> param_specs();

    std::string name() const override;
    std::string display_str() const override;
    void init(const BitVect& init_bv) override;
    void next() override;
    std::unique_ptr<Generator> copy() const override;

private:
    int type_;
    int n_;
    int m_;
    std::vector<int> p_;          // shift parameters
    std::vector<uint32_t> pu_;    // unsigned parameters

    static uint32_t SHIFT(uint32_t v, int i);
    uint32_t V(int idx) const;
    void SetV(int idx, uint32_t val);
};
