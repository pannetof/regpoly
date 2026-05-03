#pragma once
#include "transformation.h"
#include "param_spec.h"
#include <cstdint>
#include <memory>
#include <string>

class TemperMKTrans : public Transformation {
public:
    TemperMKTrans(int w, int type, int eta, int mu, int u, int l,
                  uint64_t b, uint64_t c);

    static std::unique_ptr<Transformation> from_params(
        const std::string& type_name, const Params& params);
    static std::vector<ParamSpec> param_specs();

    std::string name() const override;
    std::string display_str() const override;
    void apply(BitVect& state) const override;
    std::unique_ptr<Transformation> copy() const override;
    void update(const Params& params) override;

    int type() const { return type_; }
    int eta() const { return eta_; }
    int mu() const { return mu_; }
    int u() const { return u_; }
    int l() const { return l_; }
    uint64_t b() const { return b_; }
    uint64_t c() const { return c_; }

private:
    int type_;
    int eta_, mu_;
    int u_, l_;
    uint64_t b_, c_;
};
