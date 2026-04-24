#pragma once
#include "generateur.h"
#include "param_spec.h"
#include <memory>
#include <string>

class TGFSR : public Generateur {
public:
    TGFSR(int w, int r, int m, const BitVect& a, int L);

    static std::unique_ptr<Generateur> from_params(const Params& params, int L);
    static std::vector<ParamSpec> param_specs();

    std::string name() const override;
    std::string display_str() const override;
    void init(const BitVect& init_bv) override;
    void next() override;
    std::unique_ptr<Generateur> copy() const override;
    BitVect char_poly() const override;

    int w() const { return w_; }
    int r() const { return r_; }
    int m() const { return m_; }
    const BitVect& a() const { return a_; }

private:
    int w_, r_, m_;
    BitVect a_;
};
