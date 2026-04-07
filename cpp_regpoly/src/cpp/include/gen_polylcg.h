#pragma once
#include "generateur.h"
#include <memory>
#include <string>

class PolyLCG : public Generateur {
public:
    PolyLCG(int k, const BitVect& poly, int L);

    std::string name() const override;
    std::string display_str() const override;
    void init(const BitVect& init_bv) override;
    void next() override;
    std::unique_ptr<Generateur> copy() const override;

    const BitVect& poly() const { return poly_; }

private:
    BitVect poly_;
};
