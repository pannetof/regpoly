#pragma once
#include "transformation.h"
#include <memory>
#include <string>

class PermutationTrans : public Transformation {
public:
    PermutationTrans(int w, int p, int q);

    std::string name() const override;
    std::string display_str() const override;
    void apply(BitVect& state) const override;
    std::unique_ptr<Transformation> copy() const override;
    void update(const Params& params) override;

    int p() const { return p_; }
    int q() const { return q_; }

private:
    int p_, q_;
};
