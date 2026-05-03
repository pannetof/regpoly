#pragma once
#include "bitvect.h"
#include "params.h"
#include <string>
#include <memory>

class Transformation {
public:
    Transformation() : w_(0) {}
    virtual ~Transformation() = default;

    virtual std::string name() const = 0;
    virtual std::string display_str() const = 0;
    virtual void apply(BitVect& state) const = 0;
    virtual std::unique_ptr<Transformation> copy() const = 0;
    virtual void update(const Params& params) = 0;

    int w() const { return w_; }

protected:
    int w_;
};
