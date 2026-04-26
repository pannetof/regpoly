#pragma once
#include "generator.h"
#include "param_spec.h"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

class MarsaXorshiftGen : public Generator {
public:
    // Type 1: single word xorshift
    struct Type1Params {
        int a, b, c;
    };

    // Type 2x: two-component
    struct Type2xParams {
        std::vector<int> p;  // 3 values
        std::vector<int> q;  // 3 values
    };

    // Type 3: multi-tap
    struct Tap {
        int position;
        int shift;
    };

    // Type 4: two-component, 2 shifts each
    struct Type4Params {
        std::vector<int> p;  // 2 values
        std::vector<int> q;  // 2 values
    };

    // Type 100 (general): multi-component
    struct MiEntry {
        int position;
        std::vector<int> shifts;
    };

    MarsaXorshiftGen(int type, int w, int r, int m,
                     const Type1Params& t1,
                     const Type2xParams& t2x,
                     const std::vector<Tap>& taps,
                     const Type4Params& t4,
                     const std::vector<MiEntry>& mi,
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
    int w_;
    int r_;
    int m_;
    Type1Params t1_;
    Type2xParams t2x_;
    std::vector<Tap> taps_;
    Type4Params t4_;
    std::vector<MiEntry> mi_;

    static uint64_t ShiftR(uint64_t v, int s);
    uint64_t V(int idx) const;
    void SetV(int idx, uint64_t val);
};
