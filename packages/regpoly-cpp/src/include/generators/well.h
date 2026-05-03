#pragma once
#include "generator.h"
#include "param_spec.h"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <functional>

class WELLGen : public Generator {
public:
    struct MatrixEntry {
        int type;
        int paramsint[3];
        uint64_t paramsulong[3];
    };

    WELLGen(int w, int r, int p, int m1, int m2, int m3,
              const std::vector<MatrixEntry>& matrices, int L);

    static std::unique_ptr<Generator> from_params(const Params& params, int L);
    static std::vector<ParamSpec> param_specs();

    std::string name() const override;
    std::string display_str() const override;
    void init(const BitVect& init_bv) override;
    void next() override;
    std::unique_ptr<Generator> copy() const override;
    BitVect get_output() const override;
    void get_transition_state(uint64_t* out_words, int out_nwords) const override;

private:
    int w_;
    int r_;
    int p_;         // output mask bits
    int m1_, m2_, m3_;
    std::vector<MatrixEntry> matrices_;  // 8 entries
    int i_;         // circular buffer pointer
    int state_bits_; // w * r (full state size for circular buffer)
    uint64_t maskp_;   // UPPER mask: p least significant bits
    uint64_t umaskp_;  // LOWER mask: ~maskp
    static constexpr uint64_t M32 = 0xFFFFFFFFULL;

    static uint64_t ShiftR(uint64_t v, int s);
    static uint64_t apply_matrix(int type, uint64_t v, const int* pi, const uint64_t* pu);
    uint64_t TMAT(int j, uint64_t val) const;
    uint64_t V(int idx) const;
    void SetV(int idx, uint64_t val);

    static int type_cost(int type);
    static std::string type_display(int type, const int* pi, const uint64_t* pu);
};
