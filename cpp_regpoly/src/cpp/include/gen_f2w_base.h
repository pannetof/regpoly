#pragma once
#include "generateur.h"
#include "gf2w_arith.h"
#include <vector>
#include <cstdint>

enum GenF2wType { GENF2W_POLYLCG = 0, GENF2W_LFSR = 1 };

class GenF2wBase : public Generateur {
protected:
    int w_, r_;
    int nbcoeff_;
    std::vector<int> nocoeff_;
    std::vector<uint64_t> coeff_;
    uint64_t modM_;
    bool normal_basis_;
    std::vector<uint64_t> table_;
    uint64_t maskw_;
    int step_count_;

public:
    GenF2wBase(int w, int r, int nbcoeff,
               const std::vector<int>& nocoeff,
               const std::vector<uint64_t>& coeff,
               uint64_t modM, bool normal_basis,
               int step_count, int L);

    // Pure virtual — GenF2wBase cannot be instantiated directly
    std::string name() const override { return "Generator in F_{2^w}"; }
    void init(const BitVect& init_bv) override = 0;
    void next() override = 0;
    std::unique_ptr<Generateur> copy() const override = 0;

    std::string display_str() const override;
    virtual std::string type_name() const = 0;

    int gf2w_w() const { return w_; }
    int gf2w_r() const { return r_; }
    int nbcoeff() const { return nbcoeff_; }
    const std::vector<int>& nocoeff() const { return nocoeff_; }
    const std::vector<uint64_t>& coeff() const { return coeff_; }
    uint64_t modM() const { return modM_; }
    bool normal_basis() const { return normal_basis_; }
    int step_count() const { return step_count_; }
    const std::vector<uint64_t>& table() const { return table_; }

protected:
    uint64_t multiply(uint64_t a, uint64_t b) const;
    uint64_t V(int idx) const;
    void SetV(int idx, uint64_t val);
};
