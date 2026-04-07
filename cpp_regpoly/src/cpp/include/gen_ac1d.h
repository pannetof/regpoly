#pragma once
#include "generateur.h"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

class AC1DGen : public Generateur {
public:
    AC1DGen(int n, const std::vector<std::vector<int>>& matrix, int L);

    std::string name() const override;
    std::string display_str() const override;
    void init(const BitVect& init_bv) override;
    void next() override;
    std::unique_ptr<Generateur> copy() const override;

private:
    int n_;
    std::vector<std::vector<int>> matrix_;  // n x n binary matrix
};
