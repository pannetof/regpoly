#pragma once
#include "generateur.h"
#include <vector>
#include <memory>
#include <string>

class Tausworthe : public Generateur {
public:
    Tausworthe(int k, const std::vector<int>& Q, int s, bool quicktaus, int L);

    std::string name() const override;
    std::string display_str() const override;
    void init(const BitVect& init_bv) override;
    void next() override;
    std::unique_ptr<Generateur> copy() const override;
    BitVect char_poly() const override;
    void get_transition_state(uint64_t* out_words, int out_nwords) const override;

    int s() const { return s_; }
    bool quicktaus() const { return quicktaus_; }
    const std::vector<int>& Q() const { return Q_; }

private:
    std::vector<int> Q_;
    int NbCoeff_;
    int s_;
    bool quicktaus_;
    int gen_kms_;

    void next_quick();
    void next_general();
};
