#pragma once
#include "generateur.h"
#include "param_spec.h"
#include <vector>
#include <memory>
#include <string>

class Tausworthe : public Generateur {
public:
    Tausworthe(int k, const std::vector<int>& Q, int s, bool quicktaus, int L);

    static std::unique_ptr<Generateur> from_params(const Params& params, int L);
    static std::vector<ParamSpec> param_specs();

    // Sample a random admissible polynomial for a Tausworthe generator
    // with the given shape.  Returns the sorted exponent list
    // [0, q_1, ..., q_{t-2}, k].  Throws std::invalid_argument when the
    // combination is inadmissible (nb_terms even or < 3, k <= 0, or
    // quicktaus=true with k > L).
    //
    // `s`: if 0, the caller wants auto-derivation (s = k - q_{t-2}) —
    // the top interior exponent may be anywhere in [1, k-1].
    // If > 0, the poly must satisfy q_{t-2} <= k - s.
    static std::vector<int> random_poly(
        int k, int nb_terms, bool quicktaus, int L, int s);

    // True iff the sorted polynomial Q_sorted together with (s,
    // quicktaus, L) satisfies the quicktaus inequality
    //   0 < s <= k - q_{t-2} < k <= L
    // (q_{t-2} = Q_sorted[size-2] = largest interior exponent)
    // or, when quicktaus=false, simply 1 <= s < k.  Also enforces the
    // odd-nb_terms requirement and Q[0] == 0 / Q[-1] == k conventions.
    static bool is_admissible(
        const std::vector<int>& Q_sorted,
        int s, bool quicktaus, int L);

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
