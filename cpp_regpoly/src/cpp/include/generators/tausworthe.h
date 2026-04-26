#pragma once
#include "generator.h"
#include "gen_enumerator.h"
#include "param_spec.h"
#include "params.h"
#include <vector>
#include <memory>
#include <string>
#include <unordered_map>
#include <cstdint>

// Sampled parameter value plus any side-effect params the family's
// sampler wants to splice back into the caller's params bag.  (Used
// today by TauswortheGen's poly sampler, which also chooses `s`.)  Kept
// adjacent to the first family that needs it; hoist to a shared
// header when a second family registers a sampler.
struct RandomParamResult {
    bool is_vec = false;
    int64_t int_val = 0;
    std::vector<int> vec_val;
    std::unordered_map<std::string, int64_t> side_ints;
};

class TauswortheGen : public Generator {
public:
    TauswortheGen(int k, const std::vector<int>& Q, int s, bool quicktaus, int L);

    static std::unique_ptr<Generator> from_params(const Params& params, int L);
    static std::vector<ParamSpec> param_specs();

    // Sample a random admissible polynomial for a TauswortheGen generator
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

    // Handle rand_type values owned by TauswortheGen
    // ("tausworthe_s", "tausworthe_poly").  Reads k, nb_terms,
    // quicktaus, poly, s from `params` as needed.  For
    // "tausworthe_poly" when `s` is unset, picks an admissible s and
    // returns it via side_ints so the caller can propagate that
    // exact value to the generator constructor (the poly was sampled
    // to match it).  Throws std::invalid_argument for unsupported
    // rand_type or invalid inputs.
    static RandomParamResult generate_random(
        const std::string& rand_type,
        const std::string& rand_args,
        const Params& params, int L);

    // True iff the sorted polynomial Q_sorted together with (s,
    // quicktaus, L) satisfies the quicktaus inequality
    //   0 < s <= k - q_{t-2} < k <= L
    // (q_{t-2} = Q_sorted[size-2] = largest interior exponent)
    // or, when quicktaus=false, simply 1 <= s < k.  Also enforces the
    // odd-nb_terms requirement and Q[0] == 0 / Q[-1] == k conventions.
    static bool is_admissible(
        const std::vector<int>& Q_sorted,
        int s, bool quicktaus, int L);

    // Build the exhaustive-search enumerator for this family.  The
    // Params bag carries the structural inputs (k, nb_terms,
    // quicktaus) plus an optional `s` that, if present, fixes the
    // decimation step.  Throws std::invalid_argument with a
    // "needs_<reason>" message when a required axis is missing.
    static std::unique_ptr<GenEnumerator> make_enumerator(
        const Params& resolved, int L);

    std::string name() const override;
    std::string display_str() const override;
    void init(const BitVect& init_bv) override;
    void next() override;
    std::unique_ptr<Generator> copy() const override;
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
