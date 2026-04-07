#pragma once
#include "bitvect.h"
#include <memory>
#include <string>
#include <sstream>
#include <vector>
class Generateur {
public:
    Generateur(int k, int L) : state_(k), k_(k), L_(L) {}
    virtual ~Generateur() = default;

    // Pure virtual interface — makes Generateur abstract
    virtual std::string name() const = 0;
    virtual std::string display_str() const = 0;
    virtual void init(const BitVect& init_bv) = 0;
    virtual void next() = 0;
    virtual std::unique_ptr<Generateur> copy() const = 0;

    // Public accessors
    int k() const { return k_; }
    int L() const { return L_; }
    const BitVect& state() const { return state_; }

    // Overridable behaviour
    virtual BitVect get_output() const;
    virtual void get_transition_state(uint64_t* out_words, int out_nwords) const;

    // Algorithms (use the virtual interface above)
    BitVect char_poly() const;

    // Compute the K×K transition matrix.
    // Returns K row BitVects, where row[i] has bit j set if A[i][j] = 1.
    std::vector<BitVect> transition_matrix() const;

protected:
    BitVect state_;
    int k_;
    int L_;
};
