#pragma once

#include <memory>
#include <string>
#include <vector>

#include "params.h"

#include <NTL/ZZ.h>

// ─── Generator-space enumerator ────────────────────────────────────────────
//
// Abstract interface used by the exhaustive-search mode.  Concrete
// enumerators live next to the generator class they enumerate (for
// example TausworthePolyEnumerator in gen_tausworthe.cpp).  The
// worker and API only ever talk to the base class.
//
// Counts are exposed as decimal strings because k-dependent
// combinatorial totals routinely exceed 2^63.  Python bindings lift
// them to arbitrary-precision ints.

class GenEnumerator {
public:
    struct Axis {
        std::string name;
        std::string size_dec;    // decimal string (may exceed 2^63)
        std::string describe;    // human-readable hint for the UI
    };

    virtual ~GenEnumerator() = default;

    // Total number of admissible combinations.  Decimal representation.
    virtual std::string size_dec() const = 0;

    // Concrete parameters for the idx-th combination (idx as decimal
    // string).  Must be admissible by construction.  Deterministic and
    // stable across processes — same inputs produce the same Params.
    virtual Params at(const std::string& idx_dec) const = 0;

    // Per-axis metadata for the UI summary; outer-most axis first.
    virtual std::vector<Axis> axes() const = 0;
};

// ─── Registry (implemented in factory.cpp) ────────────────────────────────

// Returns nullptr when the family has no exhaustive enumerator.  Throws
// std::invalid_argument with a "needs_*" reason when the specific
// `resolved` inputs are insufficient (e.g. Tausworthe without
// nb_terms).
std::unique_ptr<GenEnumerator> make_gen_enumerator(
    const std::string& family, const Params& resolved, int L);

bool family_is_enumerable(const std::string& family);

// ─── Reusable primitives for future enumerators ───────────────────────────
//
// unrank_combination(n, k, idx) returns the idx-th k-subset of
// {0, 1, ..., n-1} in lexicographic order, via the combinatorial
// number system.  idx is accepted as an NTL::ZZ so callers can pass
// large offsets without truncation.
std::vector<int> unrank_combination(int n, int k, const NTL::ZZ& idx);

// Mixed-radix decode: given N axis sizes and a flat idx, returns
// per-axis indices.  The first axis varies slowest (sizes[0] is the
// outer digit).
std::vector<NTL::ZZ> mixed_radix_decode(
    const std::vector<NTL::ZZ>& sizes, const NTL::ZZ& idx);
