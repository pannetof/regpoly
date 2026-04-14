#pragma once
#include "bitvect.h"
#include <vector>
#include <climits>
#include <algorithm>

// PolVect: a basis vector = an l-tuple of polynomials over GF(2).
//
// Storage: a single flat uint64_t array of size (resolution * nwords),
// laid out as [coord_0_word_0, ..., coord_0_word_{nw-1},
//              coord_1_word_0, ..., coord_{res-1}_word_{nw-1}].
// This avoids per-coordinate heap allocations and improves cache behavior.
//
// deg  = max degree across all coordinates (INT_MIN if zero vector).
// indicemaxdeg = which coordinate achieves the max degree.

struct PolVect {
    std::vector<uint64_t> data;   // flat storage: resolution * nwords
    int deg = INT_MIN;
    int indicemaxdeg = INT_MIN;
    int nwords = 0;               // words per coordinate
    int resolution = 0;           // number of coordinates

    PolVect() = default;
    PolVect(int resolution, int degmax);

    // Pointer to word array for coordinate j
    uint64_t*       coord(int j)       { return data.data() + j * nwords; }
    const uint64_t* coord(int j) const { return data.data() + j * nwords; }

    void swap(PolVect& other);
    void zero_lazy() { deg = INT_MIN; indicemaxdeg = INT_MIN; }
};

// PolyLatBase: polynomial lattice basis for the dual lattice method.
class PolyLatBase {
public:
    PolyLatBase(int maxresolution, int degmax);

    void dual_base(const std::vector<BitVect>& polys, const BitVect& M, int res);
    void dual_base_increase(const std::vector<BitVect>& polys);
    int lenstra(int res);

    int resolution() const { return resolution_; }
    int degmax() const { return degmax_; }

private:
    std::vector<PolVect> vect_;
    std::vector<int> perm_;
    std::vector<int> invperm_;
    int resolution_;
    int maxresolution_;
    int degmax_;
    int nwords_;                  // words per coordinate = ceil((degmax+1)/64)

    static constexpr int WL = 64;

    PolVect& vb(int i) { return vect_[i]; }
    const PolVect& vb(int i) const { return vect_[i]; }

    // Element operations
    void set_element_zero(PolVect& P, int j);
    void set_element_one(PolVect& P, int j);
    void set_element_from_bv(PolVect& P, int j, const BitVect& A);
    int element_norme_equal(const PolVect& P, int j, int norme) const;

    // PolVect operations
    void update_polvect_deg(PolVect& a);
    void add_polvect(PolVect& res, PolVect& a, PolVect& b);
    void multbys(PolVect& r, const PolVect& a, int S);
    void multbys2(PolVect& r, const PolVect& a, int S);
    void add_self_polvect_mult(PolVect& res, const PolVect& a, int s);

    // Lenstra helpers
    void renumber(int m, int resolution);
    void solve_axb(int m, std::vector<int>& x);
    void permute_coord(int m);
};
