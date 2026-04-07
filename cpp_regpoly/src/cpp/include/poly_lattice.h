#pragma once
#include "bitvect.h"
#include <vector>
#include <climits>
#include <algorithm>

// PolVect: a vector of polynomials over GF(2), one per coordinate.
// Each coeffs[j] is a BitVect representing the polynomial for coordinate j.
// deg = max degree across all coordinates (INT_MIN if zero vector).
// indicemaxdeg = coordinate index achieving the max degree.
struct PolVect {
    std::vector<BitVect> coeffs;
    int deg = INT_MIN;
    int indicemaxdeg = INT_MIN;

    PolVect() = default;
    PolVect(int resolution, int degmax);

    void zero_lazy() { deg = INT_MIN; indicemaxdeg = INT_MIN; }
    void swap(PolVect& other);
};

// PolyLatBase: polynomial lattice basis for the dual lattice method.
// Wraps the Lenstra algorithm and dual basis construction.
class PolyLatBase {
public:
    PolyLatBase(int maxresolution, int degmax);

    // Initialize dual basis with char poly M at row 0, and polys at rows 1..res-1.
    void dual_base(const std::vector<BitVect>& polys, const BitVect& M, int res);

    // Grow the basis by one row, adding the next normalized polynomial.
    void dual_base_increase(const std::vector<BitVect>& polys);

    // Run Lenstra's lattice reduction up to resolution res.
    // Returns the first successive minimum (norm of shortest vector).
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

    static constexpr int WL = 64;

    // Accessors matching C macros
    PolVect& vb(int i) { return vect_[i]; }
    const PolVect& vb(int i) const { return vect_[i]; }
    int vperm(int i) const { return perm_[i]; }

    // Element operations
    void set_element_zero(PolVect& P, int j);
    void set_element_one(PolVect& P, int j);
    void set_element_from_bv(PolVect& P, int j, const BitVect& A);
    int norme_bitvect(const BitVect& A) const;
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
