#include "poly_lattice.h"
#include <cstring>

// ══════════════════════════════════════════════════════════════════════════
// PolVect
// ══════════════════════════════════════════════════════════════════════════

PolVect::PolVect(int resolution, int degmax)
    : deg(INT_MIN), indicemaxdeg(INT_MIN)
{
    coeffs.reserve(resolution);
    for (int j = 0; j < resolution; j++)
        coeffs.emplace_back(degmax + 1);
}

void PolVect::swap(PolVect& other) {
    coeffs.swap(other.coeffs);
    std::swap(deg, other.deg);
    std::swap(indicemaxdeg, other.indicemaxdeg);
}

// ══════════════════════════════════════════════════════════════════════════
// PolyLatBase construction
// ══════════════════════════════════════════════════════════════════════════

PolyLatBase::PolyLatBase(int maxresolution, int degmax)
    : resolution_(0), maxresolution_(maxresolution), degmax_(degmax)
{
    vect_.reserve(maxresolution);
    for (int i = 0; i < maxresolution; i++)
        vect_.emplace_back(maxresolution, degmax);
    perm_.resize(maxresolution);
    invperm_.resize(maxresolution);
}

// ══════════════════════════════════════════════════════════════════════════
// Element operations
// ══════════════════════════════════════════════════════════════════════════

void PolyLatBase::set_element_zero(PolVect& P, int j) {
    if (P.deg == INT_MIN) return;
    int nw = P.deg / WL;
    uint64_t* d = P.coeffs[j].data();
    for (int i = 0; i <= nw; i++)
        d[i] = 0ULL;
}

void PolyLatBase::set_element_one(PolVect& P, int j) {
    set_element_zero(P, j);
    P.coeffs[j].set_bit(0, 1);
    if (P.deg == INT_MIN) {
        P.deg = 0;
        P.indicemaxdeg = j;
    }
}

int PolyLatBase::norme_bitvect(const BitVect& A) const {
    for (int k = degmax_; k >= 0; k--) {
        if (A.get_bit(k))
            return k;
    }
    return INT_MIN;
}

void PolyLatBase::set_element_from_bv(PolVect& P, int j, const BitVect& A) {
    // Copy A into P.coeffs[j]
    int nw = std::min(A.nwords(), P.coeffs[j].nwords());
    for (int w = 0; w < nw; w++)
        P.coeffs[j].data()[w] = A.data()[w];
    for (int w = nw; w < P.coeffs[j].nwords(); w++)
        P.coeffs[j].data()[w] = 0ULL;

    int deg = norme_bitvect(P.coeffs[j]);
    if (P.deg < deg) {
        P.deg = deg;
        P.indicemaxdeg = j;
    }
}

int PolyLatBase::element_norme_equal(const PolVect& P, int j, int norme) const {
    return P.coeffs[perm_[j]].get_bit(norme) ? 1 : 0;
}

// ══════════════════════════════════════════════════════════════════════════
// update_polvect_deg — find the highest set bit across all coordinates
// ══════════════════════════════════════════════════════════════════════════

void PolyLatBase::update_polvect_deg(PolVect& a) {
    int K = a.deg / WL;
    int LL = a.deg - WL * K;

    for (int k = K; k >= 0; k--) {
        uint64_t blob = 0ULL;
        for (int j = 0; j < resolution_; j++)
            blob |= a.coeffs[j].data()[k];
        for (int ll = LL; ll >= 0; ll--) {
            uint64_t mask = 1ULL << (63 - ll);
            if (blob & mask) {
                a.deg = k * WL + ll;
                for (int j = 0; j < resolution_; j++) {
                    if (a.coeffs[j].data()[k] & mask) {
                        a.indicemaxdeg = j;
                        return;
                    }
                }
            }
        }
        LL = WL - 1;
    }
    a.deg = INT_MIN;
    a.indicemaxdeg = INT_MIN;
}

// ══════════════════════════════════════════════════════════════════════════
// add_polvect — res = a XOR b (polynomial vector addition over GF(2))
// NOTE: modifies b in-place when a.deg > b.deg (masks tail bits)
// ══════════════════════════════════════════════════════════════════════════

void PolyLatBase::add_polvect(PolVect& res, PolVect& a, PolVect& b) {
    if (a.deg == b.deg) {
        res.deg = a.deg;
        int nw = a.deg / WL;
        for (int k = nw; k >= 0; k--)
            for (int j = 0; j < resolution_; j++)
                res.coeffs[j].data()[k] = a.coeffs[j].data()[k] ^ b.coeffs[j].data()[k];
        update_polvect_deg(res);
        return;
    }
    if (a.deg > b.deg) {
        res.deg = a.deg;
        res.indicemaxdeg = a.indicemaxdeg;
        int bK = b.deg / WL;
        int bLL = b.deg - bK * WL;
        for (int j = 0; j < resolution_; j++) {
            for (int k = a.deg / WL; k > bK; k--)
                res.coeffs[j].data()[k] = a.coeffs[j].data()[k];
            // Mask out insignificant bits in b (C code does this in-place)
            b.coeffs[j].data()[bK] &= ~0ULL << (WL - bLL - 1);
            for (int k = bK; k >= 0; k--)
                res.coeffs[j].data()[k] = a.coeffs[j].data()[k] ^ b.coeffs[j].data()[k];
        }
        return;
    }
    // a.deg < b.deg: should not happen in the algorithm's usage
}

// ══════════════════════════════════════════════════════════════════════════
// multbys — r = a * z^S (left shift of polynomial vector by S positions)
//
// In the MSB-first BitVect, multiplying by z^S means the coefficient data
// shifts to higher word indices (the bits representing lower-order coefficients
// are stored at lower indices in MSB-first).  The C code's bit layout is:
//   bit 0 (MSB of word 0) = coefficient of z^0
//   bit k = coefficient of z^k
// Multiplying by z^S: coefficient of z^k in the result comes from
// coefficient of z^(k-S) in the original.  In MSB-first storage this means
// the data words shift to the RIGHT (higher indices) by S/WL words,
// with intra-word right-shift by S%WL.
// ══════════════════════════════════════════════════════════════════════════

void PolyLatBase::multbys(PolVect& r, const PolVect& a, int S) {
    int N = (a.deg + S) / WL;
    int nbblocks = S / WL;
    int nbbits = S - nbblocks * WL;

    if (nbbits != 0) {
        int WLmn = WL - nbbits;
        for (int ii = 0; ii < resolution_; ii++) {
            uint64_t* rd = r.coeffs[ii].data();
            const uint64_t* ad = a.coeffs[ii].data();
            for (int i = 0; i < nbblocks; i++)
                rd[i] = 0ULL;
            int i = nbblocks;
            rd[i++] = ad[0] >> nbbits;
            for (; i <= N; i++)
                rd[i] = (ad[i - nbblocks] >> nbbits) | (ad[i - nbblocks - 1] << WLmn);
        }
    } else {
        for (int ii = 0; ii < resolution_; ii++) {
            uint64_t* rd = r.coeffs[ii].data();
            const uint64_t* ad = a.coeffs[ii].data();
            for (int i = 0; i < nbblocks; i++)
                rd[i] = 0ULL;
            int i = nbblocks;
            rd[i++] = ad[0];
            for (; i <= N; i++)
                rd[i] = ad[i - nbblocks];
        }
    }
    r.deg = a.deg + S;
    r.indicemaxdeg = a.indicemaxdeg;
}

// ══════════════════════════════════════════════════════════════════════════
// multbys2 — r XOR= a * z^S in-place
// ══════════════════════════════════════════════════════════════════════════

void PolyLatBase::multbys2(PolVect& r, const PolVect& a, int S) {
    int nbblocks = S / WL;
    int nbbits = S - nbblocks * WL;

    if (a.deg + S == r.deg) {
        int N = r.deg / WL;
        if (nbbits != 0) {
            int WLmn = WL - nbbits;
            for (int ii = 0; ii < resolution_; ii++) {
                uint64_t* rd = r.coeffs[ii].data();
                const uint64_t* ad = a.coeffs[ii].data();
                int i = nbblocks;
                rd[i++] ^= ad[0] >> nbbits;
                for (; i <= N; i++)
                    rd[i] ^= (ad[i - nbblocks] >> nbbits) | (ad[i - nbblocks - 1] << WLmn);
            }
        } else {
            for (int ii = 0; ii < resolution_; ii++) {
                uint64_t* rd = r.coeffs[ii].data();
                const uint64_t* ad = a.coeffs[ii].data();
                int i = nbblocks;
                rd[i++] ^= ad[0];
                for (; i <= N; i++)
                    rd[i] ^= ad[i - nbblocks];
            }
        }
        r.indicemaxdeg = a.indicemaxdeg;
        update_polvect_deg(r);
        return;
    }

    if (a.deg + S < r.deg) {
        int N = (a.deg + S) / WL;
        if (nbbits != 0) {
            int WLmn = WL - nbbits;
            for (int ii = 0; ii < resolution_; ii++) {
                uint64_t* rd = r.coeffs[ii].data();
                const uint64_t* ad = a.coeffs[ii].data();
                int i = nbblocks;
                rd[i++] ^= ad[0] >> nbbits;
                for (; i < N; i++)
                    rd[i] ^= (ad[i - nbblocks] >> nbbits) | (ad[i - nbblocks - 1] << WLmn);
                rd[i] ^= (ad[i - nbblocks - 1] << WLmn);
            }
        } else {
            for (int ii = 0; ii < resolution_; ii++) {
                uint64_t* rd = r.coeffs[ii].data();
                const uint64_t* ad = a.coeffs[ii].data();
                int i = nbblocks;
                rd[i++] ^= ad[0];
                for (; i <= N; i++)
                    rd[i] ^= ad[i - nbblocks];
            }
        }
        r.deg = a.deg + S;
        r.indicemaxdeg = a.indicemaxdeg;
        update_polvect_deg(r);
        return;
    }
}

void PolyLatBase::add_self_polvect_mult(PolVect& res, const PolVect& a, int s) {
    if (res.deg == INT_MIN)
        multbys(res, a, s);
    else
        multbys2(res, a, s);
}

// ══════════════════════════════════════════════════════════════════════════
// DualBase — initialize the dual basis
// ══════════════════════════════════════════════════════════════════════════

void PolyLatBase::dual_base(const std::vector<BitVect>& polys, const BitVect& M, int res) {
    for (int i = 0; i < maxresolution_; i++) {
        perm_[i] = i;
        invperm_[i] = i;
    }
    resolution_ = res;

    // Row 0: the characteristic polynomial M
    vb(0).zero_lazy();
    set_element_from_bv(vb(0), 0, M);

    // Rows 1..res-1: normalized polynomial h_i + identity in coordinate i
    for (int i = 1; i < res; i++) {
        vb(i).zero_lazy();
        set_element_from_bv(vb(i), 0, polys[i]);
        set_element_one(vb(i), i);
    }
}

void PolyLatBase::dual_base_increase(const std::vector<BitVect>& polys) {
    resolution_++;
    int newi = resolution_ - 1;
    // Zero out the new coordinate in all existing rows
    for (int j = 0; j <= newi - 1; j++)
        set_element_zero(vb(j), newi);
    // New row: normalized polynomial + identity
    vb(newi).zero_lazy();
    set_element_from_bv(vb(newi), 0, polys[newi]);
    set_element_one(vb(newi), newi);
}

// ══════════════════════════════════════════════════════════════════════════
// Lenstra helpers
// ══════════════════════════════════════════════════════════════════════════

void PolyLatBase::renumber(int m, int resolution) {
    int min_val = INT_MAX;
    int minj = m + 1;
    for (int j = m + 1; j <= resolution; j++) {
        int nj = vb(j).deg;
        if (nj < min_val) {
            min_val = nj;
            minj = j;
        }
    }
    if (minj != m + 1)
        vect_[m + 1].swap(vect_[minj]);
}

void PolyLatBase::solve_axb(int m, std::vector<int>& x) {
    // Small system: (m+1) x (m+2) bits — use a simple dense array
    int n = m + 1;
    int nw = (n + 1 + 63) / 64;
    std::vector<std::vector<uint64_t>> A(n, std::vector<uint64_t>(nw, 0ULL));

    auto setbit = [&](int row, int col) {
        A[row][col / 64] |= 1ULL << (63 - (col % 64));
    };
    auto getbit = [&](int row, int col) -> int {
        return (A[row][col / 64] >> (63 - (col % 64))) & 1;
    };

    // Fill upper triangle: A[j][i] = element_norme_equal(vect[i], j, norm(vect[i]))
    for (int i = 0; i <= m; i++) {
        for (int j = i; j <= m; j++) {
            if (element_norme_equal(vect_[i], j, vect_[i].deg))
                setbit(j, i);
        }
    }
    // RHS column (m+1): from vect[m+1]
    int bmp1_norme = vect_[m + 1].deg;
    for (int j = 0; j <= m; j++) {
        if (element_norme_equal(vect_[m + 1], j, bmp1_norme))
            setbit(j, m + 1);
    }

    // Forward substitution (system is already upper triangular in C code)
    for (int j = 0; j <= m; j++) {
        int sum = getbit(j, m + 1);
        for (int i = 0; i < j; i++)
            sum += getbit(j, i) * x[i];
        x[j] = sum & 1;
    }
}

void PolyLatBase::permute_coord(int m) {
    int j = invperm_[vb(m + 1).indicemaxdeg];
    if (j != (m + 1)) {
        std::swap(perm_[m + 1], perm_[j]);
        std::swap(invperm_[perm_[j]], invperm_[perm_[m + 1]]);
    }
}

// ══════════════════════════════════════════════════════════════════════════
// Lenstra's algorithm
// ══════════════════════════════════════════════════════════════════════════

int PolyLatBase::lenstra(int res) {
    int resolution = res - 1;
    int m = -1;
    std::vector<int> x;
    if (resolution > 0)
        x.resize(resolution, 0);

    PolVect temp(maxresolution_, degmax_);
    PolVect sum(maxresolution_, degmax_);

    while (m < resolution) {
        renumber(m, resolution);
        if (m > -1)
            solve_axb(m, x);

        int Normebmp1 = vb(m + 1).deg;
        sum.zero_lazy();

        for (int i = 0; i <= m; i++) {
            if (x[i]) {
                int Normebi = vb(i).deg;
                int s = Normebmp1 - Normebi;
                add_self_polvect_mult(sum, vb(i), s);
            }
        }

        int swap_flag;
        if (sum.deg == INT_MIN) {
            temp.swap(vb(m + 1));
            swap_flag = 1;
        } else {
            add_polvect(temp, vb(m + 1), sum);
            swap_flag = 0;
        }

        int NormeTbmp1 = temp.deg;
        if (NormeTbmp1 == Normebmp1) {
            vb(m + 1).swap(temp);
            permute_coord(m);
            m++;
        } else if (NormeTbmp1 < Normebmp1) {
            vb(m + 1).swap(temp);
            bool OK = false;
            int newm = 0;
            for (int L = 0; L <= m; L++) {
                if (vb(L).deg <= NormeTbmp1) {
                    newm = L;
                    OK = true;
                }
            }
            if (!OK)
                m = -1;
            else
                m = newm;
        } else if (swap_flag == 1) {
            vb(m + 1).swap(temp);
        }
    }

    return vb(0).deg;
}
