#include "poly_lattice.h"
#include <cstring>

// ══════════════════════════════════════════════════════════════════════════
// PolVect — flat contiguous storage for l polynomials of degree ≤ degmax
// ══════════════════════════════════════════════════════════════════════════

PolVect::PolVect(int resolution_, int degmax)
    : deg(INT_MIN), indicemaxdeg(INT_MIN),
      nwords((degmax + 1 + 63) / 64), resolution(resolution_)
{
    data.resize((size_t)resolution_ * nwords, 0ULL);
}

void PolVect::swap(PolVect& other) {
    data.swap(other.data);
    std::swap(deg, other.deg);
    std::swap(indicemaxdeg, other.indicemaxdeg);
    std::swap(nwords, other.nwords);
    std::swap(resolution, other.resolution);
}

// ══════════════════════════════════════════════════════════════════════════
// PolyLatBase construction
// ══════════════════════════════════════════════════════════════════════════

PolyLatBase::PolyLatBase(int maxresolution, int degmax)
    : resolution_(0), maxresolution_(maxresolution), degmax_(degmax),
      nwords_((degmax + 1 + 63) / 64)
{
    vect_.reserve(maxresolution);
    for (int i = 0; i < maxresolution; i++)
        vect_.emplace_back(maxresolution, degmax);
    perm_.resize(maxresolution);
    invperm_.resize(maxresolution);
}

// ══════════════════════════════════════════════════════════════════════════
// Element operations — identical logic to original, using flat coord()
// ══════════════════════════════════════════════════════════════════════════

void PolyLatBase::set_element_zero(PolVect& P, int j) {
    if (P.deg == INT_MIN) return;
    uint64_t* d = P.coord(j);
    int nw = P.deg / WL;
    for (int i = 0; i <= nw; i++)
        d[i] = 0ULL;
}

void PolyLatBase::set_element_one(PolVect& P, int j) {
    set_element_zero(P, j);
    P.coord(j)[0] = 1ULL << 63;  // bit 0 = MSB = z^0
    if (P.deg == INT_MIN) {
        P.deg = 0;
        P.indicemaxdeg = j;
    }
}

void PolyLatBase::set_element_from_bv(PolVect& P, int j, const BitVect& A) {
    uint64_t* d = P.coord(j);
    int nw = std::min(A.nwords(), P.nwords);
    for (int w = 0; w < nw; w++)
        d[w] = A.data()[w];
    for (int w = nw; w < P.nwords; w++)
        d[w] = 0ULL;

    // Find highest degree in this coordinate.
    // MSB-first: lowest bit position = highest polynomial degree.
    for (int k = nw - 1; k >= 0; k--) {
        if (d[k]) {
            int ctz = __builtin_ctzll(d[k]);
            int deg_j = k * WL + (WL - 1 - ctz);
            if (P.deg < deg_j) {
                P.deg = deg_j;
                P.indicemaxdeg = j;
            }
            return;
        }
    }
}

int PolyLatBase::element_norme_equal(const PolVect& P, int j, int norme) const {
    const uint64_t* d = P.coord(perm_[j]);
    int w = norme / WL;
    int b = 63 - (norme % WL);
    return (d[w] >> b) & 1;
}

// ══════════════════════════════════════════════════════════════════════════
// update_polvect_deg — uses __builtin_clzll for fast scanning
// ══════════════════════════════════════════════════════════════════════════

void PolyLatBase::update_polvect_deg(PolVect& a) {
    // MSB-first: within word k, bit position p = degree k*64 + (63-p).
    // Highest degree = lowest bit position → use __builtin_ctzll.
    int K = a.deg / WL;
    int LL = a.deg - WL * K;

    for (int k = K; k >= 0; k--) {
        uint64_t blob = 0ULL;
        for (int j = 0; j < resolution_; j++)
            blob |= a.coord(j)[k];
        // In the top word, mask out bits for degrees above a.deg:
        // degrees 0..LL within this word → bit positions (63-LL)..63
        if (k == K)
            blob &= ~0ULL << (WL - 1 - LL);
        if (blob) {
            int ctz = __builtin_ctzll(blob);
            a.deg = k * WL + (WL - 1 - ctz);
            uint64_t mask = 1ULL << ctz;
            for (int j = 0; j < resolution_; j++) {
                if (a.coord(j)[k] & mask) {
                    a.indicemaxdeg = j;
                    return;
                }
            }
        }
        LL = WL - 1;
    }
    a.deg = INT_MIN;
    a.indicemaxdeg = INT_MIN;
}

// ══════════════════════════════════════════════════════════════════════════
// add_polvect — res = a XOR b  (EXACT logic from original code)
// ══════════════════════════════════════════════════════════════════════════

void PolyLatBase::add_polvect(PolVect& res, PolVect& a, PolVect& b) {
    if (a.deg == b.deg) {
        res.deg = a.deg;
        int nw = a.deg / WL;
        for (int k = nw; k >= 0; k--)
            for (int j = 0; j < resolution_; j++)
                res.coord(j)[k] = a.coord(j)[k] ^ b.coord(j)[k];
        update_polvect_deg(res);
        return;
    }
    if (a.deg > b.deg) {
        res.deg = a.deg;
        res.indicemaxdeg = a.indicemaxdeg;
        int bK = b.deg / WL;
        int bLL = b.deg - bK * WL;
        for (int j = 0; j < resolution_; j++) {
            uint64_t* rd = res.coord(j);
            const uint64_t* ad = a.coord(j);
            uint64_t* bd = b.coord(j);  // NOTE: b is modified (tail masking)
            for (int k = a.deg / WL; k > bK; k--)
                rd[k] = ad[k];
            // Mask out insignificant bits in b (same as original C code)
            bd[bK] &= ~0ULL << (WL - bLL - 1);
            for (int k = bK; k >= 0; k--)
                rd[k] = ad[k] ^ bd[k];
        }
        return;
    }
    // a.deg < b.deg: should not happen in algorithm's usage
}

// ══════════════════════════════════════════════════════════════════════════
// multbys — r = a * z^S  (EXACT logic from original code)
// ══════════════════════════════════════════════════════════════════════════

void PolyLatBase::multbys(PolVect& r, const PolVect& a, int S) {
    int N = (a.deg + S) / WL;
    int nbblocks = S / WL;
    int nbbits = S - nbblocks * WL;

    if (nbbits != 0) {
        int WLmn = WL - nbbits;
        for (int ii = 0; ii < resolution_; ii++) {
            uint64_t* rd = r.coord(ii);
            const uint64_t* ad = a.coord(ii);
            for (int i = 0; i < nbblocks; i++)
                rd[i] = 0ULL;
            int i = nbblocks;
            rd[i++] = ad[0] >> nbbits;
            for (; i <= N; i++)
                rd[i] = (ad[i - nbblocks] >> nbbits) | (ad[i - nbblocks - 1] << WLmn);
        }
    } else {
        for (int ii = 0; ii < resolution_; ii++) {
            uint64_t* rd = r.coord(ii);
            const uint64_t* ad = a.coord(ii);
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
// multbys2 — r XOR= a * z^S in-place  (EXACT logic from original code)
// ══════════════════════════════════════════════════════════════════════════

void PolyLatBase::multbys2(PolVect& r, const PolVect& a, int S) {
    int nbblocks = S / WL;
    int nbbits = S - nbblocks * WL;

    if (a.deg + S == r.deg) {
        int N = r.deg / WL;
        if (nbbits != 0) {
            int WLmn = WL - nbbits;
            for (int ii = 0; ii < resolution_; ii++) {
                uint64_t* rd = r.coord(ii);
                const uint64_t* ad = a.coord(ii);
                int i = nbblocks;
                rd[i++] ^= ad[0] >> nbbits;
                for (; i <= N; i++)
                    rd[i] ^= (ad[i - nbblocks] >> nbbits) | (ad[i - nbblocks - 1] << WLmn);
            }
        } else {
            for (int ii = 0; ii < resolution_; ii++) {
                uint64_t* rd = r.coord(ii);
                const uint64_t* ad = a.coord(ii);
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
                uint64_t* rd = r.coord(ii);
                const uint64_t* ad = a.coord(ii);
                int i = nbblocks;
                rd[i++] ^= ad[0] >> nbbits;
                for (; i < N; i++)
                    rd[i] ^= (ad[i - nbblocks] >> nbbits) | (ad[i - nbblocks - 1] << WLmn);
                rd[i] ^= (ad[i - nbblocks - 1] << WLmn);
            }
        } else {
            for (int ii = 0; ii < resolution_; ii++) {
                uint64_t* rd = r.coord(ii);
                const uint64_t* ad = a.coord(ii);
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
// DualBase
// ══════════════════════════════════════════════════════════════════════════

void PolyLatBase::dual_base(const std::vector<BitVect>& polys, const BitVect& M, int res) {
    for (int i = 0; i < maxresolution_; i++) {
        perm_[i] = i;
        invperm_[i] = i;
    }
    resolution_ = res;

    vb(0).zero_lazy();
    set_element_from_bv(vb(0), 0, M);

    for (int i = 1; i < res; i++) {
        vb(i).zero_lazy();
        set_element_from_bv(vb(i), 0, polys[i]);
        set_element_one(vb(i), i);
    }
}

void PolyLatBase::dual_base_increase(const std::vector<BitVect>& polys) {
    resolution_++;
    int newi = resolution_ - 1;
    for (int j = 0; j <= newi - 1; j++)
        set_element_zero(vb(j), newi);
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
    int n = m + 1;
    int nw = (n + 1 + 63) / 64;
    std::vector<std::vector<uint64_t>> A(n, std::vector<uint64_t>(nw, 0ULL));

    auto setbit = [&](int row, int col) {
        A[row][col / 64] |= 1ULL << (63 - (col % 64));
    };
    auto getbit = [&](int row, int col) -> int {
        return (A[row][col / 64] >> (63 - (col % 64))) & 1;
    };

    for (int i = 0; i <= m; i++) {
        for (int j = i; j <= m; j++) {
            if (element_norme_equal(vect_[i], j, vect_[i].deg))
                setbit(j, i);
        }
    }
    int bmp1_norme = vect_[m + 1].deg;
    for (int j = 0; j <= m; j++) {
        if (element_norme_equal(vect_[m + 1], j, bmp1_norme))
            setbit(j, m + 1);
    }

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
// Lenstra's algorithm — EXACT logic from original code
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
