#include "bm.h"
#include <vector>
#include <cstring>

// Polynomials C, B and the output sequence are stored as packed uint64_t
// arrays (MSB-first, matching BitVect convention).  This gives ~64×
// speedup over the element-wise int-array version for the inner loop.

namespace {

inline int pk_get(const uint64_t* data, int i) {
    return (data[i >> 6] >> (63 - (i & 63))) & 1;
}

inline void pk_set1(uint64_t* data, int i) {
    data[i >> 6] |= 1ULL << (63 - (i & 63));
}

// Extract 64 contiguous bits starting at bit position pos.
// Caller must ensure data[] has ≥ 1 word of padding at the end.
inline uint64_t pk_extract64(const uint64_t* data, int pos) {
    int w = pos >> 6;
    int b = pos & 63;
    if (b == 0) return data[w];
    return (data[w] << b) | (data[w + 1] >> (64 - b));
}

// C ^= (B right-shifted by 'shift' bits).
// Polynomial semantics: C(z) += z^shift · B(z).
void pk_shift_xor(uint64_t* C, const uint64_t* B, int shift, int nw) {
    int ws = shift >> 6;
    int bs = shift & 63;
    if (ws >= nw) return;
    if (bs == 0) {
        for (int i = nw - 1; i >= ws; i--)
            C[i] ^= B[i - ws];
    } else {
        int comp = 64 - bs;
        for (int i = nw - 1; i > ws; i--)
            C[i] ^= (B[i - ws] >> bs) | (B[i - ws - 1] << comp);
        C[ws] ^= B[0] >> bs;
    }
}

}  // namespace

int packed_bm(const Generateur& gen,
              const BitVect& init_state,
              int K,
              BitVect* out_min_poly,
              int bit_idx)
{
    auto g = gen.copy();
    g->init(init_state);

    int N_total = 2 * K;
    int nw_seq  = (N_total + 63) / 64 + 2;  // +2 padding for pk_extract64
    int nw_poly = (K + 1 + 63) / 64 + 1;    // +1 padding

    // Reversed sequence: rseq bit i stores s[N_total-1-i].
    // This makes the discrepancy window contiguous for word-level AND+popcount.
    std::vector<uint64_t> rseq(nw_seq, 0);

    std::vector<uint64_t> C(nw_poly, 0);
    std::vector<uint64_t> B(nw_poly, 0);
    std::vector<uint64_t> T(nw_poly);  // pre-allocated temp

    C[0] = 1ULL << 63;  // C[0] = 1 (bit 0 = MSB)
    B[0] = 1ULL << 63;  // B[0] = 1

    int Len = 0;
    int x = 1;

    for (int N = 0; N < N_total; N++) {
        g->next();
        int sN = g->get_output().get_bit(bit_idx);

        if (sN) pk_set1(rseq.data(), N_total - 1 - N);

        // Discrepancy: d = Σ_{j=0}^{Len} C[j] · s[N-j]
        //              = Σ_{j=0}^{Len} C[j] · rseq[base + j]
        int base = N_total - 1 - N;
        int d = 0;
        int nw_active = (Len + 64) / 64;
        for (int w = 0; w < nw_active; w++) {
            uint64_t cw = C[w];
            if (cw)
                d ^= __builtin_popcountll(cw & pk_extract64(rseq.data(), base + 64 * w));
        }
        d &= 1;

        if (d == 0) {
            x++;
        } else if (2 * Len > N) {
            pk_shift_xor(C.data(), B.data(), x, nw_poly);
            x++;
        } else {
            std::memcpy(T.data(), C.data(), nw_poly * sizeof(uint64_t));
            pk_shift_xor(C.data(), B.data(), x, nw_poly);
            Len = N + 1 - Len;
            std::memcpy(B.data(), T.data(), nw_poly * sizeof(uint64_t));
            x = 1;
        }
    }

    if (out_min_poly) {
        *out_min_poly = BitVect(K);
        for (int j = 0; j < K; j++)
            if (pk_get(C.data(), K - j))
                out_min_poly->set_bit(j, 1);
    }
    return Len;
}
