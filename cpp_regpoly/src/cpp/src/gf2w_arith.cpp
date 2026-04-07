#include "gf2w_arith.h"
#include <vector>
#include <cstdlib>

uint64_t gf2w_multiply_z(uint64_t a, int k, uint64_t modM) {
    for (int i = 0; i < k; i++) {
        if (a & 1ULL)
            a = (a >> 1) ^ modM;
        else
            a >>= 1;
    }
    return a;
}

uint64_t gf2w_multiply_poly(uint64_t a, uint64_t b, int w, uint64_t modM) {
    uint64_t res = 0, verif = 1;
    for (int i = 0; i < w; i++) {
        if (b & verif)
            res ^= gf2w_multiply_z(a, w - i - 1, modM);
        verif <<= 1;
    }
    return res;
}

uint64_t gf2w_multiply_normal(uint64_t a, uint64_t b, int w,
                               const uint64_t* table) {
    if (a == 0 || b == 0) return 0;
    uint64_t res = 0;
    uint64_t verifa = 1ULL << (w - 1);
    for (int j = 0; j < w; j++) {
        uint64_t verifb = 1ULL << (w - 1);
        for (int i = 0; i < w; i++) {
            if ((verifa & a) && (verifb & b))
                res ^= table[i * w + j];
            verifb >>= 1;
        }
        verifa >>= 1;
    }
    return res;
}

std::vector<uint64_t> gf2w_make_table(int w, uint64_t modM) {
    // Step 1: Compute x[j] = x^j in polynomial basis
    // x[0] = 1 << (w-1)  (x^0 = 1, stored MSB-first within w bits)
    std::vector<uint64_t> x(2 * w, 0);
    x[0] = 1ULL << (w - 1);
    for (int j = 0; j < 2 * w - 1; j++) {
        if (x[j] & 1ULL)
            x[j + 1] = (x[j] >> 1) ^ modM;
        else
            x[j + 1] = x[j] >> 1;
    }

    // Step 2: Compute Frobenius elements xx[j] = x^(2^j)
    // xx[0] = x[1] = x^1
    std::vector<uint64_t> xx(w, 0);
    xx[0] = x[1];
    for (int j = 0; j < w - 1; j++) {
        uint64_t verif = 1ULL << (w - 1);
        xx[j + 1] = 0;
        for (int i = 0; i < w; i++) {
            if (xx[j] & verif)
                xx[j + 1] ^= x[2 * i];
            verif >>= 1;
        }
    }

    // Step 3: Build w x w matrix C: column j = polynomial repr of xx[j]
    // C is stored as w rows, each row is a w-bit value
    std::vector<uint64_t> C(w, 0);
    for (int j = 0; j < w; j++) {
        uint64_t verif = 1ULL << (w - 1);
        for (int i = 0; i < w; i++) {
            if (xx[j] & verif)
                C[i] |= 1ULL << (w - 1 - j);
            verif >>= 1;
        }
    }

    // Step 4: Invert C via Gauss-Jordan to get InvC
    std::vector<uint64_t> InvC(w, 0);
    for (int i = 0; i < w; i++)
        InvC[i] = 1ULL << (w - 1 - i);  // identity matrix

    std::vector<uint64_t> work(C);
    for (int col = 0; col < w; col++) {
        uint64_t bit = 1ULL << (w - 1 - col);
        // Find pivot
        int pivot = -1;
        for (int i = col; i < w; i++) {
            if (work[i] & bit) {
                pivot = i;
                break;
            }
        }
        if (pivot < 0) {
            // No normal basis — this should not happen in normal use
            // Return empty table to signal error
            return std::vector<uint64_t>();
        }
        if (pivot != col) {
            std::swap(work[col], work[pivot]);
            std::swap(InvC[col], InvC[pivot]);
        }
        for (int i = 0; i < w; i++) {
            if (i != col && (work[i] & bit)) {
                work[i] ^= work[col];
                InvC[i] ^= InvC[col];
            }
        }
    }

    // Step 5: Transpose InvC -> TransInvC
    std::vector<uint64_t> TransInvC(w, 0);
    for (int i = 0; i < w; i++) {
        for (int j = 0; j < w; j++) {
            if (InvC[i] & (1ULL << (w - 1 - j)))
                TransInvC[j] |= 1ULL << (w - 1 - i);
        }
    }

    // Step 6: Build multiplication table
    // table[i*w+j] = normal-basis product of xx[i] and xx[j]
    std::vector<uint64_t> table(w * w, 0);
    for (int i = 0; i < w; i++) {
        for (int j = 0; j < w; j++) {
            uint64_t temp = gf2w_multiply_poly(xx[i], xx[j], w, modM);
            uint64_t entry = 0;
            for (int jj = 0; jj < w; jj++) {
                if ((temp >> (w - 1 - jj)) & 1ULL)
                    entry ^= TransInvC[jj];
            }
            table[i * w + j] = entry;
        }
    }

    return table;
}
