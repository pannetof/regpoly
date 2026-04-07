#include "generateur.h"
#include <vector>
#include <cstring>

// ── Default get_output ──────────────────────────────────────────────────

BitVect Generateur::get_output() const {
    BitVect out(L_);
    out.copy_part_from(state_, L_);
    return out;
}

// ── Default get_transition_state ───────────────────────────────��─────────

void Generateur::get_transition_state(uint64_t* out_words, int out_nwords) const {
    // Copy first k_ bits of state_ into out_words.
    // out_words is MSB-first: out_words[0] holds bits 0..63, etc.
    // state_ already stores bits MSB-first, so we can copy directly.
    BitVect tmp(k_);
    tmp.copy_part_from(state_, k_);
    int n = std::min(out_nwords, tmp.nwords());
    for (int i = 0; i < n; i++)
        out_words[i] = tmp.data()[i];
    for (int i = n; i < out_nwords; i++)
        out_words[i] = 0;
}

// ── Berlekamp-Massey (char_poly) ─────────────────────────────────────────

BitVect Generateur::char_poly() const {
    int K = k_;

    // Create a copy of self
    auto gen = copy();

    // Init with canonical vector (bit 0 set)
    BitVect init_bv(K);
    init_bv.set_bit(0, 1);
    gen->init(init_bv);

    // Collect 2*K output bits: MSB of output each time
    std::vector<int> seq(2 * K, 0);
    for (int i = 0; i < 2 * K; i++) {
        gen->next();
        // Get MSB of output (bit 0) — uses get_output() to handle
        // generators with circular buffers or output rotation
        seq[i] = gen->get_output().get_bit(0);
    }

    // Berlekamp-Massey algorithm (SequenceMinimalPolynomial)
    std::vector<int> C(K + 1, 0);
    std::vector<int> B(K + 1, 0);
    C[0] = B[0] = 1;
    int Len = 0;
    int x = 1;

    for (int N = 0; N < 2 * K; N++) {
        int d = seq[N];
        for (int j = 1; j <= Len; j++) {
            if (C[j])
                d ^= seq[N - j];
        }
        d &= 1;

        if (d == 0) {
            x += 1;
        } else if (2 * Len > N) {
            // C = C - x^x * B (in GF(2): XOR)
            std::vector<int> temp(K + 1, 0);
            for (int i = 0; i + x <= K; i++)
                temp[i + x] = B[i];
            for (int i = 0; i <= K; i++)
                C[i] ^= temp[i];
            x += 1;
        } else {
            std::vector<int> T = C;
            std::vector<int> temp(K + 1, 0);
            for (int i = 0; i + x <= K; i++)
                temp[i + x] = B[i];
            for (int i = 0; i <= K; i++)
                C[i] ^= temp[i];
            Len = N + 1 - Len;
            B = T;
            x = 1;
        }
    }

    // Pack into BitVect: public bit j holds C[K - j]
    BitVect bv(K);
    for (int j = 0; j < K; j++)
        bv.set_bit(j, C[K - j]);
    return bv;
}

// ── Transition matrix ────────────────────────────────────────────────────

std::vector<BitVect> Generateur::transition_matrix() const {
    int K = k_;
    auto gen = copy();

    BitVect bc(K);
    bc.set_bit(0, 1);

    int nw = (K + 63) / 64;

    // Phase 1: build columns (= rows of Aᵀ). Column i is the k-bit
    // state after one step with basis vector e_i.
    std::vector<std::vector<uint64_t>> cols(K, std::vector<uint64_t>(nw, 0));
    for (int i = 0; i < K; i++) {
        gen->init(bc);
        bc.rshift(1);
        gen->next();
        gen->get_transition_state(cols[i].data(), nw);
    }

    // Phase 2: transpose into rows. Row[row] bit i is set iff
    // cols[i] has bit row set.
    std::vector<BitVect> rows(K, BitVect(K));
    for (int i = 0; i < K; i++) {
        for (int row = 0; row < K; row++) {
            // Check bit 'row' in cols[i]
            int w = row / 64;
            int b = 63 - (row % 64);
            if ((cols[i][w] >> b) & 1)
                rows[row].set_bit(i, 1);
        }
    }
    return rows;
}
