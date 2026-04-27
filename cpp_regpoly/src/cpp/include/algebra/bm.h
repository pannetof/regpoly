#pragma once
#include "bitvect.h"
#include "generator.h"

// Packed Berlekamp-Massey over F_2.
//
// Runs gen (after init(init_state)) for 2*K steps, observing the LSB of
// get_output() at each step.  Returns the linear complexity L of the
// resulting binary sequence.  If out_min_poly is non-null, writes the
// recovered minimal polynomial into a freshly-sized BitVect of K bits,
// MSB-first: bit j holds the coefficient of z^(K-j).  (Same convention
// as Generator::char_poly() returns.)
//
// K should be an upper bound on the actual minimal polynomial degree
// (typically gen.k() — but for non-full-period generators the true
// minimal polynomial may be smaller).
int packed_bm(const Generator& gen,
              const BitVect& init_state,
              int K,
              BitVect* out_min_poly,
              int bit_idx = 0);   // which output bit to BM (default 0)
