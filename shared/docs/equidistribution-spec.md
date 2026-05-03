# Matricial Computation of Equidistribution for Black-Box F_2-Linear Generators

**Scope.** This document specifies an algorithm and implementation plan for computing
the dimension-of-equidistribution profile `k(v)`, `v = 1, …, w`, of a pseudorandom
number generator that is known to be `F_2`-linear but whose internal structure
(characteristic polynomial, invariant subspaces, factorization) is not assumed to
be known in advance. The generator may be non-full-period: no assumption is made
that the characteristic polynomial is primitive.

**Inputs required from the generator.**

1. A state of dimension `n` over `F_2` (so the full state is `n` bits; `n = w·N`
   in the usual LFSR-array formulation).
2. Word size `w` of the output (32, 64, or 128).
3. A `step()` routine advancing the internal state by one step.
4. A routine to read the current `w`-bit output word.
5. The ability to set the state to an arbitrary `n`-bit value.

No knowledge of `χ_f`, no access to `f` as a matrix, no tempering assumptions.

**Output.** The array `k(v)` for `v = 1, …, w`, together with the recovered factor
`φ(t)` of the characteristic polynomial on which equidistribution is certified,
the dimension `p = deg φ` of the invariant subspace `V`, and the certified defect
ratio `2^(-p)`.

---

## 1. Algorithm Overview

The pipeline is seven stages:

1. Recover `χ_f(t)` via Berlekamp–Massey on scalar output sub-streams.
2. Factor `χ_f` over `F_2`.
3. Select the largest-period irreducible factor `φ`; set `p = deg φ`,
   `V = Ker φ(f)`.
4. Construct a basis `B` of `V` as an `n × p` bit-matrix.
5. Run the matricial DE core on `V` to extract all `k(v)` in one sweep.
6. Initialize-state guard: ensure the user-facing initializer produces states
   with nonzero `V`-component.
7. Verify via three sanity checks and a cross-check against a known generator.

Each stage is specified below.

---

## 2. Stage 1 — Recover `χ_f` via Berlekamp–Massey

### 2.1 Procedure

Pick a nonzero initial state `s_0`. Define a scalar `F_2`-valued output functional
`ℓ : state → F_2` — the simplest viable choice is "bit 0 of the first output word
after one step." Produce the sequence

    b_j = ℓ(f^j · s_0),   j = 0, 1, …, 2n + 32.

Run Berlekamp–Massey on `(b_j)`. It returns a polynomial `μ_{s_0, ℓ}(t)` which
divides `χ_f(t)`.

Repeat with 5 independent random initial states and/or 2–3 different scalar
functionals `ℓ` (e.g. bit 0 of word 0, bit 7 of word 0, bit 0 of word 1). Compute

    χ_f(t)  :=  lcm of all recovered minimal polynomials.

### 2.2 Verification

Before proceeding, require `deg χ_f = n`. If not, repeat with more sampling
functionals. Failure to reach degree `n` after ~20 trials indicates either that
the generator has a minimal polynomial genuinely smaller than `n` (i.e. the
declared state dimension is larger than the reachable state space), or a bug in
the `step()` wiring.

### 2.3 Cost

`O(n²)` per Berlekamp–Massey run, `O(n)` generator steps per sequence. With
`n ≈ 20000`, this stage runs in well under a minute.

### 2.4 Implementation notes

- Bit-packed state, bit-packed sequences. BM over `F_2` is dominated by XORs.
- Do *not* use NTL's generic `MinPolySeq` over a prime field; use a dedicated
  `F_2` BM implementation for speed. A from-scratch C implementation is ~100
  lines.

---

## 3. Stage 2 — Factor `χ_f` over `F_2`

Use distinct-degree factorization followed by equal-degree (Cantor–Zassenhaus).
Output the full factorization

    χ_f(t) = ∏_j φ_j(t)^{e_j}

with `φ_j` irreducible, `deg φ_j = d_j`, multiplicity `e_j`.

### 3.1 Implementation

Delegate to NTL (`CanZass` on `GF2X`) or FLINT (`nmod_poly_factor` with
modulus 2, or the dedicated `F2` routines). Do not reimplement.

### 3.2 Expected runtime

For `n ≈ 20000`, seconds to low minutes. For `n ≈ 100000`, still tractable but
benefits from NTL's tuned implementations.

---

## 4. Stage 3 — Select the Invariant Subspace

### 4.1 Selection rule

Among the irreducible factors `{φ_j}` of `χ_f`, select the one of largest degree
whose period `ord(φ_j)` in `F_2[t] / φ_j` can be certified maximal (i.e. equal to
`2^{d_j} − 1`).

- **Primary path:** if the largest-degree irreducible factor has degree `p` equal
  to a Mersenne exponent, `2^p − 1` is prime and `φ` is automatically primitive
  provided `φ(0) ≠ 0`. Verify `φ(0) = 1` (for nonzero constant term) — done.
- **Secondary path:** if `p` is not a Mersenne exponent, attempt to obtain the
  prime factorization of `2^p − 1` via local tables (Cunningham, Wagstaff) and
  FactorDB fallback. With the factorization in hand, check primitivity by the
  standard Knuth TAOCP Vol. 2 §3.2.2 algorithm.
- **Fallback:** if `2^p − 1` cannot be factored, either
  - use the largest irreducible factor without certifying primitivity — the
    period on `V` divides `2^p − 1` and equidistribution bounds remain valid
    with defect ratio `2^{−p}`; the claimed period is a lower bound, or
  - use the largest factor for which primitivity *can* be certified, even at
    smaller `p`.

Set

    φ  := the selected factor,
    p  := deg φ,
    V  := Ker φ(f)  ⊂  F_2^n,    dim V = p.

### 4.2 Note on multiplicity

If `φ` appears in `χ_f` with multiplicity `e > 1`, `Ker φ(f)` still has dimension
`p`, not `e·p` — it is the simple kernel, not the generalized kernel. This is the
subspace we want.

---

## 5. Stage 4 — Construct a Basis `B` of `V`

Two routes. Use **Route B** by default; fall back to Route A if Route B's rank
check fails.

### 5.1 Route B (orbit-based — preferred)

Let `ψ(t) := χ_f(t) / φ(t)`, a polynomial of degree `n − p`.

1. Pick a random state `s ∈ F_2^n`.
2. Compute `s' := ψ(f) · s` by evaluating `ψ` at `f` via Horner's scheme:
   iterate `step()` from `s` and accumulate the `ψ`-coefficients. Cost:
   `n − p` generator steps, `n − p` state-XORs.
3. Verify `s' ≠ 0` and `φ(f) · s' = 0`. If `s' = 0`, retry with new `s`
   (this happens with probability `2^{−p}`, so effectively never).
4. Collect the orbit `{f^k · s' : k = 0, 1, …, p − 1}` as columns of an
   `n × p` matrix `B`.
5. Check `rank_{F_2}(B) = p`. If yes, done: the columns of `B` span `V`.
6. If rank `< p` (generically impossible but possible for pathological
   generators where the minimal polynomial of `s'` is a proper divisor of `φ`),
   repeat from step 1.

### 5.2 Route A (annihilator-based — fallback)

1. Pick random `s`, compute `s' := ψ(f) · s` as above. Now `s' ∈ V` by
   construction.
2. Repeat until `p` linearly independent such `s'` have been collected.
3. Gauss-eliminate to confirm dim = `p` and to obtain a clean basis.

Route A is more robust but costs `~p` applications of `ψ(f)` instead of `~1`.

### 5.3 Representation of `B`

Store `B` as a bit-packed `n × p` matrix, rows as `p/W`-word bitstrings
(`W = 64`). Memory: `n · p / 8` bytes. For `n = p = 20000`: 50 MB. Fine.

---

## 6. Stage 5 — The Matricial DE Core

### 6.1 Object maintained

For each output phase `i` (see §6.4), maintain a single echelon form `E_i` over
`F_2` with `p` columns. Rows are successively inserted; at most `p` pivots exist.

Also maintain, for each accuracy `v ∈ {1, …, w}`:
- `c_i(v)`: current rank of the `v`-restricted submatrix for phase `i`;
- `alive_i(v)`: boolean flag, true while `c_i(v) = k·v` has kept up.

### 6.2 Key structural observation

Order the rows of the output matrix bit-by-bit within each output word, MSB
first:

    output 0 bit 1, output 0 bit 2, …, output 0 bit w,
    output 1 bit 1, output 1 bit 2, …, output 1 bit w,
    …

Then the matrix `M_k(v)` for accuracy `v` is obtained from `M_k(w)` by keeping,
within each block of `w` consecutive rows, only the first `v`. **A single echelon
form on the full-precision rows, maintained in this order, gives every `k(v)`
simultaneously.**

### 6.3 The inner loop (single phase)

```
initialize E empty; c(v) ← 0, alive(v) ← true for v = 1..w
for k = 1, 2, 3, …:
    compute this step's p dual vectors for bits 1..w of the output word,
        represented as p-bit rows over the basis B of V
    for j = 1, …, w:
        r ← dual vector of bit j
        attempt F_2-echelon insertion of r into E
        if r was linearly independent:
            for every v ≥ j with alive(v):
                c(v) ← c(v) + 1
        for every v ≥ j with alive(v):
            if c(v) < k·v:            // only possible if r was dependent
                alive(v) ← false
                k(v) ← k − 1
    if no v is still alive: break
    if k·w ≥ p: break                   // upper bound k(v) ≤ ⌊p/v⌋ saturated
```

### 6.4 Phase handling (for `w_output < w_state`)

If the generator emits `w`-bit words but is designed around a larger internal
word (e.g. SFMT emits 32-bit words from a 128-bit internal word, giving four
phases `i ∈ {0, 1, 2, 3}` per internal step), maintain **four independent
echelon forms** `E_0, …, E_3`, one per phase, sharing the same `V`-basis `B`
and the same advance of the generator state. Each `E_i` sees only the output
bits of phase `i`.

The reported `k(v)` for the generator is

    k(v) = min_i k_i(v).

For generators with no phase mismatch, use a single phase (`i = 0`).

### 6.5 Lazy output generation (the trick that makes this affordable)

Never materialize `f` as an `n × n` matrix. Instead:

1. Pack the `p` basis vectors of `V` (the columns of `B`) as `p` "virtual
   registers," each an `n`-bit state of the generator.
2. At each `k`-step, advance all `p` virtual registers by one call to `step()`
   each. This is `p` independent generator steps — parallelizable, and each
   is as cheap as a single generator step.
3. Read the `w`-bit output from each virtual register. Transpose the resulting
   `p × w` slab to obtain `w` dual vectors, each of length `p` bits — exactly
   the rows we need to insert.

Cost per `k`-step: `p` generator steps. For SFMT19937 with `p ≈ 20000`: a few
hundred thousand XOR/shift operations, i.e. ~millisecond.

### 6.6 Bit-packed echelon reduction

Each row of `E` is `p` bits, stored as `⌈p/64⌉` `uint64_t` words. Insertion:

- Find the highest set bit of the incoming row `r`.
- If no pivot exists at that position: `r` is independent, store as new pivot.
- If a pivot exists: XOR the pivot into `r`, repeat.

Worst case: `O(p² / 64)` word-ops per insertion. Total insertions: at most `p`
(rank is capped). Total echelon work per phase: **`O(p³ / 64)`** word-ops.

### 6.7 Total complexity per phase

    Echelon maintenance:    O(p³ / W)     ← dominant
    Generator advance:      O(p · n / W) per k-step × O(p/w) k-steps
                          = O(p² · n / (w · W))

For `p = n = 19937`, `w = 32`, `W = 64`: echelon is `~10^{11}` word-ops, i.e.
seconds to a minute. Generator advance is `~10^{10}`, i.e. comparable but
slightly faster.

Four phases (SFMT case) multiply by 4. Total: minutes on a single core.

### 6.8 Optional refinement: weighted-norm prioritization

The Couture–L'Ecuyer 1993 weighted-norm trick (and the Harase refinements on
top of it) reorder *which* new dual vectors to insert, prioritizing by a weight
reflecting how close each `v` is to its bound `⌊p/v⌋`. Savings are a 3–10×
constant factor. **Defer to a second pass** — implement §6.3 first, validate,
then optimize if the base version is too slow for a parameter-search loop.

---

## 7. Stage 6 — Initial-State Guard

To use the generator correctly, the public `seed()` routine must produce states
with nonzero `V`-component. Analogue of SFMT's `6d736d6d` MSB trick.

### 7.1 Procedure

1. After seeding, compute the `V`-projection: `s_V := (a specific projector)(s)`
   where the projector is expressible in terms of `φ` and `ψ = χ_f / φ`. A clean
   choice: `s_V := ψ(f) · s / normalizer`, but for guarding purposes we only
   need `s_V ≠ 0`, which holds iff `ψ(f) · s ≠ 0`.
2. If `ψ(f) · s = 0`, perturb `s` (e.g. flip one bit in a predetermined
   location) and retry. This happens with probability `2^{−p}`.

For generators where the state has a natural "magic value" position (cf. SFMT),
hardcoding that value into the initializer is sufficient and cheaper than
runtime projection.

---

## 8. Stage 7 — Verification

Three mandatory sanity checks, in order. A failure in any check means `k(v)` is
not to be trusted.

### 8.1 Period check

From a random state with `s_V ≠ 0`, iterate the generator and confirm the
orbit on `V` does not close before `ord(φ)` steps. For Mersenne-exponent `p`,
confirm `f^{(2^p − 1)/q} · s_V ≠ s_V` for every small prime `q` dividing `2^p − 1`
(standard primitivity test on the orbit, Knuth Vol. 2 §3.2.2).

### 8.2 Trivial DE bounds

Confirm:

- `k(w) ≥ 1` (a single `w`-bit output should be `w`-uniform; failure here means
  output tempering is catastrophically broken or `B` is wrong).
- `k(1) ≤ p` (upper bound always holds; violation means the echelon
  bookkeeping is buggy).
- `k(v) ≤ ⌊p/v⌋` for all `v` (same: upper bound must hold).

### 8.3 Cross-check against a known generator

Run the *entire* pipeline — stages 1 through 5 — on MT19937, for which the
`k(v)` profile is published (cf. Matsumoto–Nishimura 1998, Table 2; also SFMT
paper Table 3). Require exact agreement on all 32 values of `k(v)`. This is
the single most important test; it catches basis-construction bugs,
phase-handling bugs, and off-by-one errors in the echelon loop.

A secondary cross-check against SFMT19937 (Table 3 of Saito–Matsumoto 2008) is
recommended once the four-phase handling is implemented.

---

## 9. Data Structures and File Layout (for `regpoly`)

Proposed module layout within the existing `MinimalCode` structure:

```
src/
  de/
    bm.c, bm.h              ← Berlekamp–Massey over F_2
    factor.c, factor.h      ← wrapper around NTL/FLINT factorization
    basis.c, basis.h        ← Route A / Route B basis construction
    echelon.c, echelon.h    ← bit-packed F_2 echelon maintenance
    de_core.c, de_core.h    ← §6 main loop, phase handling
    guard.c, guard.h        ← initial-state guard
    verify.c, verify.h      ← §8 sanity checks
  de_driver.c               ← CLI: takes a generator, emits k(v) table
tests/
  de/
    mt19937_kv.c            ← cross-check against MT19937 published k(v)
    sfmt19937_kv.c          ← cross-check against SFMT19937 Table 3
    small_known.c           ← degenerate small-n generators with
                                hand-computed k(v)
```

**Key types.**

```c
typedef struct {
    uint64_t *words;   // ⌈p/64⌉ words
    size_t    p;
} bitrow_t;

typedef struct {
    bitrow_t *pivots;  // indexed by leading-bit position
    size_t    p;
    size_t    rank;
} echelon_t;

typedef struct {
    size_t   n, w, p;
    bitrow_t *B_cols; // p columns of B, each n bits — stored col-major
    // ... generator-specific callbacks:
    void (*step)(void *state);
    void (*read)(const void *state, uint8_t *out_w_bytes);
    void (*set_state)(void *state, const uint8_t *in_n_bytes);
    size_t state_size;
} de_context_t;
```

---

## 10. Limits of This Specification

- **Nonlinear generators** are out of scope. The entire framework rests on
  `F_2`-linearity of `f` and of the output functionals.
- **Tempered generators with nonlinear tempering** (e.g. XOR-shift-multiply on
  output) break the linearity of the output functionals. This pipeline computes
  `k(v)` for the *pre-tempering* sequence only; post-tempering DE requires
  separate analysis.
- **Generators with `n > 10^5`** are in principle handled but factorization
  (stage 2) and `ψ(f)` evaluation (stage 4) become the bottleneck. Consider
  whether the generator's design already exposes `φ` — if so, skip stages 1–3.

---

## 11. References

1. M. Saito, M. Matsumoto. *SIMD-oriented Fast Mersenne Twister: a 128-bit
   Pseudorandom Number Generator.* In *Monte Carlo and Quasi-Monte Carlo
   Methods 2006*, Springer, 2008.
2. R. Couture, P. L'Ecuyer, S. Tezuka. *On the distribution of k-dimensional
   vectors for simple and combined Tausworthe sequences.* Math. Comp. 60 (1993),
   749–761.
3. M. Matsumoto, T. Nishimura. *Mersenne Twister: A 623-dimensionally
   equidistributed uniform pseudorandom number generator.* ACM TOMACS 8 (1998),
   3–30.
4. S. Harase. *On the F_2-linear relations of Mersenne Twister pseudorandom
   number generators.* Math. Comput. Simulation 100 (2014), 103–113, and
   subsequent refinements.
5. D. E. Knuth. *The Art of Computer Programming, Vol. 2: Seminumerical
   Algorithms.* 3rd ed., Addison-Wesley, 1997. §3.2.2.
6. F. Panneton. (author's own paper — fill in citation locally: "Random Number
   Generators Based on Linear Recurrences in `GF(2^w)`," with P. L'Ecuyer.)

---

*Document version: draft 1. Amend in-repo as implementation reveals gaps.*
