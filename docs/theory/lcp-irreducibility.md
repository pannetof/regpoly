# LCP-Based Fast Rejection Filter and Primitivity Testing for $\mathrm{GF}(2)$-Linear Generators

## Implementation Specification

---

## 1. Purpose and Scope

This document describes two complementary tools for testing whether the characteristic
polynomial $P(x)$ of a $\mathrm{GF}(2)$-linear generator is **primitive**:

1. **LCP fast rejection filter (BM-based):** cheap $O(k^2)$ test that can *reject* bad
   parameters with certainty, but *cannot* certify irreducibility or primitivity on its own.
2. **Exact irreducibility test (Rabin / `IterIrredTest`):** $O(k^{2.585})$ test that certifies
   irreducibility exactly. For Mersenne prime exponent $k$, this simultaneously certifies
   primitivity.

The intended workflow is to use the LCP filter to cheaply discard obviously bad parameters,
then apply the exact test only to candidates that survive.

> **⚠ Critical limitation of the LCP test**
>
> $L = k$ observed by BM does **not** imply $P(x)$ is irreducible — not even when $k$ is a
> Mersenne prime exponent. A reducible $P(x)$ with large balanced factors (e.g., degrees
> 100 + 507 when $k = 607$) will produce $L = k$ for almost every random seed, and the LCP
> test will not detect the reducibility. The LCP test is a rejection filter only. The exact
> irreducibility test (`IterIrredTest` in NTL, or Rabin's algorithm) is required to certify
> irreducibility.

---

## 2. Theoretical Foundation

### 2.1 $\mathrm{GF}(2)$-linear generators and minimal polynomials

A $\mathrm{GF}(2)$-linear generator evolves its $k$-bit state by a linear map $A$ over $\mathrm{GF}(2)^k$:

$$\mathbf{x}_{n+1} = A \cdot \mathbf{x}_n$$

Every output bit is a fixed linear function of the state: $s_n = c^T \mathbf{x}_n$. The sequence
$\{s_n\}$ satisfies a linear recurrence over $\mathrm{GF}(2)$ whose **minimal polynomial** $m(x)$ divides
the **characteristic polynomial** $P(x) = \det(xI - A)$:

$$m(x) \mid P(x) \quad \text{always}$$
$$\deg(m(x)) \le \deg(P(x)) = k \quad \text{always}$$

### 2.2 What Berlekamp-Massey computes

Given $N \ge 2k$ bits of $\{s_n\}$, BM finds:

- $m(x)$ — the minimal polynomial of the sequence for the given seed $\mathbf{x}_0$
- $L$ — the linear complexity $= \deg(m(x))$

The relationship between $L$ and the structure of $P(x)$ depends on the seed:

$$\begin{aligned}
L = k &\iff m(x) = P(x) \text{ for this particular seed } \mathbf{x}_0 \\
L < k &\iff \mathbf{x}_0 \text{ lies in a proper } A\text{-invariant subspace (eigenspace of a factor of } P(x)\text{)}
\end{aligned}$$

### 2.3 The correct interpretation of $L = k$

**$L = k$ does NOT imply $P(x)$ is irreducible.**

If $P(x) = f(x) \cdot g(x)$ with $\deg(f) = d$ and $\deg(g) = k - d$, then:
- Seeds in $\ker(g(A))$ give $L = d < k$ — detects reducibility
- Seeds in $\ker(f(A))$ give $L = k - d < k$ — detects reducibility
- Generic seeds give $L = k$ — does NOT detect reducibility

The fraction of seeds that detect reducibility is only $\approx 2^{-d_\text{min}}$ where $d_\text{min}$ is the
degree of the smallest factor. For large balanced factors ($d \approx k/2$), this fraction is
$\approx 2^{-k/2} \approx 0$, so the LCP test essentially never detects the reducibility.

**Empirical confirmation:** the $a = 0$ MELG607 generator has $k = 607$ (a Mersenne prime
exponent), yet its characteristic polynomial is reducible. NTL's `IterIrredTest` returns
false. Its factors have degrees like 3+604 or 100+507 — neither factor degree divides 607.
BM on a random seed gives $L = 607$, failing to detect the reducibility.

### 2.4 What the LCP test IS good for

The failure direction is exact and useful:

$$L < k \text{ for any seed} \implies P(x) \text{ is definitely reducible (certain, no caveats)}$$

This is cheap ($O(k^2)$ per seed) and allows fast rejection of bad parameters before running
the expensive exact test. In a parameter search where most candidates are bad, this filter
eliminates the majority of candidates at low cost.

### 2.5 The Mersenne prime exponent property (correctly stated)

When $k$ is a Mersenne prime exponent, $2^k - 1$ is a Mersenne prime. This gives:

**Correct claim:** if $P(x)$ is **irreducible** of degree $k$ and $k$ is a Mersenne prime
exponent, then $P(x)$ is **primitive**.

**Proof:** any root $\alpha$ of an irreducible degree-$k$ polynomial has order dividing $2^k - 1$.
Since $2^k - 1$ is prime, the order is either 1 (impossible for $\alpha$ nonzero) or $2^k - 1$
(maximal). Therefore $P(x)$ is primitive. ∎

**What this does NOT say:**
- It does not say that every degree-$k$ polynomial over $\mathrm{GF}(2)$ is irreducible when $k$ is a
  Mersenne prime exponent. Reducible degree-$k$ polynomials exist for all $k > 1$.
- It does not help the LCP test detect reducibility any better.

**What this DOES say:** once irreducibility is certified by `IterIrredTest`, primitivity
follows automatically for Mersenne prime exponent degrees — no separate order check is
needed.

```
For Mersenne prime exponent k:

  IterIrredTest(P(x)) == true   ⟹   P(x) is primitive

  No order check required.
```

This eliminates the need for the Brillhart-Lehmer-Selfridge-Tuckerman-Wagstaff (BLSTW)
factor tables and the associated order verification, which is the main computational
bottleneck for large $k$ in the general case.

---

## 3. The Two-Stage Workflow

```
For each candidate parameter set:

  Stage 1 — LCP fast rejection filter (cheap: O(k²) per seed, t seeds)
  ├─ Run BM on t random seeds
  ├─ If any seed gives L < k:
  │    → REJECT immediately (certain — P(x) is reducible)
  └─ If all seeds give L = k:
       → PASS to Stage 2 (P(x) might be irreducible — not yet confirmed)

  Stage 2 — Exact irreducibility test (expensive: O(k^2.585))
  ├─ Run IterIrredTest (NTL) or Rabin's algorithm
  ├─ If irreducible:
  │    → If k is Mersenne prime exponent: ACCEPT as PRIMITIVE
  │    → Otherwise: run order check to confirm primitivity
  └─ If not irreducible:
       → REJECT (P(x) is reducible — LCP filter missed it)
```

In a typical parameter search, Stage 1 rejects most candidates cheaply (the fraction of
primitive polynomials among all degree-$k$ polynomials is $\varphi(2^k - 1) / (k \cdot 2^k) \approx 1/k$ for
large $k$, so roughly $(1 - 1/k)$ of candidates are bad). Stage 2 is only run on the small
fraction that survive Stage 1.

---

## 4. Inputs

```
generator          — the concrete generator instance (e.g. MELG607)
k                  — degree of characteristic polynomial (e.g. 607)
t                  — number of seeds for LCP filter (recommended: 20-40)
N                  — number of output bits per seed for BM (use N = 2*k + 64)
```

---

## 5. Stage 1: LCP Fast Rejection Filter

### 5.1 Generate the output bit stream

For each seed trial:

```
initialise_generator_with_seed(seed)
bit_buffer = empty
for i = 0 to N-1:
    if bit_buffer is empty:
        word = generator_next_word()
        load bit_buffer with bits of word, LSB first
    s[i] = pop_bit(bit_buffer)
```

> **Note on tempering:** for generators with a $\mathrm{GF}(2)$-linear output transformation (MT
> tempering, MELG lagged tempering), the minimal polynomial of the tempered sequence equals
> that of the raw state — run the generator normally, no special handling needed.
>
> **Note on nonlinear scramblers:** for xorshift\*, xoshiro\*\*, xoroshiro+, the scrambled
> output does not satisfy a $\mathrm{GF}(2)$ linear recurrence. Apply BM to the **raw state bits**
> before the scrambler. See Section 8.

### 5.2 Berlekamp-Massey over $\mathrm{GF}(2)$

BM maintains:

| Variable | Type | Meaning |
|---|---|---|
| `L` | int | Current linear complexity $= \deg(m(x))$ so far |
| `C[0..k]` | bit array | Current connection polynomial, `C[0] = 1` |
| `B[0..k]` | bit array | Previous connection polynomial |
| `m` | int | Steps since last length increase |

**Initialise:** `L = 0`, `C = [1]`, `B = [1]`, `m = 1`.

**For each bit `s[n]`, `n = 0, 1, …, N-1`:**

```
1. Compute discrepancy:
       d = s[n] XOR popcount(C[1..L] AND s[n-1..n-L]) mod 2

2. If d == 0:
       m += 1

3. If d == 1:
       T = copy(C)
       for i = m to deg(B)+m:
           C[i] ^= B[i-m]          // C = C XOR x^m * B

       if 2*L <= n:
           L = n + 1 - L
           B = T
           m = 1
       else:
           m += 1
```

All arithmetic is XOR and AND over $\mathrm{GF}(2)$. No modular inverses needed.

**Polynomial shift:** `C[i] ^= B[i-m]` for `i = m .. deg(B)+m`, with `B[j] = 0`
for `j < 0`.

### 5.3 Filter decision

```
if L < k:
    → REJECT: P(x) is definitely reducible for this parameter set
    → (optional) report the factor degree L — it equals deg(m(x))
    → return REJECT immediately, no need to try more seeds

if L == k:
    → This seed does not witness reducibility
    → Continue to next seed
```

After all t seeds give L = k:

```
→ No witness of reducibility found
→ P(x) MAY be irreducible — proceed to Stage 2 (IterIrredTest)
→ DO NOT conclude irreducibility or primitivity from this alone
```

### 5.4 Intermediate plateau detection (diagnostic only)

If during a BM run $L$ stabilises at $\ell$ with $0 < \ell < k$ for more than $2\ell$ consecutive steps:

```
→ Strong evidence of an irreducible factor of degree ℓ
→ Useful for understanding which part of the parameter search went wrong
```

---

## 6. Stage 2: Exact Irreducibility Test

### 6.1 Using NTL

```cpp
#include <NTL/GF2X.h>
using namespace NTL;

// P must be a degree-k polynomial over GF(2) stored as GF2X
bool is_irreducible = IterIrredTest(P);

if (is_irreducible) {
    if (k_is_mersenne_prime_exponent(k)) {
        // P(x) is primitive — no order check needed
        accept(P);
    } else {
        // Must additionally verify order is 2^k − 1
        // Requires factoring 2^k − 1 and checking x^((2^k-1)/q) != 1 mod P(x)
        // for each prime factor q of 2^k − 1
        if (order_check(P, k))
            accept(P);
        else
            reject(P);   // irreducible but not primitive
    }
} else {
    reject(P);           // reducible — LCP filter missed it
}
```

### 6.2 Complexity of the exact test

| Algorithm | Cost |
|---|---|
| Rabin / `IterIrredTest` (naive) | $O(k^3)$ |
| Rabin / `IterIrredTest` (Karatsuba, NTL default for GF2X) | $O(k^{2.585})$ |
| Order check for Mersenne exponent $k$ | **Free — not needed** |
| Order check for general $k$ | $O(k^3) +$ cost of factoring $2^k - 1$ |

For $k = 607$ with NTL's Karatsuba GF2X: approximately $607^{2.585} \approx 5 \times 10^7$ bit operations.

---

## 7. Combined Complexity

```
Stage 1 (LCP filter, t seeds):   O(t · k²) bit operations
Stage 2 (IterIrredTest):          O(k^2.585) bit operations  [NTL Karatsuba]

Total per accepted candidate:     O(t · k²) + O(k^2.585)
Total per rejected candidate:     O(t · k²) only [Stage 2 skipped]
```

The filter's value is that for most candidates (reducible or wrong period), Stage 2 is never
reached. In a parameter search with rejection rate $r$:

$$\text{Expected cost per search step} = t \cdot k^2 + (1 - r) \cdot k^{2.585}$$

For $r = 0.99$ (99% of candidates rejected by Stage 1) and $k = 607$:
- LCP filter: $40 \times 607^2 \approx 1.5 \times 10^7$
- IterIrredTest (1% of the time): $0.01 \times 5 \times 10^7 \approx 5 \times 10^5$
- Total per step: $\approx 1.5 \times 10^7$

Compared to running IterIrredTest alone on every candidate: $5 \times 10^7$ per step.
**Speedup from the filter: approximately 3× in this scenario**, and better when the
rejection rate is higher.

---

## 8. Applying BM to the Raw State (Bypassing Scramblers)

For generators with nonlinear output scramblers (xorshift\*, xoshiro\*\*, xoroshiro+),
the scrambled output does not satisfy a $\mathrm{GF}(2)$ linear recurrence. BM on the scrambled
output will not find $P(x)$. Apply BM to the raw state instead:

```cpp
uint64_t next_raw_state_bit(generator_state* g, int bit_position) {
    advance_state(g);                               // state update (GF(2)-linear)
    return (g->state_word >> bit_position) & 1ULL;  // no scrambler
}
```

For generators with linear scramblers (MT tempering, MELG lagged tempering), either the
raw state or the output can be used — both yield the same minimal polynomial.

---

## 9. C++ Implementation Notes

### Polynomial representation

```cpp
static const int WORDS = (k + 63) / 64;
uint64_t C[WORDS] = {};   C[0] = 1;
uint64_t B[WORDS] = {};   B[0] = 1;
```

### Discrepancy computation

```cpp
int discrepancy(const uint64_t* C, int L, const uint64_t* seq_buf, int n) {
    int d = bit_at(seq_buf, n);
    for (int i = 1; i <= L; i++)
        d ^= C_bit(C, i) & bit_at(seq_buf, n - i);
    return d & 1;
}
```

For better performance, pack into 64-bit words with AND + `__builtin_popcountll`, XOR
the parities.

### Polynomial shift

```cpp
void xor_shifted(uint64_t* C, const uint64_t* B, int m, int degB) {
    int wshift = m / 64, bshift = m % 64;
    for (int i = 0; i <= (degB + 63) / 64; i++) {
        C[i + wshift] ^= B[i] << bshift;
        if (bshift > 0)
            C[i + wshift + 1] ^= B[i] >> (64 - bshift);
    }
}
```

### Circular sequence buffer

BM needs bits $s[n]$ and $s[n-1..n-L]$ with $L \le k$. A circular buffer of $k+1$ bits suffices:

```cpp
uint64_t seq_buf[(k + 64) / 64] = {};
// write s[n]: set bit at (n % (k+1))
// read s[n-j]: read bit at ((n-j+k+1) % (k+1))
```

### Bit extraction

```cpp
uint64_t word = genrand64_int64();
for (int bit = 0; bit < 64; bit++) {
    int s = (word >> bit) & 1;
    bm_step(s);   // feed into BM
}
```

Extracting all 64 bits per word reduces generator calls by 64×.

### PCLMUL acceleration

On x86 with CLMUL, enable with `-mpclmul` and `NTL_PCLMUL=on`. At $k = 607$ (~10 machine
words), Karatsuba + PCLMUL is faster than any FFT-based approach; FFT only wins for
$k \gtrsim 10{,}000$.

---

## 10. Seed Generation for the LCP Filter

Seeds must be nonzero, independent, and uniformly random over $\mathrm{GF}(2)^k \setminus \{0\}$. Use a
separate generator (e.g., splitmix64 seeded from hardware clock) to produce $k$ bits per
seed. Reject all-zero seeds.

The purpose of $t = 20\text{–}40$ seeds is to increase the probability of catching reducibility
when $d_\text{min}$ is small (few-bit factors). For large balanced factors the filter provides no
benefit — those are caught by Stage 2. The filter is most effective when the parameter
search naturally produces polynomials with occasional small factors, which is the common
case in practice.

---

## 11. Decision Table

| Observation | Conclusion | Certainty |
|---|---|---|
| BM gives $L < k$ for any seed | $P(x)$ **reducible** | **Certain** |
| BM gives intermediate plateau at $\ell$ | Factor of degree $\ell$ exists | **Certain** |
| BM gives $L = k$ for all $t$ seeds | $P(x)$ **may or may not** be irreducible | **No conclusion** |
| `IterIrredTest` returns true | $P(x)$ is **irreducible** | **Certain** |
| `IterIrredTest` true + $k$ Mersenne exponent | $P(x)$ is **primitive** | **Certain** |
| `IterIrredTest` true + order check passes | $P(x)$ is **primitive** | **Certain** |

---

## 12. Recommended Parameters

| $k$ | Mersenne prime exponent | $t$ (LCP seeds) | $N$ (bits) | Order check after IterIrredTest? |
|---|---|---|---|---|
| 607 | Yes | 20–40 | 1278 | Not needed |
| 1279 | Yes | 20–40 | 2622 | Not needed |
| 2281 | Yes | 20–40 | 4626 | Not needed |
| 4253 | Yes | 20–40 | 8570 | Not needed |
| 19937 | Yes | 10–20 | 39938 | Not needed |
| 44497 | Yes | 10–20 | 89058 | Not needed |
| other | No | 20–40 | $2k + 64$ | Required |

---

## 13. References

- Berlekamp, E. R. (1968). *Algebraic Coding Theory*. McGraw-Hill.
- Massey, J. L. (1969). Shift-register synthesis and BCH decoding. *IEEE Trans. Inf. Theory*,
  15(1), 122–127.
- Matsumoto, M., & Kurita, Y. (1992). Twisted GFSR generators. *ACM TOMACS*, 2(3), 179–194.
- Matsumoto, M., & Nishimura, T. (1998). Mersenne twister. *ACM TOMACS*, 8(1), 3–30.
- L'Ecuyer, P., & Panneton, F. (2009). F2-linear random number generators. In *Advancing
  the Frontiers of Simulation*, Springer.
- Harase, S., & Kimoto, T. (2018). Implementing 64-bit maximally equidistributed F2-linear
  generators. *ACM TOMS*, 44(3).
- Shoup, V. NTL: A Library for doing Number Theory. https://shoup.net/ntl/