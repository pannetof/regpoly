# Tausworthe Generators and QuickTaus

L'Ecuyer-style framework for Tausworthe random number generators: the general
bit-level linear recurrence, its matrix form, the QuickTaus fast implementation
for primitive trinomials, and the generalization of QuickTaus to primitive
polynomials with more than three nonzero terms.

---

## 1. The linear recurrence

A Tausworthe generator (Tausworthe, 1965) produces a bit sequence
$\{b_n\}_{n \ge 0}$ over $\mathrm{GF}(2)$ via

$$
b_n \;=\; a_1 b_{n-1} \oplus a_2 b_{n-2} \oplus \cdots \oplus a_k b_{n-k},
\qquad n \ge k,
$$

with $a_i \in \{0,1\}$ and $a_k = 1$. The **characteristic polynomial** is

$$
P(z) \;=\; z^k + a_1 z^{k-1} + \cdots + a_k \;\in\; \mathrm{GF}(2)[z].
$$

The sequence attains its maximal period $\rho = 2^k - 1$ iff $P(z)$ is
**primitive** over $\mathrm{GF}(2)$. Non-primitive choices give shorter periods,
sometimes much shorter; in practice $P$ is always chosen primitive.

---

## 2. Output construction

Real-valued outputs are formed by grouping the bit stream into $L$-bit
fractional words with step size $s$:

$$
u_n \;=\; \sum_{j=1}^{L} b_{ns + j - 1}\, 2^{-j} \;\in\; [0, 1).
$$

Parameters:

- $L$ — output precision (typically the machine word size, e.g. 32 or 64).
- $s$ — step size / decimation rate; the number of recurrence steps between
  the first bit of $u_n$ and the first bit of $u_{n+1}$.

Two requirements:

1. $\gcd(s,\, 2^k - 1) = 1$, so that decimating by $s$ preserves the full period.
2. $L \le k$, so each output word fits inside the current state window.

---

## 3. F2-linear matrix form

The Tausworthe generator is a special case of the general $\mathrm{GF}(2)$-linear
generator

$$
\mathbf{x}_n \;=\; A\,\mathbf{x}_{n-1} \;\in\; \mathrm{GF}(2)^k,
\qquad
\mathbf{y}_n \;=\; B\,\mathbf{x}_n \;\in\; \mathrm{GF}(2)^L,
\qquad
u_n \;=\; \sum_{j=1}^{L} y_{n,j}\, 2^{-j},
$$

with $A \in \mathrm{GF}(2)^{k \times k}$ and $B \in \mathrm{GF}(2)^{L \times k}$.
The period is $2^k - 1$ iff the minimal polynomial of $A$ is primitive of
degree $k$.

For a Tausworthe generator specifically:

- **State:** $\mathbf{x}_n = (b_{ns},\, b_{ns+1},\, \ldots,\, b_{ns+k-1})^\top$,
  a sliding window of $k$ consecutive bits.
- **Transition:** $A = C^s$, where $C$ is the companion matrix of $P(z)$.
  Advancing the output index by $1$ corresponds to advancing the bit index
  by $s$, i.e. applying $C$ $s$ times.
- **Output:** $B = [\,I_L \mid 0\,]$, extracting the first $L$ coordinates of
  the state window.

The minimal polynomial of $A = C^s$ equals that of $C$ (namely $P$) iff
$\gcd(s,\, 2^k - 1) = 1$, recovering the period condition stated in §2.

---

## 4. Equidistribution (brief)

A generator is **$(t, \ell)$-equidistributed** if, partitioning $[0,1)^t$ into
$2^{t\ell}$ identical cubic cells, every cell contains exactly $2^{k - t\ell}$
of the $2^k - 1$ output tuples $(u_n, u_{n+1}, \ldots, u_{n+t-1})$ taken over
one full period (plus the zero state, giving exact equality).

Equivalently, the $k \times t\ell$ matrix over $\mathrm{GF}(2)$ whose rows come
from the first $\ell$ output bits at $t$ successive times — built from
$B, BA, BA^2, \ldots, BA^{t-1}$ restricted to their first $\ell$ rows — has
full rank $t\ell$.

A generator is **maximally equidistributed (ME)** if it is
$(t, \lfloor k/t \rfloor)$-equidistributed for every $t \ge 1$. Single-component
Tausworthe generators with sparse characteristic polynomial are never ME —
this is the main motivation for the combined constructions of §7.

---

## 5. The trinomial case

Take a **primitive trinomial**

$$
P(z) \;=\; z^k + z^q + 1, \qquad 0 < q < k.
$$

By the reciprocal-polynomial argument ($P$ primitive $\iff z^k P(1/z)$
primitive), we may assume $q < k/2$, i.e. $0 < 2q < k$, without loss of
generality. The corresponding bit-level recurrence is

$$
b_n \;=\; b_{n-k+q} \;\oplus\; b_{n-k}.
$$

### 5.1 The QuickTaus step

Store the $k$-bit state left-justified in a machine word $z$ of size
$L_w \ge k$, with the oldest bit at the MSB side:

$$
z \;=\; \underbrace{b_{n-k}\; b_{n-k+1}\; \cdots\; b_{n-1}}_{k\ \text{bits}}\;\;
       \underbrace{0\;0\;\cdots\;0}_{L_w - k\ \text{bits}}.
$$

Advancing the output by one step (i.e., the bit index by $s$) is then
implemented by:

```c
b = ((z << q) ^ z) >> (L_w - s);         /* s new bits, right-justified */
z = ((z & mask) << s) ^ b;                /* shift state, inject new bits */
```

with `mask` chosen to zero out the $s$ bits about to be shifted off the top.
The first line computes $b_n, b_{n+1}, \ldots, b_{n+s-1}$ in parallel via a
single shift-XOR-shift: shifting $z$ left by $q$ aligns the $b_{\cdot-k+q}$
tap with the $b_{\cdot-k}$ tap, XOR applies the recurrence to all $s$
positions at once, and the final right-shift extracts the $s$ newest bits.

### 5.2 The conditions

The step above is valid iff

$$
\boxed{\;0 < 2q < k, \qquad 0 < s \le k - q, \qquad k \le L_w.\;}
$$

- $\mathbf{0 < 2q < k}$: normalization — $q$ is the smaller of the two taps,
  $k - q$ the larger. (A trinomial has two interior "gaps" of size $q$ and
  $k-q$; picking $q < k/2$ makes $k - q$ the binding gap for what follows.)
- $\mathbf{s \le k - q}$: ensures all $s$ new bits can be produced from the
  current state window without referring to bits that have not yet been
  computed. To compute $b_{n+i}$ we need $b_{n+i-k+q}$, which sits at offset
  $q + i$ in the state window; requiring $q + (s-1) \le k - 1$ gives
  $s \le k - q$. The other operand $b_{n+i-k}$ (at offset $i$) is
  automatically in range.
- $\mathbf{k \le L_w}$: the whole state fits in one machine word, so shift
  operations act on the full state.

### 5.3 Equidistribution is weak

A single trinomial QuickTaus component is never maximally equidistributed.
Because $P$ has only three nonzero terms, its lattice in
$\mathrm{GF}(2, x)^k$ is anisotropic and the generator has predictable defects
in high dimensions (L'Ecuyer, 1996). This is the structural reason behind the
combined constructions of §7.

---

## 6. Generalization to polynomials with $t > 3$ terms

Let

$$
P(z) \;=\; z^k + z^{q_{t-2}} + z^{q_{t-3}} + \cdots + z^{q_1} + 1,
\qquad 0 < q_1 < q_2 < \cdots < q_{t-2} < k,
$$

be a primitive polynomial of degree $k$ with $t$ nonzero terms ($t-2$ interior
exponents). The recurrence becomes

$$
b_n \;=\; b_{n-k+q_{t-2}} \,\oplus\, b_{n-k+q_{t-3}} \,\oplus\, \cdots \,\oplus\,
         b_{n-k+q_1} \,\oplus\, b_{n-k}.
$$

### 6.1 The general QuickTaus step

```c
b = (z ^ (z << q_1) ^ (z << q_2) ^ ... ^ (z << q_{t-2})) >> (L_w - s);
z = ((z & mask) << s) ^ b;
```

One shift-XOR per interior tap, then a final right-shift to extract the new
bits, then the state update. The trinomial case ($t = 3$) is recovered with
a single interior tap $q_1 = q$.

### 6.2 The generalized condition

$$
\boxed{\;0 < s \le k - q_{t-2} < k \le L_w,\;}
$$

with $P$ primitive.

**The binding tap is the largest interior exponent $q_{t-2}$** — the one
closest to $k$. All smaller interior taps $q_j$ ($j < t-2$) are automatically
in range once the largest is.

**Derivation.** To compute $b_{n+i}$ for $i = 0, \ldots, s-1$ we need
$b_{n+i-k+q_j}$ for each $j = 1, \ldots, t-2$, located at state offset
$q_j + i$. The tightest constraint comes from the largest $q_j$:

$$
q_{t-2} + (s - 1) \;\le\; k - 1 \quad\iff\quad s \;\le\; k - q_{t-2}.
$$

The tap at offset $i$ (coming from $b_{n-k}$) is in range automatically since
$i \le s - 1 \le k - q_{t-2} - 1 \le k - 1$.

### 6.3 Comments on the trinomial normalization

The "$2q < k$" normalization of the trinomial case has **no clean analogue**
for $t > 3$: there is no single swap that plays the role of the reciprocal
polynomial trick. One can still WLOG reduce to $q_{t-2} \le k/2$ by comparing
$P(z)$ with its reciprocal $z^k P(1/z)$ (both primitive or both not), and
choosing whichever has the smaller largest interior exponent. Doing so makes
$k - q_{t-2}$ as large as possible, maximizing the achievable step size $s$.

### 6.4 Practical cost

Each step requires $t - 1$ shift-XOR pairs and produces at most
$s \le k - q_{t-2}$ bits. Adding terms to $P$ tends to push $q_{t-2}$ closer
to $k$ (denser polynomials have less room at the top), which shrinks $s$ and
hurts throughput.

This trade-off is why L'Ecuyer's designs stick to sparse polynomials
(trinomials, occasionally pentanomials) and achieve equidistribution by
**combining** components (§7) rather than by using a single dense polynomial.
Dense $A$ matrices with cheap word-level implementations are instead obtained
by dropping the "polynomial with taps" framing entirely — this is the route
taken by the WELL family (Panneton–L'Ecuyer–Matsumoto, 2006), where $A$
factors into a short product of sparse elementary matrices without any of its
individual taps being visible at the polynomial level.

---

## 7. Combined Tausworthe generators

L'Ecuyer's practical designs XOR-combine $J$ independent Tausworthe
components, each implemented via QuickTaus on a primitive trinomial:

$$
u_n \;=\; u_n^{(1)} \,\oplus\, u_n^{(2)} \,\oplus\, \cdots \,\oplus\, u_n^{(J)}.
$$

The combined output is equivalent (Tezuka–L'Ecuyer, 1991) to a single
Tausworthe generator with **reducible** characteristic polynomial

$$
P(z) \;=\; P_1(z)\, P_2(z) \,\cdots\, P_J(z)
$$

of degree $k = k_1 + k_2 + \cdots + k_J$. The component periods
$2^{k_j} - 1$ are chosen pairwise coprime so the combined period equals their
product.

Key instances from L'Ecuyer's tables:

| Name | $J$ | word $L_w$ | component degrees $(k_j)$ | period |
|------|:---:|:----------:|---------------------------|--------|
| Combined-3 (Tezuka–L'Ecuyer, 1991) | 3 | 32 | (31, 29, 28) | $\approx 2^{88}$ |
| **LFSR113** (L'Ecuyer, 1999) | 4 | 32 | (31, 29, 28, 25) | $\approx 2^{113}$ |
| **LFSR258** (L'Ecuyer, 1999) | 5 | 64 | (63, 55, 52, 47, 41) | $\approx 2^{258}$ |

LFSR113 and LFSR258 are maximally equidistributed, found by exhaustive search
over trinomial parameter tuples $(k_j, q_j, s_j)$ satisfying the QuickTaus
constraints of §5.2, subject to a coprimality condition on the $k_j$ and the
ME criterion on the combined lattice.

---

## 8. Back to the matrix form

A QuickTaus component with parameters $(k, q, s)$ has matrix representation

$$
A \;=\; C_{P}^{\,s}, \qquad B \;=\; \bigl[\, I_L \mid 0 \,\bigr],
$$

where $C_P$ is the companion matrix of $P(z) = z^k + z^q + 1$ and $L$ is the
output precision. $A$ is full-rank and has primitive minimal polynomial $P$;
$B$ is the top-$L$-row selector.

A combined generator with $J$ QuickTaus components has

$$
A \;=\; \mathrm{diag}(A_1, \ldots, A_J) \;\in\; \mathrm{GF}(2)^{k \times k},
\qquad
B \;=\; [\,B_1 \mid B_2 \mid \cdots \mid B_J\,] \;\in\; \mathrm{GF}(2)^{L \times k},
$$

with state $\mathbf{x} = (\mathbf{x}^{(1)}, \ldots, \mathbf{x}^{(J)})$
concatenated across components and $B$ implementing the component-wise XOR.
This is the form in which equidistribution tests, lattice-structure
computations in $\mathrm{GF}(2, x)^k$, and observability-type analyses (e.g.
for locally antithetic output design) are all carried out.

---

## References

- R. C. Tausworthe. *Random numbers generated by linear recurrence modulo two.*
  Math. Comp. **19** (1965), 201–209.
  [DOI](https://doi.org/10.1090/S0025-5718-1965-0184406-1)
- S. Tezuka and P. L'Ecuyer. *Efficient and portable combined Tausworthe random
  number generators.* ACM Trans. Model. Comput. Simul. **1** (1991), 99–112.
  [DOI](https://doi.org/10.1145/122124.122127).  Published as
  [`taus88`](/library/taus88) in the library.
- P. L'Ecuyer. *Maximally equidistributed combined Tausworthe generators.*
  Math. Comp. **65** (1996), 203–213.
  [DOI](https://doi.org/10.1090/S0025-5718-96-00696-5)
- P. L'Ecuyer. *Tables of maximally equidistributed combined LFSR generators.*
  Math. Comp. **68** (1999), 261–269.
  [JSTOR](https://www.jstor.org/stable/2585109).  Published as
  [`lfsr113`](/library/lfsr113) in the library.
- F. Panneton, P. L'Ecuyer, and M. Matsumoto. *Improved long-period generators
  based on linear recurrences modulo 2.* ACM Trans. Math. Softw. **32** (2006),
  1–16. [WELL family.]
  [DOI](https://doi.org/10.1145/1132973.1132974)
- P. L'Ecuyer and F. Panneton. *F<sub>2</sub>-linear random number generators.*
  In *Advancing the Frontiers of Simulation* (2009), Springer, 169–193.
  [DOI](https://doi.org/10.1007/978-1-4419-0652-0_8).
  [Survey of the general framework.]