# Dual Lattice Method for Equidistribution Testing

## Couture and L'Ecuyer (2000)

### References

- Couture, R. and L'Ecuyer, P. (2000). "Lattice Computations for Random Numbers."
  *Mathematics of Computation*, 69(230):757--765.
- L'Ecuyer, P. and Panneton, F. (2009). "F2-Linear Random Number Generators."
  In: *Advancing the Frontiers of Simulation*, Springer, pp. 169--193.

---

## 1. Problem Statement

Given a combined F2-linear generator with total state dimension k bits and
output resolution L (bits per output word), compute the **resolution gaps**
(dimension defects):

$$\delta_\ell = \left\lfloor k / \ell \right\rfloor - t_\ell, \qquad \ell = 1, \ldots, L$$

where $t_\ell$ is the **dimension of equidistribution with $\ell$-bit accuracy**:
the largest integer $t$ such that the set of $t$-tuples
$(u_i, u_{i+1}, \ldots, u_{i+t-1})$ truncated to $\ell$ bits each covers all
$2^{t\ell}$ cells equally (each cell visited $2^{k - t\ell}$ times over one full period).

The total dimension defect is $\Delta = \sum_{\ell=1}^{L} \delta_\ell$.
A generator is **maximally equidistributed (ME)** when $\Delta = 0$.

### Equivalent matrix characterisation

Form the $t\ell \times k$ binary matrix $M_{t,\ell}$ whose rows are the first
$\ell$ output-bit functions applied to states $A^0 x_0, A^1 x_0, \ldots,
A^{t-1} x_0$. Then $t_\ell = \max\{t : \operatorname{rank}(M_{t,\ell}) = t\ell\}$.

---

## 2. Algebraic Setting

### Field of formal Laurent series

Let $\mathbb{K} = \mathbb{F}_2((z^{-1}))$ be the field of formal Laurent
series with coefficients in $\mathbb{F}_2$:

$$x = \sum_{i=-\infty}^{n} \alpha_i z^i, \qquad \alpha_n \ne 0$$

The **norm** (or degree) of $x \ne 0$ is $|x| = n$ (the exponent of the
leading term).  The norm of zero is $-\infty$.

### Generating polynomials

For each output bit position $j = 0, \ldots, L-1$, there exists a polynomial
$g_j(z) \in \mathbb{F}_2[z]$ with $\deg(g_j) < k$ such that the formal power
series of the $j$-th output bit sequence is:

$$y_j(z) = \frac{g_j(z)}{P(z)}$$

where $P(z)$ is the characteristic polynomial of degree $k$.

**Construction of $g_j$:**

1. Initialize the generator with the canonical state $e_0 = (1, 0, \ldots, 0)$.
2. Run for $k$ steps, collecting the $j$-th output bit at each step $i = 0, \ldots, k-1$.
   Call these bits $a_j[i]$.
3. If the generator has a tempering transformation, apply it before reading bits.
4. Compute:

$$g_j(z) = \left(\sum_{i=0}^{k-1} a_j[i] \cdot z^{i+1}\right) \cdot P(z) \bmod z^k$$

In practice: for each $i$ where $a_j[i] = 1$, XOR $P(z) \cdot z^{i+1}$ into an
accumulator, then mask to the first $k$ bits (coefficients of $z^0$ through $z^{k-1}$).

### Normalization

Compute $g_0^{-1} \bmod P(z)$ using the extended GCD algorithm.  Define:

$$h_0(z) = 1, \qquad h_j(z) = g_j(z) \cdot g_0^{-1} \bmod P(z), \quad j = 1, \ldots, L-1$$

**Degenerate case:** if $\gcd(g_0, P) \ne 1$, divide out the common factor from
$P$ and all $g_j$ before inverting.  This reduces the effective degree $k$.

---

## 3. The Dual Lattice

For resolution $\ell$ (with $\ell \le L$), the **dual lattice** $L_\ell'$ is an
$\ell$-dimensional lattice over $\mathbb{F}_2[z]$ with basis:

$$H_\ell = \begin{pmatrix}
P(z) & 0      & 0      & \cdots & 0 \\
h_1(z) & 1    & 0      & \cdots & 0 \\
h_2(z) & 0    & 1      & \cdots & 0 \\
\vdots &       &        & \ddots & \vdots \\
h_{\ell-1}(z) & 0 & 0 & \cdots & 1
\end{pmatrix}$$

Row $i$ is the $i$-th basis vector $\mathbf{b}_i$, an $\ell$-tuple of polynomials.
The first column contains $P$ and the normalized generating polynomials $h_j$;
the remaining columns form an identity block.

### Key theorem (Couture-L'Ecuyer 2000, Theorem 1)

The dimension of equidistribution with $\ell$-bit accuracy is:

$$t_\ell = \min\left(\lambda_1(L_\ell'),\; \lfloor k/\ell \rfloor\right)$$

where $\lambda_1(L_\ell')$ is the **first successive minimum** of $L_\ell'$:
the degree (norm) of the shortest nonzero vector in $L_\ell'$.

Therefore the resolution gap is:

$$\delta_\ell = \lfloor k/\ell \rfloor - \min\left(\lambda_1(L_\ell'),\; \lfloor k/\ell \rfloor\right)$$

### Incremental basis construction

The basis for resolution $\ell$ is grown to resolution $\ell + 1$ by:

1. Zero out coordinate $\ell$ in all existing rows.
2. Append a new row: $(h_\ell(z),\; 0,\; \ldots,\; 0,\; 1)$ with the $1$ in
   coordinate $\ell$.
3. Re-run Lenstra's reduction on the expanded basis.

---

## 4. Data Structures

### PolVect (polynomial vector = one basis row)

A basis vector $\mathbf{b} = (p_0(z), p_1(z), \ldots, p_{\ell-1}(z))$ where each
$p_j(z)$ is a polynomial over $\mathbb{F}_2$ of degree $\le k$.  Stored as:

| Field | Type | Meaning |
|-------|------|---------|
| `coeffs[j]` | BitVect (k+1 bits) | Polynomial for coordinate $j$: bit $i$ = coefficient of $z^i$ |
| `deg` | int | $\|\mathbf{b}\| = \max_j \deg(p_j)$, or $-\infty$ if zero |
| `indicemaxdeg` | int | Which coordinate $j$ achieves the max degree |

### PolyLatBase (the full basis)

| Field | Type | Meaning |
|-------|------|---------|
| `vect_[i]` | PolVect | $i$-th basis vector $\mathbf{b}_i$ |
| `perm_[j]` | int | Permutation of coordinates (tracks pivots) |
| `invperm_[j]` | int | Inverse permutation |
| `resolution_` | int | Current number of coordinates $\ell$ |
| `degmax_` | int | Maximum possible degree $= k$ |

---

## 5. Lenstra's Lattice Reduction Algorithm

**Input:** Basis $\{\mathbf{b}_0, \ldots, \mathbf{b}_{\ell-1}\}$ of $L_\ell'$.

**Output:** $\|\mathbf{b}_0\|$ = first successive minimum = $\lambda_1(L_\ell')$.

**Invariant:** at step $m$, vectors $\mathbf{b}_0, \ldots, \mathbf{b}_m$ are
reduced (sorted by non-decreasing degree, with leading coefficients eliminated).

```
LENSTRA(basis, ell):
    m := -1

    WHILE m < ell - 1:
        // Step 1: RENUMBER — bring the shortest unreduced vector to position m+1
        j* := argmin_{j in {m+1, ..., ell-1}} ||b_j||
        swap b_{m+1} with b_{j*}

        // Step 2: SOLVE — find coefficients to cancel the leading term of b_{m+1}
        IF m >= 0:
            // Build the (m+1) x (m+2) binary system:
            //   A[j][i] = coefficient of z^{||b_i||} in coordinate perm[j] of b_i
            //   RHS[j]  = coefficient of z^{||b_{m+1}||} in coordinate perm[j] of b_{m+1}
            // The system is upper triangular (by construction from permute_coord).
            // Solve by forward substitution: x[0], x[1], ..., x[m] in F_2.
            solve_axb(m, x)

        // Step 3: FORM SUM — accumulate z^s * b_i for each x[i] = 1
        sum := 0 (zero polynomial vector)
        FOR i = 0 TO m:
            IF x[i] = 1:
                s := ||b_{m+1}|| - ||b_i||
                sum := sum + z^s * b_i      // XOR with shifted vector

        // Step 4: REDUCE
        IF sum = 0:
            temp := b_{m+1}                  // no reduction possible
            swap_flag := 1
        ELSE:
            temp := b_{m+1} + sum            // XOR = addition over F_2
            swap_flag := 0

        // Step 5: CHECK RESULT
        IF ||temp|| = ||b_{m+1}||:
            // Same norm — reduction consumed the leading term.
            // Accept and advance.
            b_{m+1} := temp
            permute_coord(m)                 // record which coordinate has the new max
            m := m + 1

        ELSE IF ||temp|| < ||b_{m+1}||:
            // Norm decreased — restart from an earlier m.
            b_{m+1} := temp
            newm := max{L in {0,...,m} : ||b_L|| <= ||temp||}
            IF no such L exists:
                m := -1
            ELSE:
                m := newm

        ELSE:  // ||temp|| > ||b_{m+1}|| (only when swap_flag = 1)
            b_{m+1} := temp                  // restore (swap back)

    RETURN ||b_0||
```

### Sub-procedures

**renumber(m, ell):** Find the vector of smallest norm among $\mathbf{b}_{m+1},
\ldots, \mathbf{b}_{\ell-1}$ and swap it into position $m+1$.

**solve_axb(m, x):** Solve a small $(m+1) \times (m+2)$ binary system by forward
substitution.  Row $j$, column $i$: the bit at degree $\|\mathbf{b}_i\|$ of
$\mathbf{b}_i$ in permuted coordinate $\operatorname{perm}[j]$.  The last column
(RHS) comes from $\mathbf{b}_{m+1}$.

**permute_coord(m):** After reducing $\mathbf{b}_{m+1}$, record which coordinate
achieves the maximum degree.  Update `perm_` and `invperm_`.

### Polynomial vector operations

| Operation | Formula | Cost |
|-----------|---------|------|
| `multbys(r, a, S)` | $\mathbf{r} = z^S \cdot \mathbf{a}$ (shift all coordinates by $S$) | $O(\ell \cdot k/64)$ |
| `add_polvect(r, a, b)` | $\mathbf{r} = \mathbf{a} + \mathbf{b}$ (XOR each coordinate) | $O(\ell \cdot k/64)$ |
| `multbys2(r, a, S)` | $\mathbf{r} \mathrel{+}= z^S \cdot \mathbf{a}$ (accumulated shift-XOR) | $O(\ell \cdot k/64)$ |

---

## 6. Full Algorithm

```
TestMELat(generators, transformations, L):
    // Step 1: Combined characteristic polynomial
    P(z) := product of individual characteristic polynomials over F_2[z]
    k    := deg(P)

    // Step 2: Generating polynomials
    FOR j = 0 TO min(L, k) - 1:
        Initialize all generators with canonical state e_0
        Run for k steps, collect bit j of the (tempered) combined output at each step
        g_j(z) := sum_{i: bit[j][i]=1} z^{i+1} * P(z)  (masked to k bits)

    // Step 3: Normalize
    Compute g_0^{-1} mod P(z) using extended GCD
    h_j(z) := g_j(z) * g_0^{-1} mod P(z)  for all j

    // Step 4: Initialize dual basis for resolution 1
    //   Row 0: (P(z))
    DualBase(basis, polys, P, 1)

    // Step 5: Loop over resolutions
    se := 0
    FOR ell = 1 TO min(L, k):
        lambda := Lenstra(basis, ell)              // first successive minimum
        t_ell  := min(lambda, floor(k / ell))
        delta[ell] := floor(k / ell) - t_ell
        se := se + delta[ell]

        // Early termination (optional)
        IF delta[ell] > delta_max[ell] OR se > se_max:
            BREAK

        IF ell < min(L, k):
            DualBaseIncrease(basis, polys)          // grow to resolution ell+1

    RETURN delta[], se
```

---

## 7. Complexity

### Per resolution $\ell$

- Lenstra's algorithm performs $O(\ell)$ iterations of the while loop (with
  possible restarts).  Each iteration involves:
  - One renumber: $O(\ell)$ comparisons.
  - One solve_axb: $O(m^2)$ bit operations.
  - One shift + XOR of polynomial vectors: $O(\ell \cdot k / 64)$ word operations.
  - One norm update: $O(\ell \cdot k / 64)$ word operations.
- Total per Lenstra call: $O(\ell^2 \cdot k / 64)$ word operations.

### All resolutions $\ell = 1, \ldots, L$

$$\text{Total} = O\!\left(\sum_{\ell=1}^{L} \ell^2 \cdot \frac{k}{64}\right) = O\!\left(\frac{L^3 \cdot k}{64}\right)$$

### Preprocessing (Steps 1--3)

- Characteristic polynomial: $O(k^2)$ if using BM, or $O(w^3 \cdot \text{poly\_ops})$
  for algebraic formulas.
- Generating polynomials: $O(L \cdot k^2 / 64)$ word operations (the main cost:
  $k$ generator steps collecting $L$ bits each, then $L$ polynomial products).
- Normalization: $O(k \log k)$ for the extended GCD.

### Comparison with the matrix method

| Method | Cost | Memory |
|--------|------|--------|
| Matrix (Gaussian elimination) | $O(k^2 \cdot t_\ell)$ per resolution | $O(k \cdot t_\ell)$ |
| Lattice (Couture-L'Ecuyer) | $O(\ell^2 \cdot k / 64)$ per resolution | $O(\ell \cdot k)$ |

For $k = 19937$, $L = 64$: the matrix method builds a $19937 \times \sim 312$ matrix
per resolution and row-reduces it.  The lattice method works with a $64 \times 64$
basis of degree-19937 polynomials.  The lattice method is dramatically faster for
large $k$.
