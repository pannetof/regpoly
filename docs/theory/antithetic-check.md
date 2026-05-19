# Algorithm: Checking Local Antitheticity of a Linear RNG Point Set

## Overview

Given a linear random number generator defined by matrices $A$ and $B$ over
$\mathbb{F}_2$, this algorithm determines whether the generated point set is locally
antithetic in dimension $d$.

---

## Symbol Definitions

| Symbol | Type | Description |
|--------|------|-------------|
| $k$ | positive integer | State size in bits |
| $L$ | positive integer | Output size in bits, $L \le k$ |
| $d$ | positive integer | Dimension to test, $d \ge 1$ |
| $A$ | matrix in $\mathbb{F}_2^{k \times k}$ | Transition matrix of the recurrence $\mathbf{x}_n = A \cdot \mathbf{x}_{n-1}$ |
| $B$ | matrix in $\mathbb{F}_2^{L \times k}$ | Output matrix, $\mathbf{y}_n = B \cdot \mathbf{x}_n$ |
| $C_j$ | matrix in $\mathbb{F}_2^{L \times k}$ | Output matrix for coordinate $j$, defined as $C_j = B \cdot A^{j-1}$ |
| $\mathbf{v}$ | column vector in $\mathbb{F}_2^L$ | Target reflection vector, defined as $\mathbf{v} = (1, 1, \ldots, 1)^T$ (all ones) |
| $\delta_j$ | column vector in $\mathbb{F}_2^k$ | Antithetic offset for coordinate $j$ |
| $\mathbf{0}_L$ | column vector in $\mathbb{F}_2^L$ | Zero vector of length $L$ |
| $D_\text{max}$ | positive integer | Maximum feasible dimension, defined as $D_\text{max} = \lfloor \sqrt{k} \rfloor$ |

All arithmetic is performed in $\mathbb{F}_2$: addition is XOR (`^`), multiplication
is AND (`&`). All matrices and vectors have entries in $\{0, 1\}$.

---

## What Local Antitheticity Means

The point set is **locally antithetic in dimension $d$** if and only if for
each coordinate index $j$ in $\{1, \ldots, d\}$ there exists a nonzero vector
$\delta_j$ in $\mathbb{F}_2^k$ such that:

$$\begin{aligned}
C_j \cdot \delta_j &= \mathbf{v} \quad &\text{(condition A: coordinate } j \text{ is reflected)} \\
C_i \cdot \delta_j &= \mathbf{0}_L \quad &\text{(condition B: all other coordinates } i \ne j \text{ are unchanged)}
\end{aligned}$$

for all $i$ in $\{1, \ldots, d\}$, $i \ne j$.

Here $C_j = B \cdot A^{j-1}$ is computed over $\mathbb{F}_2$:
- $C_1 = B$
- $C_2 = B \cdot A$
- $C_3 = B \cdot A^2$
- etc.

The vector $\mathbf{v}$ is the all-ones vector of length $L$:

$$\mathbf{v} = (1, 1, \ldots, 1)^T \quad \text{in } \mathbb{F}_2^L.$$

---

## Input

- `A`: a $k \times k$ matrix over $\mathbb{F}_2$, stored as an array of $k$ rows, each
  row being an array of $k$ bits (or equivalently $k$ integers in $\{0, 1\}$)
- `B`: an $L \times k$ matrix over $\mathbb{F}_2$, stored as an array of $L$ rows, each
  row being an array of $k$ bits
- `d`: a positive integer, the dimension to test; must satisfy $1 \le d \le D_\text{max}$
  where $D_\text{max} = \lfloor \sqrt{k} \rfloor$

---

## Output

- `is_antithetic`: boolean, `True` if the point set is locally antithetic in
  dimension $d$, `False` otherwise
- `offsets`: if `is_antithetic` is `True`, a list $[\delta_1, \ldots, \delta_d]$
  where each $\delta_j$ is a column vector in $\mathbb{F}_2^k$ (array of $k$ bits);
  if `is_antithetic` is `False`, return `None`

---

## Step-by-Step Algorithm

### Step 1: Validate Input

Check that:
- $A$ has shape $k \times k$
- $B$ has shape $L \times k$
- $L \le k$
- $1 \le d \le \lfloor \sqrt{k} \rfloor$

If any check fails, raise an error with a descriptive message.

---

### Step 2: Precompute the Output Matrices $C_1, \ldots, C_d$

Define $C_1 = B$ (copy the matrix $B$).

For $j = 2, 3, \ldots, d$:

$$C_j = C_{j-1} \cdot A \quad \text{(matrix multiplication over } \mathbb{F}_2\text{)}$$

Matrix multiplication over $\mathbb{F}_2$ of an $L \times k$ matrix $M$ by a $k \times k$
matrix $A$ produces an $L \times k$ matrix $P$ where:

$$P[i][m] = \bigoplus_{s \in \{0, \ldots, k-1\}} (M[i][s] \mathbin{\text{AND}} A[s][m])$$

for $i$ in $\{0, \ldots, L-1\}$ and $m$ in $\{0, \ldots, k-1\}$.

Store all $d$ matrices: $C[1], C[2], \ldots, C[d]$, each of shape $L \times k$.

---

### Step 3: For Each Coordinate $j$, Build and Solve a Linear System

For each $j$ in $\{1, 2, \ldots, d\}$:

#### 3a: Assemble the System Matrix and Right-Hand Side

Build the system matrix $S_j$ of shape $(dL) \times k$ over $\mathbb{F}_2$ by stacking
the output matrices vertically:

```
S_j = [ C[1]   ]    <- rows 0       to L-1
      [ C[2]   ]    <- rows L       to 2L-1
      [  ...   ]
      [ C[d]   ]    <- rows (d-1)*L to d*L-1
```

That is, $S_j[iL : (i+1)L, :] = C[i+1]$ for $i$ in $\{0, \ldots, d-1\}$.

Note that $S_j$ is the same matrix for all $j$ — it does not depend on $j$.
Call it $S$ (no subscript needed). It has shape $(dL) \times k$.

Build the right-hand side vector $\mathrm{rhs}_j$ of length $dL$ over $\mathbb{F}_2$ as
follows:

$$\mathrm{rhs}_j[iL : (i+1)L] = \begin{cases} \mathbf{v} & \text{if } i+1 = j \text{ (the } j\text{-th block is all ones)} \\ \mathbf{0}_L & \text{otherwise (all other blocks are zero)} \end{cases}$$

That is, $\mathrm{rhs}_j$ is zero everywhere except in the block of $L$ entries
corresponding to coordinate $j$, where it equals $\mathbf{v} = (1, 1, \ldots, 1)^T$.

#### 3b: Solve the Linear System Over $\mathbb{F}_2$

Solve the system:

$$S \cdot \delta_j = \mathrm{rhs}_j$$

where:
- $S$ has shape $(dL) \times k$
- $\delta_j$ is the unknown column vector of length $k$
- $\mathrm{rhs}_j$ is the right-hand side column vector of length $dL$
- All operations are over $\mathbb{F}_2$

This is an overdetermined system (more equations than unknowns when $dL > k$,
underdetermined when $dL < k$). Solve it using **Gaussian elimination over
$\mathbb{F}_2$** on the augmented matrix $[S \mid \mathrm{rhs}_j]$ of shape $(dL) \times (k+1)$.

**Gaussian elimination over $\mathbb{F}_2$:**

```
Augment: M = [S | rhs_j]   shape (d*L) x (k+1)

pivot_row = 0
for col = 0 to k-1:
    # Find a pivot in column col at or below pivot_row
    found = False
    for row = pivot_row to d*L - 1:
        if M[row][col] == 1:
            swap rows M[pivot_row] and M[row]
            found = True
            break
    if not found:
        continue   # No pivot in this column, move to next

    # Eliminate all other 1s in column col
    for row = 0 to d*L - 1:
        if row != pivot_row and M[row][col] == 1:
            M[row] = M[row] XOR M[pivot_row]   # row operation over GF(2)

    pivot_row = pivot_row + 1
```

After elimination:

**Consistency check:** For each row $r$ in $\{0, \ldots, dL - 1\}$ where all
entries $M[r][0], \ldots, M[r][k-1]$ are zero: if $M[r][k] = 1$, the system
is **inconsistent** (no solution exists).

**If inconsistent:** set `is_antithetic = False`, return immediately with
`offsets = None`. Do not proceed to the next coordinate $j$ or the next
dimension.

**If consistent:** extract any particular solution $\delta_j$ from the
row-reduced augmented matrix by back-substitution. Free variables (columns
with no pivot) may be set to 0.

#### 3c: Check Non-Zero Condition

Verify that the solution $\delta_j$ is not the zero vector:

$$\delta_j \ne \mathbf{0}_k \quad \text{(at least one entry of } \delta_j \text{ is 1)}$$

If $\delta_j = \mathbf{0}_k$, the system is consistent but the only solution is
trivial. This should not happen if $A$ has maximal period and $\mathbf{v} \ne 0$, but
check defensively. If it occurs, set `is_antithetic = False` and return.

---

### Step 4: All Systems Consistent

If all $d$ systems (one per coordinate $j$) were consistent and produced
nonzero solutions, then the point set is locally antithetic in dimension $d$.

Set:
```
is_antithetic = True
offsets = [delta_1, delta_2, ..., delta_d]
```

---

### Step 5: Return

Return `(is_antithetic, offsets)`.

---

## Worked Example for Verification

With $k = 4$, $L = 2$, $d = 1$:

```
A = [[0, 1, 0, 0],    (4 x 4 over GF(2))
     [0, 0, 1, 0],
     [0, 0, 0, 1],
     [1, 1, 0, 0]]

B = [[1, 0, 0, 0],    (2 x 4 over GF(2))
     [0, 1, 0, 0]]

v = [1, 1]^T          (L=2 ones)
```

$C_1 = B$ (shape $2 \times 4$).

For $d = 1$, $j = 1$:
- $S = C_1$ (shape $2 \times 4$)
- $\mathrm{rhs}_1 = \mathbf{v} = [1, 1]^T$

Augmented matrix:
```
[1, 0, 0, 0 | 1]
[0, 1, 0, 0 | 1]
```
Already in row echelon form. Solution: $\delta_1 = [1, 1, 0, 0]^T$
(free variables $\delta_1[2] = 0$, $\delta_1[3] = 0$).

Check: $B \cdot \delta_1 = [1, 1]^T = \mathbf{v}$. Correct.

`is_antithetic = True`, `offsets = [[1,1,0,0]^T]`.

---

## Implementation Notes

- All arithmetic must be in $\mathbb{F}_2$: use bitwise XOR for addition, bitwise AND
  for multiplication. No modular arithmetic with integers is needed if
  matrices are stored as arrays of 0/1 integers or as Python `int` bitmasks.
- For efficiency with large $k$ (e.g., $k = 800$), represent each row of a
  matrix as a Python `int` (bitmask of $k$ bits) and use bitwise operations
  directly. XOR two rows: `row_a ^ row_b`. Check if a bit is set:
  `(row >> col) & 1`. Set a bit: `row | (1 << col)`.
- The system matrix $S$ is the same for all $j$: assemble it once and reuse,
  only changing $\mathrm{rhs}_j$ for each $j$.
- The row-reduced form of $[S \mid \mathrm{rhs}_j]$ for different $j$ only differs in
  the last column. Therefore, row-reduce $S$ once and apply the same row
  operations to each $\mathrm{rhs}_j$ separately, saving significant computation.

---

## Complexity

| Step | Cost |
|------|------|
| Precompute $C_1, \ldots, C_d$ | $d$ matrix multiplications, each $O(L \cdot k^2 / 64)$ with bitmasks |
| Assemble $S$ | $O(d \cdot L \cdot k / 64)$ |
| Row-reduce $S$ (once) | $O(d \cdot L \cdot k^2 / 64)$ |
| Apply row ops to each $\mathrm{rhs}_j$ | $O(d \cdot d \cdot L / 64)$ per $j$, total $O(d^2 \cdot L / 64)$ |
| Total | $O(d \cdot L \cdot k^2 / 64)$ dominated by row reduction |