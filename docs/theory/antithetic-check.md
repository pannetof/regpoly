# Algorithm: Checking Local Antitheticity of a Linear RNG Point Set

## Overview

Given a linear random number generator defined by matrices `A` and `B` over
GF(2), this algorithm determines whether the generated point set is locally
antithetic in dimension `d`.

---

## Symbol Definitions

| Symbol | Type | Description |
|--------|------|-------------|
| `k` | positive integer | State size in bits |
| `L` | positive integer | Output size in bits, `L <= k` |
| `d` | positive integer | Dimension to test, `d >= 1` |
| `A` | matrix in GF(2)^{k x k} | Transition matrix of the recurrence `x_n = A * x_{n-1}` |
| `B` | matrix in GF(2)^{L x k} | Output matrix, `y_n = B * x_n` |
| `C_j` | matrix in GF(2)^{L x k} | Output matrix for coordinate `j`, defined as `C_j = B * A^{j-1}` |
| `v` | column vector in GF(2)^L | Target reflection vector, defined as `v = (1, 1, ..., 1)^T` (all ones) |
| `delta_j` | column vector in GF(2)^k | Antithetic offset for coordinate `j` |
| `0_L` | column vector in GF(2)^L | Zero vector of length `L` |
| `D_max` | positive integer | Maximum feasible dimension, defined as `D_max = floor(sqrt(k))` |

All arithmetic is performed in GF(2): addition is XOR (`^`), multiplication
is AND (`&`). All matrices and vectors have entries in {0, 1}.

---

## What Local Antitheticity Means

The point set is **locally antithetic in dimension `d`** if and only if for
each coordinate index `j` in `{1, ..., d}` there exists a nonzero vector
`delta_j` in GF(2)^k such that:

```
C_j * delta_j = v          (condition A: coordinate j is reflected)
C_i * delta_j = 0_L        (condition B: all other coordinates i != j are unchanged)
```

for all `i` in `{1, ..., d}`, `i != j`.

Here `C_j = B * A^{j-1}` is computed over GF(2):
- `C_1 = B`
- `C_2 = B * A`
- `C_3 = B * A^2`
- etc.

The vector `v` is the all-ones vector of length `L`:
```
v = (1, 1, ..., 1)^T   in GF(2)^L
```

---

## Input

- `A`: a `k x k` matrix over GF(2), stored as an array of `k` rows, each
  row being an array of `k` bits (or equivalently `k` integers in {0,1})
- `B`: an `L x k` matrix over GF(2), stored as an array of `L` rows, each
  row being an array of `k` bits
- `d`: a positive integer, the dimension to test; must satisfy `1 <= d <= D_max`
  where `D_max = floor(sqrt(k))`

---

## Output

- `is_antithetic`: boolean, `True` if the point set is locally antithetic in
  dimension `d`, `False` otherwise
- `offsets`: if `is_antithetic` is `True`, a list `[delta_1, ..., delta_d]`
  where each `delta_j` is a column vector in GF(2)^k (array of `k` bits);
  if `is_antithetic` is `False`, return `None`

---

## Step-by-Step Algorithm

### Step 1: Validate Input

Check that:
- `A` has shape `k x k`
- `B` has shape `L x k`
- `L <= k`
- `1 <= d <= floor(sqrt(k))`

If any check fails, raise an error with a descriptive message.

---

### Step 2: Precompute the Output Matrices C_1, ..., C_d

Define `C_1 = B` (copy the matrix `B`).

For `j = 2, 3, ..., d`:
```
C_j = C_{j-1} * A   (matrix multiplication over GF(2))
```

Matrix multiplication over GF(2) of an `L x k` matrix `M` by a `k x k`
matrix `A` produces an `L x k` matrix `P` where:
```
P[i][m] = XOR over s in {0,...,k-1} of (M[i][s] AND A[s][m])
```
for `i` in `{0, ..., L-1}` and `m` in `{0, ..., k-1}`.

Store all `d` matrices: `C[1], C[2], ..., C[d]`, each of shape `L x k`.

---

### Step 3: For Each Coordinate j, Build and Solve a Linear System

For each `j` in `{1, 2, ..., d}`:

#### 3a: Assemble the System Matrix and Right-Hand Side

Build the system matrix `S_j` of shape `(d*L) x k` over GF(2) by stacking
the output matrices vertically:

```
S_j = [ C[1]   ]    <- rows 0       to L-1
      [ C[2]   ]    <- rows L       to 2L-1
      [  ...   ]
      [ C[d]   ]    <- rows (d-1)*L to d*L-1
```

That is, `S_j[i*L : (i+1)*L, :]  = C[i+1]` for `i` in `{0, ..., d-1}`.

Note that `S_j` is the same matrix for all `j` — it does not depend on `j`.
Call it `S` (no subscript needed). It has shape `(d*L) x k`.

Build the right-hand side vector `rhs_j` of length `d*L` over GF(2) as
follows:

```
rhs_j[i*L : (i+1)*L] = v      if i+1 == j   (the j-th block is all ones)
                      = 0_L   otherwise       (all other blocks are zero)
```

That is, `rhs_j` is zero everywhere except in the block of `L` entries
corresponding to coordinate `j`, where it equals `v = (1,1,...,1)^T`.

#### 3b: Solve the Linear System Over GF(2)

Solve the system:
```
S * delta_j = rhs_j
```
where:
- `S` has shape `(d*L) x k`
- `delta_j` is the unknown column vector of length `k`
- `rhs_j` is the right-hand side column vector of length `d*L`
- All operations are over GF(2)

This is an overdetermined system (more equations than unknowns when `d*L > k`,
underdetermined when `d*L < k`). Solve it using **Gaussian elimination over
GF(2)** on the augmented matrix `[S | rhs_j]` of shape `(d*L) x (k+1)`.

**Gaussian elimination over GF(2):**

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

**Consistency check:** For each row `r` in `{0, ..., d*L - 1}` where all
entries `M[r][0], ..., M[r][k-1]` are zero: if `M[r][k] == 1`, the system
is **inconsistent** (no solution exists).

**If inconsistent:** set `is_antithetic = False`, return immediately with
`offsets = None`. Do not proceed to the next coordinate `j` or the next
dimension.

**If consistent:** extract any particular solution `delta_j` from the
row-reduced augmented matrix by back-substitution. Free variables (columns
with no pivot) may be set to 0.

#### 3c: Check Non-Zero Condition

Verify that the solution `delta_j` is not the zero vector:
```
delta_j != 0_k   (at least one entry of delta_j is 1)
```

If `delta_j == 0_k`, the system is consistent but the only solution is
trivial. This should not happen if `A` has maximal period and `v != 0`, but
check defensively. If it occurs, set `is_antithetic = False` and return.

---

### Step 4: All Systems Consistent

If all `d` systems (one per coordinate `j`) were consistent and produced
nonzero solutions, then the point set is locally antithetic in dimension `d`.

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

With `k = 4`, `L = 2`, `d = 1`:

```
A = [[0, 1, 0, 0],    (4 x 4 over GF(2))
     [0, 0, 1, 0],
     [0, 0, 0, 1],
     [1, 1, 0, 0]]

B = [[1, 0, 0, 0],    (2 x 4 over GF(2))
     [0, 1, 0, 0]]

v = [1, 1]^T          (L=2 ones)
```

`C_1 = B` (shape 2 x 4).

For `d=1`, `j=1`:
- `S = C_1` (shape 2 x 4)
- `rhs_1 = v = [1, 1]^T`

Augmented matrix:
```
[1, 0, 0, 0 | 1]
[0, 1, 0, 0 | 1]
```
Already in row echelon form. Solution: `delta_1 = [1, 1, 0, 0]^T`
(free variables `delta_1[2] = 0`, `delta_1[3] = 0`).

Check: `B * delta_1 = [1,1]^T = v`. Correct.

`is_antithetic = True`, `offsets = [[1,1,0,0]^T]`.

---

## Implementation Notes

- All arithmetic must be in GF(2): use bitwise XOR for addition, bitwise AND
  for multiplication. No modular arithmetic with integers is needed if
  matrices are stored as arrays of 0/1 integers or as Python `int` bitmasks.
- For efficiency with large `k` (e.g., `k = 800`), represent each row of a
  matrix as a Python `int` (bitmask of `k` bits) and use bitwise operations
  directly. XOR two rows: `row_a ^ row_b`. Check if a bit is set:
  `(row >> col) & 1`. Set a bit: `row | (1 << col)`.
- The system matrix `S` is the same for all `j`: assemble it once and reuse,
  only changing `rhs_j` for each `j`.
- The row-reduced form of `[S | rhs_j]` for different `j` only differs in
  the last column. Therefore, row-reduce `S` once and apply the same row
  operations to each `rhs_j` separately, saving significant computation.

---

## Complexity

| Step | Cost |
|------|------|
| Precompute `C_1,...,C_d` | `d` matrix multiplications, each `O(L * k^2 / 64)` with bitmasks |
| Assemble `S` | `O(d * L * k / 64)` |
| Row-reduce `S` (once) | `O(d * L * k^2 / 64)` |
| Apply row ops to each `rhs_j` | `O(d * d * L / 64)` per `j`, total `O(d^2 * L / 64)` |
| Total | `O(d * L * k^2 / 64)` dominated by row reduction |