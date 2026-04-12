# AC1DGen

**C++ class:** `AC1DGen`
**Name in code:** `"AC-1D"`
**Legacy aliases:** `AC1D`

Additive Congruential 1-Dimensional generator. A pure binary linear
transformation generator whose state evolves by multiplication with an
n x n binary matrix over GF(2). The transition matrix defines the
recurrence completely.

## Mathematical recurrence

The state is an n-bit vector s[0], s[1], ..., s[n-1]. At each step,
the new state is computed by multiplying the current state by the
binary transition matrix A:

```
For i = 0, ..., n-1:
    new_state[i] = XOR over all j where A[i][j] = 1 of state[j]
state = new_state
```

In matrix notation:

```
s(t+1) = A * s(t)     (over GF(2))
```

where A is the n x n binary transition matrix and all arithmetic is
modulo 2 (XOR for addition, AND for multiplication).

The generator has maximal period 2^n - 1 when the characteristic
polynomial of A is primitive over GF(2).

## Parameters

| Name     | Type      | Role       | Rand type | Description |
|----------|-----------|------------|-----------|-------------|
| `n`      | `int`     | structural | --        | State dimension (number of bits). The transition matrix is n x n. |
| `matrix` | `int_vec` | search     | `none`    | Flattened n x n binary matrix in row-major order. Entry matrix[i*n + j] is 1 if position (i,j) of the transition matrix is set, 0 otherwise. Total length: n*n. |

### Matrix encoding

The matrix is specified as a flat vector of n*n integers (0 or 1) in
row-major order. For an n=3 matrix:

```
A = | a00 a01 a02 |
    | a10 a11 a12 |
    | a20 a21 a22 |
```

The flat encoding is: `[a00, a01, a02, a10, a11, a12, a20, a21, a22]`.

**Note:** Random generation is not supported for the matrix parameter
(`rand_type=none`). The transition matrix must be provided explicitly.

## State size (period exponent)

```
k = n
```

The period of the generator is 2^n - 1 when the characteristic
polynomial of the transition matrix A is primitive over GF(2).

## Search examples

Since random generation is not supported for this family, the examples
below show verification of known generators.

### Verify a small n=5 companion matrix

```yaml
search:
  family: AC1DGen
  L: 5
  limit:
    max_tries: 1

structural:
  n: 5

search_params:
  matrix: [
    0, 1, 0, 0, 0,
    0, 0, 1, 0, 0,
    0, 0, 0, 1, 0,
    0, 0, 0, 0, 1,
    1, 0, 1, 0, 0
  ]
```

This encodes the companion matrix of the primitive polynomial
z^5 + z^2 + 1, which gives a maximal-length sequence with period
2^5 - 1 = 31.

### Verify an 8x8 companion matrix

```yaml
search:
  family: AC1DGen
  L: 8
  limit:
    max_tries: 1

structural:
  n: 8

search_params:
  matrix: [
    0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 1,
    1, 0, 0, 1, 1, 1, 0, 1
  ]
```

This encodes the companion matrix of the primitive polynomial
z^8 + z^4 + z^3 + z^2 + 1 (0x11D), which gives a maximal-length
sequence with period 2^8 - 1 = 255.
