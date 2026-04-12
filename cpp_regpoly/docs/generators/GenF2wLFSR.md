# GenF2wLFSR

**C++ class:** `GenF2wLFSR`
**Name in code:** `"Generator in F_{2^w}"` / `"LFSR in F_{2^w}"`
**Legacy aliases:** `genf2w` (with type=lfsr)

Linear Feedback Shift Register in the extension field GF(2^w). The
state consists of r elements of GF(2^w). Unlike the Polynomial LCG
variant (GenF2wPolyLCG), which right-shifts and feeds back from the
last position, this LFSR variant left-shifts the state and computes
the new element as a linear combination of the current state words.

## Mathematical recurrence

The state is a vector V[0], V[1], ..., V[r-1] of elements in GF(2^w).
At each step:

```
res = 0
For each non-zero coefficient j = 0, ..., nbcoeff-1:
    res ^= multiply_GF2w(V[nocoeff[j]], coeff[j])
Left-shift the state by one position:
    V[0] = V[1], V[1] = V[2], ..., V[r-2] = V[r-1]
V[r-1] = res
```

The multiplication `multiply_GF2w(a, b)` is performed in GF(2^w) with
the irreducible polynomial `modM` defining the field representation.
When `normal_basis` is true, multiplication uses a precomputed table
for the normal basis representation; otherwise it uses polynomial basis
arithmetic.

The feedback polynomial (characteristic polynomial) has the form:

```
P(z) = z^r + coeff[0]*z^{r - nocoeff[0]} + ... + coeff[t-1]*z^{r - nocoeff[t-1]}
```

where each `coeff[j]` is an element of GF(2^w) and the `nocoeff[j]`
identify which state positions are tapped.

The generator has maximal period 2^k - 1 when P(z) is primitive over
GF(2^w).

**Note:** During initialization and after each `next()` call, the LFSR
extends its output by computing additional words beyond the r-word
state, using the same recurrence. This allows L (output resolution) to
exceed k when k is small.

## Parameters

| Name            | Type       | Role       | Rand type      | Description |
|-----------------|------------|------------|----------------|-------------|
| `w`             | `int`      | structural | --             | Word size in bits (degree of the field extension GF(2^w)) |
| `r`             | `int`      | structural | --             | Number of state words (degree of the characteristic polynomial) |
| `modM`          | `int`      | structural | --             | Irreducible polynomial over GF(2) of degree w, defining the field GF(2^w). Encoded as an integer bitmask. |
| `normal_basis`  | `bool`     | structural | --             | If true, use normal basis representation for GF(2^w) arithmetic. Default: `false`. |
| `step`          | `int`      | structural | --             | Number of recurrence steps per call to `next()`. Default: `1`. |
| `nocoeff`       | `int_vec`  | structural | --             | Positions (indices into the state vector) of non-zero feedback taps. Each value is in the range [0, r-1]. |
| `coeff`         | `uint_vec` | search     | `bitmask_vec`  | Coefficient values in GF(2^w), one per entry in `nocoeff`. Each is a w-bit mask. Randomized as independent w-bit bitmasks. |

### Coefficient encoding

The feedback is specified by `nocoeff` and `coeff` vectors. The entry
`nocoeff[j]` identifies which state element V[nocoeff[j]] is tapped,
and `coeff[j]` gives the GF(2^w) multiplier for that tap.

For example, with `r=5`, `nocoeff=[0, 3]`, `coeff=[0x1b, 0x47]`, the
new state element is:

```
V_new[4] = multiply(V[0], 0x1b) XOR multiply(V[3], 0x47)
```

## State size (period exponent)

```
k = w * r
```

The period of the generator is 2^k - 1 when the characteristic
polynomial is primitive over GF(2^w).

## Search examples

### Simple search -- find primitive LFSR generators with w=8, r=3

```yaml
search:
  family: GenF2wLFSR
  L: 8
  limit:
    max_tries: 10000

structural:
  w: 8
  r: 3
  modM: 0x169         # x^8 + x^6 + x^5 + x^3 + 1
  nocoeff: [0, 2]

search_params:
  coeff:               # randomized: two independent 8-bit GF(2^w) elements
```

### Fixed coefficients -- verify a known generator

```yaml
search:
  family: GenF2wLFSR
  L: 32
  limit:
    max_tries: 1

structural:
  w: 32
  r: 3
  modM: 0x100000285   # irreducible of degree 32
  nocoeff: [0, 1]

search_params:
  coeff: [0xa34f09b2, 0x7c51e3d4]   # fixed coefficients
```

### LFSR with step > 1 and output file

```yaml
search:
  family: GenF2wLFSR
  L: 8
  output: genf2w_lfsr_results.yaml
  limit:
    max_tries: 50000
    max_seconds: 30

structural:
  w: 8
  r: 7
  modM: 0x169
  step: 3
  nocoeff: [0, 2, 5]

search_params:
  coeff:               # randomized: three 8-bit coefficients
```
