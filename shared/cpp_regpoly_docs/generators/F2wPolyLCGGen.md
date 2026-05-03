# GenF2wPolyLCG

**C++ class:** `GenF2wPolyLCG`
**Name in code:** `"Generator in F_{2^w}"` / `"Polynomial LCG in F_{2^w}[z]/P(z)"`
**Legacy aliases:** `genf2w`

Polynomial Linear Congruential Generator in the extension field
GF(2^w). The state consists of r elements of GF(2^w), viewed as
coefficients of a polynomial in GF(2^w)[z] modulo a characteristic
polynomial. At each step the state polynomial is multiplied by z in
the quotient ring GF(2^w)[z] / P(z), where the coefficients of P(z)
are elements of GF(2^w).

## Mathematical recurrence

The state is a vector V[0], V[1], ..., V[r-1] of elements in GF(2^w).
At each step:

```
VP = V[r-1]
Right-shift the state by one position:
    V[r-1] = V[r-2], V[r-2] = V[r-3], ..., V[1] = V[0], V[0] = 0
If VP != 0, for each non-zero coefficient j = 0, ..., nbcoeff-1:
    V[nocoeff[j]] ^= multiply_GF2w(VP, coeff[j])
```

The multiplication `multiply_GF2w(a, b)` is performed in GF(2^w) with
the irreducible polynomial `modM` defining the field representation.
When `normal_basis` is true, multiplication uses a precomputed table
for the normal basis representation; otherwise it uses polynomial basis
arithmetic.

The characteristic polynomial of the recurrence has the form:

```
P(z) = z^r + coeff[0]*z^{nocoeff[0]} + coeff[1]*z^{nocoeff[1]} + ...
```

where each `coeff[j]` is an element of GF(2^w).

The generator has maximal period 2^k - 1 when P(z) is primitive over
GF(2^w).

## Parameters

| Name            | Type       | Role       | Rand type      | Description |
|-----------------|------------|------------|----------------|-------------|
| `w`             | `int`      | structural | --             | Word size in bits (degree of the field extension GF(2^w)) |
| `r`             | `int`      | structural | --             | Number of state words (degree of the characteristic polynomial) |
| `modM`          | `int`      | structural | --             | Irreducible polynomial over GF(2) of degree w, defining the field GF(2^w). Encoded as an integer bitmask (e.g., 0x169 for a degree-8 polynomial). |
| `normal_basis`  | `bool`     | structural | --             | If true, use normal basis representation for GF(2^w) arithmetic. Default: `false`. |
| `step`          | `int`      | structural | --             | Number of recurrence steps per call to `next()`. Default: `1`. |
| `nocoeff`       | `int_vec`  | structural | --             | Positions (exponents) of non-zero coefficients in the characteristic polynomial. Each value is in the range [0, r-1]. |
| `coeff`         | `uint_vec` | search     | `bitmask_vec`  | Coefficient values in GF(2^w), one per entry in `nocoeff`. Each is a w-bit mask. Randomized as independent w-bit bitmasks. |

### Coefficient encoding

The characteristic polynomial is fully determined by the `nocoeff` and
`coeff` vectors. Position `nocoeff[j]` gives the exponent of z, and
`coeff[j]` gives the corresponding GF(2^w) coefficient. The leading
term z^r is implicit and always present.

For example, with `r=3`, `nocoeff=[1, 0]`, `coeff=[0x5a, 0x3f]`, the
characteristic polynomial is z^3 + (0x5a)*z + (0x3f).

## State size (period exponent)

```
k = w * r
```

The period of the generator is 2^k - 1 when the characteristic
polynomial is primitive over GF(2^w) (equivalently, when the minimal
polynomial over GF(2) is primitive of degree w*r).

## Search examples

### Simple search -- find primitive generators with w=8, r=3

```yaml
search:
  family: GenF2wPolyLCG
  L: 8
  limit:
    max_tries: 10000

structural_params:
  w: 8
  r: 3
  modM: 0x169         # x^8 + x^6 + x^5 + x^3 + 1
  nocoeff: [1, 0]

fixed_params:
  coeff:               # randomized: two independent 8-bit GF(2^w) elements
```

### Fixed coefficients -- verify a known generator

```yaml
search:
  family: GenF2wPolyLCG
  L: 32
  limit:
    max_tries: 1

structural_params:
  w: 32
  r: 3
  modM: 0x100000285   # irreducible of degree 32
  nocoeff: [2, 0]

fixed_params:
  coeff: [0x5a3cc04e, 0x72f82b1a]   # fixed coefficients
```

### Normal basis with step > 1

```yaml
search:
  family: GenF2wPolyLCG
  L: 8
  output: genf2w_polylcg_nb_results.yaml
  limit:
    max_tries: 50000
    max_seconds: 30

structural_params:
  w: 8
  r: 5
  modM: 0x169
  normal_basis: true
  step: 2
  nocoeff: [3, 1, 0]

fixed_params:
  coeff:               # randomized: three 8-bit coefficients
```
