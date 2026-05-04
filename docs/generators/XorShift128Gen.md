# XorShift128 — Marsaglia 128-bit Xor-shift

**C++ class:** `XorShift128Gen`
**Name in code:** `"XorShift128Gen"` (alias: `"XorShift128"`)

`XorShift128Gen` is the 128-bit-state Xor-shift generator of
Marsaglia (2003): four 32-bit words rotated through a triple-shift
recurrence that produces 32-bit output. Two recurrence shapes are
supported via the `pattern` field, matching MTToolBox's
`samples/XORSHIFT/xorshift-2.cpp` (canonical Marsaglia 2003 form,
`pattern=0`) and `samples/XORSHIFT/xorshift-5.cpp` (search-result
form, `pattern=1`).

## Mathematical recurrence

The state is four 32-bit words `(x, y, z, w)`; the seed convention
used by MTToolBox `calc_equidist` is `w = ~v; x = y = z = v`.

For `pattern == 0` (xorshift-2 / Marsaglia 2003 canonical):

```
t = x ^ (x << a)
x = y;  y = z;  z = w
w = (w ^ (w >> c)) ^ (t ^ (t >> b))
out = w
```

For `pattern == 1` (xorshift-5 / search-result form):

```
t = x ^ (x << b)
t ^= t >> c
x = y;  y = z;  z = w
w = (w ^ (w << a)) ^ t
out = w
```

The shift triple `(a, b, c)` is the search axis.

### Characteristic polynomial

For shift triples in Marsaglia's original list of 81 (and the
xorshift-5-form analogues from `xorshift-5.cpp`), $\chi_f$ is the
primitive polynomial of degree 128 over $\mathrm{GF}(2)$ giving the
maximum-length period $2^{128} - 1$. regpoly recovers the minimal
polynomial at runtime via Berlekamp–Massey on the bit-stream
produced by `next()`.

## Parameters

| Name      | Type  | Role       | Rand type | Rand args | Description |
|-----------|-------|------------|-----------|-----------|-------------|
| `a`       | `int` | search     | `range`   | `1,31`    | First shift count of the triple. |
| `b`       | `int` | search     | `range`   | `1,31`    | Second shift count of the triple. |
| `c`       | `int` | search     | `range`   | `1,31`    | Third shift count of the triple. |
| `pattern` | `int` | structural | --        | --        | Recurrence shape: `0` = xorshift-2, `1` = xorshift-5. |

Roles:

- **structural** — fixed by the user; defines the shape of $A$.
- **search** — randomized by the search loop. The shift triple
  `(a, b, c)` is the canonical search axis; `pattern` selects
  between the two recurrence shapes.

## State size (period exponent)

```
k = 128
```

There is one canonical XorShift128 size; the family does not vary
$k$ across published parametrisations, only the shift triple
`(a, b, c)` and the recurrence shape `pattern`.

## Search examples

### Search for xorshift-2-form triples

```yaml
search:
  family: XorShift128Gen
  L: 32
  limit:
    max_tries: 10000

structural_params:
  pattern: 0

fixed_params:
  a:        # randomized in [1, 31]
  b:        # randomized in [1, 31]
  c:        # randomized in [1, 31]
```

### Verify a published Marsaglia (5, 14, 1) triple

```yaml
search:
  family: XorShift128Gen
  L: 32
  limit:
    max_tries: 1

structural_params:
  pattern: 0

fixed_params:
  a: 5
  b: 14
  c: 1
```

### Verify an xorshift-5-form (20, 11, 7) triple

```yaml
search:
  family: XorShift128Gen
  L: 32
  limit:
    max_tries: 1

structural_params:
  pattern: 1

fixed_params:
  a: 20
  b: 11
  c: 7
```

## Implementation notes

- **State word ordering.** The state is stored MSB-first across the
  four 32-bit slots so that `state[0]` ↔ `w`, `state[1]` ↔ `z`,
  `state[2]` ↔ `y`, `state[3]` ↔ `x`. The `load_state` /
  `store_state` helpers absorb this convention; the recurrence
  pseudo-code above uses the natural `(x, y, z, w)` order.
- **Pattern 1 is not Marsaglia 2003.** Triples found by
  `xorshift-5.cpp` are valid full-period generators but use a
  different recurrence shape; the two patterns are not
  interchangeable parameter spaces.

---

## References

- G. Marsaglia. *Xorshift RNGs.* Journal of Statistical Software
  **8(14)** (2003), 1–6.
  [DOI](https://doi.org/10.18637/jss.v008.i14).
  Published as
  [`marsaglia-2003-xorshift`](/library/marsaglia-2003-xorshift)
  in the library.
