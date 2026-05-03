# MersenneTwister

**C++ class:** `MersenneTwister`
**Name in code:** `"Mersenne Twister"`
**Legacy aliases:** `MT`

The Mersenne Twister algorithm (MT19937 family). A generalization of
the TGFSR that uses a partial word split controlled by a separation
parameter p, allowing the period to be a Mersenne prime 2^k - 1 where
k = w*r - p.

## Mathematical recurrence

The state consists of r words V[0], V[1], ..., V[r-1], each w bits
wide, accessed through a circular buffer with index i.

At each step:

```
upper_mask = bits w-1 down to p    (the top w-p bits)
lower_mask = bits p-1 down to 0    (the bottom p bits)

Y    = (V[i % r] & upper_mask) | (V[(i+1) % r] & lower_mask)
V[i] = V[(i+m) % r] XOR (Y >> 1) XOR (if Y is odd then a else 0)
i    = (i + 1) % r
```

The twist matrix A acts on Y in the same way as for TGFSR:

```
A * Y = (Y >> 1) XOR (if LSB(Y) == 1 then a else 0)
```

The partial-word split (controlled by p) allows the total state size
to be w*r - p bits rather than a full w*r, making it possible to
target Mersenne primes.

### Characteristic polynomial

```
P(t) = sum_{j=0}^{p-1} a_j * (t^{r-1} + t^{m-1})^{p-j-1} * (t^r + t^m)^{w-p}
     + sum_{j=p}^{w-1} a_j * (t^r + t^m)^{w-j-1}
     + (t^{r-1} + t^{m-1})^p * (t^r + t^m)^{w-p}
```

where a_j is bit j of the twist coefficient a.

## Parameters

| Name | Type  | Role       | Rand type | Rand args | Description |
|------|-------|------------|-----------|-----------|-------------|
| `w`  | `int` | structural | --        | --        | Word size in bits |
| `r`  | `int` | structural | --        | --        | Number of words in the state |
| `p`  | `int` | structural | --        | --        | Separation point (default: 0). Determines how many low bits come from the next word. |
| `m`  | `int` | search     | `range`   | `1,r-1`   | Feedback offset |
| `a`  | `int` | search     | `bitmask` | `w`       | Twist coefficient (w-bit mask) |

## State size (period exponent)

```
k = w * r - p
```

Classical examples:
- **MT19937:** w=32, r=624, p=31 => k = 32*624 - 31 = 19937
- **MT11213:** w=32, r=351, p=19 => k = 32*351 - 19 = 11213
- **MT44497:** w=32, r=1391, p=15 => k = 32*1391 - 15 = 44497

## Search examples

### Search for MT19937-class generators

```yaml
search:
  family: MersenneTwister
  L: 32
  limit:
    max_tries: 100000

structural_params:
  w: 32
  r: 624
  p: 31

fixed_params:
  m:     # randomized in [1, 623]
  a:     # randomized 32-bit mask
```

### Fixed offset, randomize twist coefficient (smaller period)

```yaml
search:
  family: MersenneTwister
  L: 32
  output: mt_small_results.yaml
  limit:
    max_tries: 50000
    max_seconds: 30

structural_params:
  w: 32
  r: 351
  p: 19

fixed_params:
  m: 175            # fixed near r/2
  a:                # randomized 32-bit mask
```

### No separation (p=0), equivalent to a TGFSR with MT recurrence

```yaml
search:
  family: MersenneTwister
  L: 32
  limit:
    max_tries: 10000

structural_params:
  w: 32
  r: 3
  p: 0

fixed_params:
  m:     # randomized in [1, 2]
  a:     # randomized 32-bit mask
```

---

## References

- M. Matsumoto and T. Nishimura. *Mersenne Twister: A 623-dimensionally
  equidistributed uniform pseudo-random number generator.*
  ACM Trans. Model. Comput. Simul. **8** (1998), 3–30.
  [DOI](https://doi.org/10.1145/272991.272995).  Published as
  [`mt19937`](/library/mt19937) in the library.
- T. Nishimura. *Tables of 64-bit Mersenne Twisters.*
  ACM Trans. Model. Comput. Simul. **10** (2000), 348–357.
  [DOI](https://doi.org/10.1145/369534.369540).
