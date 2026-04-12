# TGFSRGen

**C++ class:** `TGFSRGen`
**Name in code:** `"TGFSR"`
**Legacy aliases:** `tgfsr`

Twisted Generalized Feedback Shift Register. The state consists of r
words of w bits each. The feedback uses a "twist" matrix A that
conditionally XORs a coefficient based on the least-significant bit of
the shifted word.

## Mathematical recurrence

The state at time i is a sequence of r words V[0], V[1], ..., V[r-1],
each w bits wide. The recurrence updates the oldest word:

```
Y       = V[i] >> 1
V_new   = V[i-r+m] XOR Y XOR (if LSB(V[i]) == 1 then a else 0)
```

The state register then shifts: V[i] is replaced by V_new, and the
index advances.

The twist matrix A acts on a w-bit word as:

```
        ( 0  1  0 ... 0 )
        ( 0  0  1 ... 0 )
A * x = ( .  .  . ... . ) * x  =  (x >> 1) XOR (if x is odd then a else 0)
        ( 0  0  0 ... 1 )
        (a_{w-1} a_{w-2} ... a_0)
```

where a = (a_{w-1}, a_{w-2}, ..., a_0) is the twist coefficient.

### Characteristic polynomial

The characteristic polynomial of the generator over GF(2) is:

```
P(t) = (t^r + t^m)^w + sum_{j: bit j of a is set} (t^r + t^m)^j
```

The generator has maximal period 2^k - 1 when P(t) is primitive.

## Parameters

| Name | Type  | Role       | Rand type | Rand args | Description |
|------|-------|------------|-----------|-----------|-------------|
| `w`  | `int` | structural | --        | --        | Word size in bits |
| `r`  | `int` | structural | --        | --        | Number of words in the state |
| `m`  | `int` | search     | `range`   | `1,r-1`   | Feedback offset (distance in the state array) |
| `a`  | `int` | search     | `bitmask` | `w`       | Twist coefficient (w-bit mask) |

## State size (period exponent)

```
k = w * r
```

The period is 2^(w*r) - 1 when the characteristic polynomial is
primitive.

## Search examples

### Simple search for a small TGFSR

```yaml
search:
  family: TGFSRGen
  L: 32
  limit:
    max_tries: 10000

structural:
  w: 32
  r: 3

search_params:
  m:     # randomized in [1, 2]
  a:     # randomized 32-bit mask
```

### Fixed feedback offset, randomize twist coefficient

```yaml
search:
  family: TGFSRGen
  L: 32
  output: tgfsr_w32r7_results.yaml
  limit:
    max_tries: 50000
    max_seconds: 30

structural:
  w: 32
  r: 7

search_params:
  m: 4               # fixed offset
  a:                  # randomized 32-bit mask
```

### Larger state space

```yaml
search:
  family: TGFSRGen
  L: 32
  limit:
    max_tries: 200000

structural:
  w: 32
  r: 25

search_params:
  m:     # randomized in [1, 24]
  a:     # randomized 32-bit mask
```
