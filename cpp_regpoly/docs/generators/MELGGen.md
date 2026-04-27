# MELG

**C++ class:** `MELG`
**Name in code:** `"MELG"`

Maximally Equidistributed Long-period Generator of Harase & Kimoto
(2018).  A 64-bit F₂-linear generator family designed to improve on
MT19937-64 by reaching full multi-dimensional equidistribution for more
resolutions.

## State and recurrence

Let $N$ be the state length (in $w$-bit words, $w = 64$) and $r$ the
upper bit-mask width of the first state word.  With the $N+1$ words
$x_0, \ldots, x_N$ arranged in a circular buffer and $a$ a twist
constant, one step is

$$
x_{n+N} = x_{n+M}
        \oplus \bigl( (x_n^u \mathbin| x_{n+1}^l) \gg 1 \bigr)
        \oplus (a \cdot \mathrm{lsb}(x_{n+1}^l)),
$$

followed by two additional XORs with shifts of $\sigma_1$ and
$\sigma_2$ bits to complete the transition, and a tempering stage
defined by Harase & Kimoto.

## Parameters

| Name     | Kind        | Description |
|----------|-------------|-------------|
| `w`      | structural  | word width (default 64) |
| `N`      | structural  | number of state words |
| `r`      | structural  | upper-bit mask width |
| `M`      | randomised  | recurrence offset, range `[1, N-2]` |
| `sigma1` | randomised  | post-twist shift, range `[1, w-1]` |
| `sigma2` | randomised  | post-twist shift, range `[1, w-1]` |
| `a`      | randomised  | twist bitmask of width `w` |

## Published parametrisations

A 64-bit / 19937-state configuration is available in the library as
[`melg19937-64`](/library/melg19937-64) once its YAML is committed.

---

## References

- S. Harase and T. Kimoto. *Implementing 64-bit maximally equidistributed
  F₂-linear generators with Mersenne prime period.* ACM Trans. Math. Softw.
  **44** (2018), 30:1–30:11.
  [DOI](https://doi.org/10.1145/3159444)
