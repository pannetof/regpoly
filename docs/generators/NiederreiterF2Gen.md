# NiederreiterF2Gen — Niederreiter digital sequence in $\mathbb{F}_2$

**C++ class:** `NiederreiterF2Gen`
**Registered name:** `"NiederreiterF2Gen"`
**Aliases:** `niederreiter`, `Niederreiter`
**Base class:** `DigitalNet` → `Generator`

The Niederreiter digital sequence (Niederreiter 1988) in base 2.
For each coordinate $j$, the $m \times m$ generating matrix $C_j$ is built
from the $j$-th monic irreducible polynomial $p_j \in \mathbb{F}_2[x]$ via
the Laurent-series formula recorded in
[`tvalue-spec.md`](../theory/tvalue-spec.md). Irreducibles are enumerated
on the fly at construction (no checked-in data file) — milliseconds even
for $s_{\max} = 100$.

## Architectural fit

Same `DigitalNet` contract as `SobolNet`: the equidistribution kernel
runs unmodified, and the [t-value test](../theory/tvalue-spec.md) is the
natural companion.

## Generating-matrix formula (Niederreiter 1988)

Let $e_j = \deg p_j$. For each $r \ge 1$ and $0 \le k < e_j$, take the
Laurent expansion in $\mathbb{F}_2(\!(x^{-1})\!)$:

$$
\frac{x^k}{p_j(x)^r}
\;=\;
\sum_{\ell = 0}^{\infty} a^{(j)}(r, k, \ell)\, x^{-\ell - 1}.
$$

The entries of $C_j$ (1-based row index $i \in \{1, \ldots, m\}$,
0-based column index $\ell \in \{0, \ldots, m - 1\}$) are

$$
c_{j, i, \ell} \;=\; a^{(j)}(Q + 1,\, k,\, \ell),
\qquad \text{where } i - 1 = Q \cdot e_j + k,\;\; 0 \le k < e_j.
$$

i.e. Euclidean-divide $i - 1$ by $e_j$ to get $(Q, k)$, then read off the
$\ell$-th Laurent coefficient of $x^k / p_j(x)^{Q + 1}$.

## Ordering of irreducibles

Canonical (matches MinT/Salzburg, Dick–Pillichshammer 2010, QMCPy):

| $j$ | polynomial $p_j$         | binary value (leading-1 MSB) |
|---|----------------------------|------------------------------|
| 1 | $x$                        | $10$    |
| 2 | $x + 1$                    | $11$    |
| 3 | $x^2 + x + 1$              | $111$   |
| 4 | $x^3 + x + 1$              | $1011$  |
| 5 | $x^3 + x^2 + 1$            | $1101$  |
| 6 | $x^4 + x + 1$              | $10011$ |
| 7 | $x^4 + x^3 + 1$            | $11001$ |
| 8 | $x^4 + x^3 + x^2 + x + 1$  | $11111$ |
| … | …                          | …       |

Within each degree, lex order by binary value with the leading $1$ as
the MSB.

## Parameters

| Name    | Type  | Role       | Description |
|---------|-------|------------|-------------|
| `L`     | `int` | structural | Output bit width per coordinate; equals input width $m$. |
| `s_max` | `int` | structural | Number of generating matrices. Must satisfy $s_{\max} \ge m$. |

## v1 limits

- $m \le 30$. The Laurent recurrence uses `uint64_t` polynomial
  arithmetic and the highest polynomial degree touched is
  $m + e_{\max} - 1$; the 30-bit cap leaves generous headroom.
- $s_{\max} \le 100$ (irreducible enumeration runs in milliseconds even
  at this scale, but constructor cost grows roughly linearly).

## Usage

```python
import sys

from regpoly.core.generator import Generator
from regpoly.core.combination import Combination
from regpoly.analyses.tvalue_test import TValueTest

gen = Generator.create("NiederreiterF2Gen", L=10, s_max=10)
C = Combination.single(gen)

tvt = TValueTest(s_max=8, max_t_sum=sys.maxsize)
print(tvt.run(C).display())
```

## References

- H. Niederreiter, "Low-discrepancy and low-dispersion sequences,"
  *J. Number Theory* **30**, 51–70, 1988.
- P. Bratley, B. L. Fox & H. Niederreiter, "Implementation and tests of
  low-discrepancy sequences," *ACM Trans. Modeling and Computer
  Simulation* **2**(3), 195–213, 1992.
- P. Bratley, B. L. Fox & H. Niederreiter, "Algorithm 738: Programs to
  generate Niederreiter's low-discrepancy sequences," *ACM Trans. Math.
  Software* **20**(4), 494–495, 1994.
- J. Dick & F. Pillichshammer, *Digital Nets and Sequences: Discrepancy
  Theory and Quasi-Monte Carlo Integration*, Cambridge Univ. Press,
  2010. §8.1.
