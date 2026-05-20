# SobolNet — Sobol digital net in base 2

**C++ class:** `SobolNet`
**Registered name:** `"SobolNet"`
**Aliases:** `sobol`, `Sobol`
**Base class:** `DigitalNet` → `Generator`

The Sobol low-discrepancy net (Sobol 1967) in base 2, with direction
numbers from Joe & Kuo (2003). Each coordinate $j$ carries a fixed
$m \times m$ generating matrix $C_j$ over $\mathbb{F}_2$; the point set is

$$
\mathcal{P} \;=\; \bigl\{\, (C_1 \cdot \mathbf{i},\; C_2 \cdot \mathbf{i},\; \ldots,\;
C_{s_{\max}} \cdot \mathbf{i}) \;:\; \mathbf{i} \in \mathbb{F}_2^m \,\bigr\}
$$

with $\mathbf{i}$ the $m$-bit input index. The first coordinate is the
van der Corput / radical-inverse coordinate ($C_1 = I_m$); subsequent
coordinates are built via the Antonov–Saleev recurrence on direction
numbers driven by primitive polynomials over $\mathbb{F}_2$.

## Architectural fit

`SobolNet` inherits from `DigitalNet`, which itself inherits from
`Generator`. This means the existing **equidistribution machinery
runs on a SobolNet without modification** — the matricial kernel
sees only `Generator`'s `init` / `next` / `get_output` interface,
and `SobolNet` reinterprets those as "set input index → advance
coordinate counter → return $C_j \cdot \mathbf{i}$".

The [t-value test](../theory/tvalue-spec.md) is the natural
companion analysis, computing $t(s)$ for $s = 2, \ldots, s_{\max}$.

## Direction-number recurrence

For coordinate $j \ge 2$ with primitive polynomial of degree $e_j$
whose interior coefficients $a_1, \ldots, a_{e_j - 1}$ are packed in the
integer $a_j$ (bit $e_j - 1 - i$ of $a_j$ encodes $a_i$), the direction
numbers $m_{j,k}$ satisfy the Antonov–Saleev recurrence:

$$
m_{j,k} \;=\; 2 a_1\, m_{j, k - 1}
\;\oplus\; 4 a_2\, m_{j, k - 2}
\;\oplus\; \cdots
\;\oplus\; 2^{e_j - 1} a_{e_j - 1}\, m_{j, k - e_j + 1}
\;\oplus\; 2^{e_j}\, m_{j, k - e_j}
\;\oplus\; m_{j, k - e_j}, \qquad k > e_j,
$$

seeded by $m_{j, 1}, \ldots, m_{j, e_j}$ from the table. Column $k$
of $C_j$ (1-based, $k \in \{1, \ldots, m\}$) is the $m$-bit MSB-first
encoding of $m_{j, k} \ll (m - k)$.

## Parameters

| Name    | Type  | Role       | Description |
|---------|-------|------------|-------------|
| `L`     | `int` | structural | Output bit width per coordinate; equals input width $m$. |
| `s_max` | `int` | structural | Number of stored generating matrices. **Must satisfy $s_{\max} \ge m$** (the equidistribution kernel passes $\mathrm{indice\_max} = k_g = m$ to `GaussMatrix::prepare` and steps the generator $m - 1$ times). |

## v1 limits

- $m \le 64$ (`DigitalNet` base cap).
- $s_{\max} \le 8$ when using the embedded Joe–Kuo 2003 Table 1
  (v1 ships trusted direction numbers for coordinates $j = 2, \ldots, 8$;
  coordinate $j = 1$ is the built-in identity). Constructing
  `SobolNet(m, s_max)` with $s_{\max} > 8$ raises
  `std::invalid_argument`. Extending to the full 32-dimensional
  `new-joe-kuo-6.21201` table is a follow-up requiring user-authorized
  fetch of that file.

## Usage

```python
import sys

from regpoly.core.generator import Generator
from regpoly.core.combination import Combination
from regpoly.analyses.equidistribution_test import EquidistributionTest
from regpoly.analyses.tvalue_test import TValueTest

# Build a SobolNet with m = 4 bits / coordinate, s_max = 4 coordinates.
gen = Generator.create("SobolNet", L=4, s_max=4)
C = Combination.single(gen)

# Existing equidistribution test runs unchanged:
eqt = EquidistributionTest(L=4, delta=[sys.maxsize] * 5, mse=sys.maxsize)
print(eqt.run(C).display_table(C, by="l"))

# t-value test:
tvt = TValueTest(s_max=4, max_t_sum=sys.maxsize)
print(tvt.run(C).display())
```

## YAML

```yaml
generator:
  family: SobolNet
  L: 4
  s_max: 4
tests:
  - type: equidistribution
    method: matricial
  - type: tvalue
    s_max: 4
    max_t_sum: 99
```

## References

- I. M. Sobol, "On the distribution of points in a cube and the
  approximate evaluation of integrals," *USSR Comput. Math. and Math.
  Phys.* **7**(4), 86–112, 1967.
- I. A. Antonov & V. M. Saleev, "An economic method of computing
  L.P.τ-sequences," *USSR Comput. Math. and Math. Phys.* **19**(1),
  252–256, 1979.
- S. Joe & F. Y. Kuo, "Remark on Algorithm 659: Implementing Sobol's
  quasirandom sequence generator," *ACM Trans. Math. Software*
  **29**(1), 49–57, 2003.
- S. Joe & F. Y. Kuo, "Constructing Sobol sequences with better
  two-dimensional projections," *SIAM J. Sci. Comput.* **30**(5),
  2635–2654, 2008.
  ([direction-number tables](https://web.maths.unsw.edu.au/~fkuo/sobol/))
