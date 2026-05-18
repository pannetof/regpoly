# F₂-linear generators

This page is the umbrella reference for the family of pseudo-random number generators REGPOLY targets: **F₂-linear** generators, where the state evolves under a linear recurrence over the two-element field $\mathbb{F}_2 = \mathrm{GF}(2)$. Almost every modern high-quality bit-source for simulation — Mersenne Twister, WELL, MELG, SFMT, Tausworthe, xorshift — falls in this family.

## State, recurrence, and output

A generator with **state size** $k$ keeps a vector $\mathbf{x}_n \in \mathbb{F}_2^k$ at step $n$. Its evolution is governed by a linear map $A \in \mathbb{F}_2^{k \times k}$:

$$
\mathbf{x}_{n+1} = A \,\mathbf{x}_n .
$$

The output at step $n$ is read off the state through a second linear map $B \in \mathbb{F}_2^{w \times k}$ that selects $w$ output bits:

$$
\mathbf{y}_n = B \,\mathbf{x}_n .
$$

Finally, the integer or floating-point sample the user actually consumes is the binary expansion of $\mathbf{y}_n$, optionally postprocessed by a **tempering** linear map $T \in \mathbb{F}_2^{w \times w}$ (see [Tempering](#tempering-as-a-linear-map) below).

The matrices $A$, $B$, and $T$ together fully specify the generator. In implementation, $A$ is rarely written down explicitly — the family-specific recurrence (a TGFSR twist, a MELG lagged feedback, a Tausworthe LFSR) IS the description of $A$, encoded in word-sized operations that are equivalent to one matrix-vector product over $\mathbb{F}_2$.

## Characteristic polynomial and full period

The **characteristic polynomial** of the generator is

$$
\chi(t) = \det(tI - A) \in \mathbb{F}_2[t] ,
$$

a polynomial of degree $k$ over $\mathbb{F}_2$. Its factorization controls the cycle structure of the state space:

- **Full period.** If $\chi(t)$ is **primitive** over $\mathbb{F}_2$ (irreducible AND its roots generate $\mathbb{F}_{2^k}^*$), then $A$ has order $2^k - 1$. Every non-zero state lies on a single cycle of length $2^k - 1$. This is what families like MT19937 are designed for: the state size $k$ is chosen so $2^k - 1$ is a Mersenne prime, which simplifies primitivity checks.
- **Reducible $\chi(t)$.** If $\chi(t)$ factors, the state space splits into invariant subspaces, and cycle lengths divide the orders of the irreducible factors. REGPOLY's [equidistribution machinery](equidistribution-spec.md) is designed to handle this case: it does not assume full period.

Primitivity (and full period) is checked in REGPOLY by the `is_full_period` driver in `regpoly_cpp::algebra::primitivity`, backed by an embedded Cunningham factor table.

## Equidistribution

The defining quality metric for an F₂-linear generator is **equidistribution**: how uniformly its output points fill the unit hypercube $[0, 1)^L$.

For a single resolution $\ell \in \{1, \ldots, w\}$, the **dimension of equidistribution** $d(\ell)$ is the largest $L$ such that the set of all $L$-tuples of consecutive outputs, truncated to $\ell$ leading bits each, covers every $\ell L$-bit pattern exactly $2^{k - \ell L}$ times. Equivalently, the matrix that maps the state to those $\ell L$ output bits has rank $\ell L$.

The **dimension gap** at resolution $\ell$ is

$$
\delta(\ell) = \left\lfloor k / \ell \right\rfloor - d(\ell) ,
$$

with $\sum_\ell \delta(\ell)$ used as the scalar score for ranking parameter sets in a search. A generator is **maximally equidistributed (ME)** if every gap is zero. See [Equidistribution spec](equidistribution-spec.md) for the matricial method REGPOLY uses to compute $d(\ell)$ on generators that may not be full-period.

Two related metrics are also computed:

- **Antithetic structure** — see [Antithetic check](antithetic-check.md).
- **$t$-tuples test** — Couture & L'Ecuyer's $\Delta(t_1, \ldots, t_d)$ over the lattice [`couture-lecuyer`](lattice_method_couture_lecuyer.md), and the Harase variant in [`lattice_method_harase`](lattice_method_harase.md).

## Combined generators

REGPOLY's main subject is the combined generator

$$
\mathbf{y}_n = \mathbf{y}^{(1)}_n \oplus \mathbf{y}^{(2)}_n \oplus \cdots \oplus \mathbf{y}^{(J)}_n ,
$$

where each component $j$ is itself an F₂-linear generator of state size $k_j$, and the outputs are XORed bitwise. The combined state has size $k_g = \sum_j k_j$ and characteristic polynomial $\chi_g(t) = \prod_j \chi_j(t)$ (provided the $\chi_j$ are pairwise coprime).

Combined generators improve equidistribution by widening the state space without paying the per-step cost of running a single huge component, and by letting the combination "average out" the equidistribution defects of the individual components. The `Combination` C++ class builds these and is the live object every search loop iterates.

## Tempering as a linear map

Many families wrap a **tempering** linear bijection $T : \mathbb{F}_2^w \to \mathbb{F}_2^w$ around the raw output:

$$
\mathbf{z}_n = T \,\mathbf{y}_n .
$$

Common operations in $T$:

- shift-XOR: $T_s(\mathbf{y}) = \mathbf{y} \oplus (\mathbf{y} \gg s)$ — a one-bit-shifted XOR with itself,
- shift-and-mask-XOR: $T_{s,m}(\mathbf{y}) = \mathbf{y} \oplus ((\mathbf{y} \ll s) \mathbin{\&} m)$ — used by MT.

A tempering chain is a composition of such steps. Crucially, $T$ is invertible (it has a triangular structure on $\mathbb{F}_2$), so it does NOT change the cycle length or the characteristic polynomial — but it dramatically improves the equidistribution of low-order bits, which an LFSR otherwise produces poorly. REGPOLY's [tempering optimizer](tempering_optimization.md) searches for $T$ that minimizes $\sum_\ell \delta(\ell)$ for a fixed component $A, B$.

## See also

- Per-family details: [`generators/`](../generators/index.md) — one page per family (MT, SFMT, MELG, WELL, Tausworthe, …).
- Algorithmic specs:
    - [Equidistribution computation](equidistribution-spec.md) — the design of record for matricial equidistribution on possibly-non-primitive generators.
    - [Antithetic check](antithetic-check.md) — local antitheticity of a linear point set.
    - [Couture–L'Ecuyer lattice method](lattice_method_couture_lecuyer.md) and the [Harase variant](lattice_method_harase.md).
    - [Tempering optimization](tempering_optimization.md).
    - [Search YAML format](search_format.md).
