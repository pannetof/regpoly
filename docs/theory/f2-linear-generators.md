# F₂-linear generators

This page is the umbrella reference for the family of pseudo-random number generators REGPOLY targets: **F₂-linear** generators, where the state evolves under a linear recurrence over the two-element field $\mathbb{F}_2$ (also written $\mathrm{GF}(2)$ in the older literature). Almost every modern high-quality bit-source for simulation — Mersenne Twister, WELL, MELG, SFMT, Tausworthe, xorshift — falls in this family.

## State, recurrence, and output

A generator with **state size** $k$ keeps a vector $\mathbf{x}_n \in \mathbb{F}_2^k$ at step $n$. Its evolution is governed by a linear **state-transition map**

$$
f : \mathbb{F}_2^k \to \mathbb{F}_2^k, \qquad f(\mathbf{x}) := A \mathbf{x},
$$

where $A \in \mathbb{F}_2^{k \times k}$ is the **state-transition matrix**. One step of the generator is one application of $f$:

$$
\mathbf{x}_{n+1} = f(\mathbf{x}_n) = A \,\mathbf{x}_n .
$$

The output at step $n$ is read off the state through a second linear map $B \in \mathbb{F}_2^{L \times k}$ that selects $L$ output bits:

$$
\mathbf{y}_n = B \,\mathbf{x}_n .
$$

Finally, the integer or floating-point sample the user actually consumes is the binary expansion of $\mathbf{y}_n$, optionally postprocessed by a **tempering** linear map $T \in \mathbb{F}_2^{L \times L}$ (see [Tempering](#tempering-as-a-linear-map) below).

Throughout this and every other theory page, **bit-vector quantities are bold** ($\mathbf{x}, \mathbf{y}, \mathbf{u}, \mathbf{0}, \mathbf{1}$) and scalars / matrices are plain italic ($k, \ell, t, A, B, T$). See [Notation](notation.md) for the full convention.

The matrices $A$, $B$, and $T$ together fully specify the generator. In implementation, $A$ is rarely written down explicitly — the family-specific recurrence (a TGFSR twist, a MELG lagged feedback, a Tausworthe LFSR) IS the description of $A$, encoded in word-sized operations that are equivalent to one matrix-vector product over $\mathbb{F}_2$.

## Characteristic polynomial and full period

The **characteristic polynomial** of the state-transition map $f$ is

$$
\chi_f(t) = \det(t I_k - A) \in \mathbb{F}_2[t] ,
$$

where $I_k$ is the $k \times k$ identity matrix and $t$ is here the formal polynomial variable (distinct from the tuple-dimension scalar $t$ used in the equidistribution section below — the two never appear in the same expression). $\chi_f(t)$ is a polynomial of degree $k$ over $\mathbb{F}_2$, and its factorization controls the cycle structure of the state space:

- **Full period.** If $\chi_f(t)$ is **primitive** over $\mathbb{F}_2$ (irreducible AND its roots generate the multiplicative group $\mathbb{F}_{2^k}^* = \mathbb{F}_{2^k} \setminus \{0\}$ of the extension field $\mathbb{F}_{2^k}$), then $A$ has order $2^k - 1$. Every non-zero state lies on a single cycle of length $2^k - 1$. This is what families like MT19937 are designed for: the state size $k$ is chosen so $2^k - 1$ is a Mersenne prime, which simplifies primitivity checks.
- **Reducible $\chi_f(t)$.** If $\chi_f(t)$ factors, the state space splits into invariant subspaces, and cycle lengths divide the orders of the irreducible factors. REGPOLY's [equidistribution machinery](equidistribution-spec.md) is designed to handle this case: it does not assume full period.

Primitivity (and full period) is checked in REGPOLY by the `is_full_period` driver in `regpoly_cpp::algebra::primitivity`, backed by an embedded Cunningham factor table.

## Equidistribution

The defining quality metric for an F₂-linear generator is **equidistribution**: how uniformly its output points fill the unit hypercube $[0, 1)^t$, where $t$ is the **tuple dimension** (a positive integer).

For a single **bit resolution** $\ell \in \{1, \ldots, L\}$, the **dimension of equidistribution** $d(\ell)$ is the largest $t$ such that the set of all $t$-tuples of consecutive outputs, truncated to $\ell$ leading bits each, covers every $\ell t$-bit pattern exactly $2^{k - \ell t}$ times. Equivalently, the matrix that maps the state to those $\ell t$ output bits has rank $\ell t$.

The **dimension gap** at resolution $\ell$ is

$$
\delta(\ell) = \left\lfloor k / \ell \right\rfloor - d(\ell) ,
$$

and the **sum of dimension gaps**

$$
\text{SE} = \sum_{\ell=1}^{L} \delta(\ell)
$$

is the scalar score REGPOLY's search loops minimise when ranking parameter sets. A generator is **maximally equidistributed (ME)** if $\delta(\ell) = 0$ for every $\ell \in \{1, \ldots, L\}$. See [Equidistribution spec](equidistribution-spec.md) for the matricial method REGPOLY uses to compute $d(\ell)$ on generators that may not be full-period.

Two related metrics are also computed:

- **Antithetic structure** — see [Antithetic check](antithetic-check.md).
- **$t$-tuples test** — Couture & L'Ecuyer's $\Delta(t_1, \ldots, t_r)$ figure of merit (where $r$ is the number of stacked tuple dimensions $t_1, \ldots, t_r$) over the lattice [`couture-lecuyer`](lattice_method_couture_lecuyer.md), and the Harase variant in [`lattice_method_harase`](lattice_method_harase.md).

## Combined generators

REGPOLY's main subject is the combined generator, formed from $J \ge 2$ component generators (the **number of components**, a positive integer). Let $j \in \{1, \ldots, J\}$ index a component, with component-$j$ state size $k_j$, raw output $\mathbf{y}^{(j)}_n \in \mathbb{F}_2^L$, and characteristic polynomial $\chi_j(t) \in \mathbb{F}_2[t]$. The combined raw output is

$$
\mathbf{y}_n = \mathbf{y}^{(1)}_n \oplus \mathbf{y}^{(2)}_n \oplus \cdots \oplus \mathbf{y}^{(J)}_n ,
$$

i.e. the components' outputs XORed bitwise. The **combined state size** is $k_g = \sum_{j=1}^{J} k_j$ and the **combined characteristic polynomial** is $\chi_g(t) = \prod_{j=1}^{J} \chi_j(t)$ (provided the $\chi_j$ are pairwise coprime). Where the umbrella symbol $\chi_f$ from the Notation reference applies, the combined generator's $\chi_f$ is exactly this $\chi_g$.

Combined generators improve equidistribution by widening the state space without paying the per-step cost of running a single huge component, and by letting the combination "average out" the equidistribution defects of the individual components. The `Combination` C++ class builds these and is the live object every search loop iterates.

## Tempering as a linear map

Many families wrap a **tempering** linear bijection $T : \mathbb{F}_2^L \to \mathbb{F}_2^L$ around the raw output. The tempered output (canonical symbol from the [Notation](notation.md) reference) is

$$
\mathbf{u}_n = T \,\mathbf{y}_n .
$$

Common operations in $T$, for a non-negative integer **shift amount** $s$ and a bitmask $\mathbf{m} \in \mathbb{F}_2^L$:

- shift-XOR: $T_s(\mathbf{y}) = \mathbf{y} \oplus (\mathbf{y} \gg s)$ — XOR with the bit-vector $\mathbf{y}$ right-shifted by $s$ bits,
- shift-and-mask-XOR: $T_{s,\mathbf{m}}(\mathbf{y}) = \mathbf{y} \oplus \bigl((\mathbf{y} \ll s) \mathbin{\&} \mathbf{m}\bigr)$ — used by MT.

A tempering chain is a composition of such steps. Crucially, $T$ is invertible (it has a triangular structure on $\mathbb{F}_2$), so it does NOT change the cycle length or the characteristic polynomial $\chi_f$ — but it dramatically improves the equidistribution of low-order bits, which an LFSR otherwise produces poorly. REGPOLY's [tempering optimizer](tempering_optimization.md) searches for $T$ that minimizes the sum of dimension gaps $\text{SE}$ for a fixed component $A, B$.

## See also

- Per-family details: [`generators/`](../generators/index.md) — one page per family (MT, SFMT, MELG, WELL, Tausworthe, …).
- Algorithmic specs:
    - [Equidistribution computation](equidistribution-spec.md) — the design of record for matricial equidistribution on possibly-non-primitive generators.
    - [Antithetic check](antithetic-check.md) — local antitheticity of a linear point set.
    - [Couture–L'Ecuyer lattice method](lattice_method_couture_lecuyer.md) and the [Harase variant](lattice_method_harase.md).
    - [Tempering optimization](tempering_optimization.md).
    - [Search YAML format](search_format.md).
