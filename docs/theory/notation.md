# Notation

Canonical mathematical symbols used across REGPOLY's theory pages and
docstrings. Every page reads with the conventions below unless it
explicitly redefines a symbol locally (the BM and MELG sections both
do this for `L` — see [Page-local overrides](#page-local-overrides)).

When you read a new docstring or theory section, this is the table to
check first.

## Typographic conventions

- **Bit vectors are bold.** Every element of $\mathbb{F}_2^k$ or
  $\mathbb{F}_2^L$ is written with `\mathbf{}` — e.g. the state
  $\mathbf{x}$, the output word $\mathbf{y}$, the tempered output
  $\mathbf{u}$, the all-ones vector $\mathbf{1}$, the zero vector
  $\mathbf{0}$. Scalars stay in plain italic ($k$, $\ell$, $t$).
- **Matrices are plain capitals** in italic: $A$, $B$, $T$.
- **Polynomials and field elements** are plain italic ($\chi_f$,
  $\varphi$, $p$).
- **Code identifiers** use a monospace font, never math. A function
  parameter named `L` in C++ is `L` in code-font even though the
  underlying math symbol is the same $L$.

## Field and vector spaces

| Symbol | Meaning |
|--------|---------|
| $\mathbb{F}_2$ | The two-element field $\{0, 1\}$ with addition $=$ XOR and multiplication $=$ AND. |
| $\mathbb{F}_2^k$ | The $k$-dimensional vector space over $\mathbb{F}_2$ — the **state space** of a generator with state size $k$. |
| $\mathbb{F}_2^{m \times n}$ | The set of $m \times n$ matrices over $\mathbb{F}_2$. |

The historical alternative `\mathrm{GF}(2)` is the same field; the
docs use `\mathbb{F}_2`.

## Generator structure

| Symbol | Type | Meaning |
|--------|------|---------|
| $k$ | positive integer | **State size** in bits — the dimension of the F₂-linear state vector. |
| $L$ | positive integer | **Output bit width** — number of bits emitted per `next()` call. Matches the Python `L=32` constructor argument and the C++ `int L` parameter throughout. |
| $\mathbf{x}_n \in \mathbb{F}_2^k$ | column vector | Generator **state** at step $n$. |
| $\mathbf{y}_n \in \mathbb{F}_2^L$ | column vector | Raw **output word** at step $n$. |
| $\mathbf{u}_n \in \mathbb{F}_2^L$ | column vector | **Tempered output** at step $n$, $\mathbf{u}_n = T \cdot \mathbf{y}_n$. |
| $A \in \mathbb{F}_2^{k \times k}$ | matrix | **State-transition matrix**. $\mathbf{x}_{n+1} = A \cdot \mathbf{x}_n$. |
| $B \in \mathbb{F}_2^{L \times k}$ | matrix | **Output matrix**. $\mathbf{y}_n = B \cdot \mathbf{x}_n$. |
| $T \in \mathbb{F}_2^{L \times L}$ | matrix | **Tempering matrix**. $\mathbf{u}_n = T \cdot \mathbf{y}_n$. |
| $f$ | linear map $\mathbb{F}_2^k \to \mathbb{F}_2^k$ | The F₂-linear **state-transition map**, defined by $f(\mathbf{x}) := A \mathbf{x}$. Iterating $\mathbf{x}_{n+1} = f(\mathbf{x}_n)$ is equivalent to multiplying by $A$. |
| $n$ | integer | **Time step** index in $\mathbf{x}_n$, $\mathbf{y}_n$, $f^n$. Not to be confused with the state size — state size is always $k$. |

## Characteristic polynomial and primitivity

| Symbol | Meaning |
|--------|---------|
| $\chi_f(t)$ | **Characteristic polynomial** of the state-update map $f$. Polynomial of degree $k$ over $\mathbb{F}_2$. |
| $\chi_f$ | Same polynomial, when the variable is clear from context. |
| $\varphi(t)$ | **Selected primitive factor** of $\chi_f$ used for invariant-subspace equidistribution (see [equidistribution spec](equidistribution-spec.md)). |
| $\varphi$ | Same factor, when the variable is clear. |
| $p$ | $\deg \varphi$ — degree of the selected primitive factor. |
| $V$ | The **invariant subspace** $V := \ker \varphi(f) \subseteq \mathbb{F}_2^k$, dimension $p$. |
| $\psi(t)$ | The **complementary factor**: $\psi(t) := \chi_f(t) / \varphi(t)$, degree $k - p$. |
| $2^k - 1$ | Maximum period of a full-period generator (when $\chi_f$ is primitive). |

## Equidistribution

| Symbol | Meaning |
|--------|---------|
| $\ell$ | **Bit resolution** — the number of leading bits kept from each output word in the equidistribution test, $\ell \in \{1, \ldots, L\}$. |
| $t$ | **Tuple dimension** (free variable) — the number of consecutive output words stacked into one $t$-dimensional point. |
| $d(\ell)$ | **Dimension of equidistribution at resolution $\ell$** — the largest $t$ such that the set of all $t$-tuples of consecutive outputs (truncated to $\ell$ leading bits each) covers every $\ell t$-bit pattern equally often. |
| $\lfloor k / \ell \rfloor$ | **Upper bound** on $d(\ell)$ — at most $\lfloor k / \ell \rfloor$ independent $\ell$-bit coordinates can be packed into a $k$-bit state. |
| $\delta(\ell)$ | **Dimension gap** at resolution $\ell$: $\delta(\ell) = \lfloor k / \ell \rfloor - d(\ell) \ge 0$. |
| $\text{SE}$ | **Sum of dimension gaps**, $\text{SE} = \sum_{\ell=1}^{L} \delta(\ell)$. The scalar REGPOLY's search loops minimise. |
| ME | A generator with $\delta(\ell) = 0$ for every $\ell \in \{1, \ldots, L\}$ is **maximally equidistributed**. |

## Resolution sets used by the analyses

| Symbol | Meaning |
|--------|---------|
| $\Psi_{12}$ | **Resolution set** used by the matricial equidistribution method — the resolutions $\ell$ at which $d(\ell)$ is measured. |
| $\Phi_4$ | **Collision-free resolution set**. |

## Tempering optimisation

| Symbol | Meaning |
|--------|---------|
| $b$, $c$ | **Tempering bitmasks** — the optimisable parameters of the tempering map $T$. |
| $S_\ell$ | **Safe mask** at resolution $\ell$: the bitmask of bits that may be flipped without invalidating cached lattice-method state for resolutions $< \ell$. |

## Page-local overrides

A few self-contained pages reuse a symbol with a paper-specific
meaning. Where they do, the page introduces the override at the top:

- [LCP irreducibility](lcp-irreducibility.md) uses $L$ for the
  **linear complexity** (the Berlekamp–Massey output), not the
  output bit width. Local context throughout that page.
- [Antithetic check](antithetic-check.md) uses $d$ for the
  **antithetic test dimension**, not for $d(\ell)$. The two pages
  do not appear together in any analysis chain, so the symbol
  conflict is benign.
- In the MELG `lag` tempering, the C++ constructor parameter
  `int L` is a **word-index lag** (not an output bit width). In
  prose, this is written $\ell_\text{lag}$ to disambiguate; the
  code identifier `L` stays in code font.

## C++ / Python ↔ math symbol map

For readers cross-referencing the API: the math symbol on the left
appears as the parameter name on the right in the implementation.

| Math | C++ parameter | Python parameter |
|------|---------------|------------------|
| $k$  | `int k`       | `gen.k`          |
| $L$  | `int L`       | `gen.L`, `L=32` kwarg |
| $\ell$ | `int l` in loop indices | local integer |
| $t$  | tuple-dim local | local integer |
| $d(\ell)$ | `MatricialEquidResult::ecart[ℓ]` element | `res["gaps"][ℓ]` |
| $\text{SE}$ | `MatricialEquidResult::se` | `res["se"]` |
| $\chi_f$ | `BitVect Generator::char_poly()` | `gen.char_poly()` |

## See also

- [F₂-linear generators](f2-linear-generators.md) — full framework
  context for the symbols above.
- [Equidistribution in 5 minutes](equidistribution-in-5-minutes.md) —
  $d(\ell)$, $\delta(\ell)$, $\text{SE}$, and ME explained without
  algorithmic detail.
- [Equidistribution (design of record)](equidistribution-spec.md) —
  the matricial algorithm in full, including $\chi_f$, $\varphi$,
  and $V$.
