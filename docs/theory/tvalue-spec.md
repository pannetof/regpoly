# $(t, m, s)$-net t-value computation

This page is the design of record for the t-value test in the
regpoly monorepo. It mirrors the structure of
[equidistribution-spec.md](equidistribution-spec.md) and assumes
that page's notation and the conventions in [notation.md](notation.md).

## Setup

Let $m \ge 1$ and $s \ge 1$. An $\mathbb{F}_2$ **digital $(t, m, s)$-net** is
defined by $s$ generating matrices $C_1, \ldots, C_s \in \mathbb{F}_2^{m \times m}$
and produces the $2^m$-point set

$$
\mathcal{P} \;=\; \bigl\{\, (C_1 \cdot \mathbf{i},\; C_2 \cdot \mathbf{i},\; \ldots,\; C_s \cdot \mathbf{i})
\;:\; \mathbf{i} \in \mathbb{F}_2^m \,\bigr\}
\;\subset\; \bigl(\{0, 1\}^m\bigr)^s
$$

where each coordinate $C_j \cdot \mathbf{i}$ is interpreted as an $m$-bit
binary fraction in $[0, 1)$. The **t-value** of $\mathcal{P}$ is the smallest
non-negative integer $t$ such that for every choice of non-negative integers
$(\ell_1, \ldots, \ell_s)$ with $\ell_1 + \cdots + \ell_s = m - t$, the
$2^{m - t}$ points obtained by taking the first $\ell_j$ bits of each
coordinate fall in distinct $\ell_1 \times \cdots \times \ell_s$-dimensional
boxes.

Equivalently (Niederreiter 1987 Theorem 4.10):

$$
t \;=\; m - d_{\max},
\qquad
d_{\max} \;=\; \max\, \Bigl\{\, d \ge 0
\;:\; \exists\, (\ell_1, \ldots, \ell_s) \text{ composition of } d \text{ with }
\operatorname{rank}_{\mathbb{F}_2}\! \bigl(M_{(\ell_1, \ldots, \ell_s)}\bigr) = d \,\Bigr\},
$$

where $M_{(\ell_1, \ldots, \ell_s)}$ stacks the first $\ell_j$ rows of each
$C_j$ vertically. Smaller $t$ is better; $t = 0$ is the strongest possible
(the $(0, m, s)$-net — perfect equidistribution at every resolution).

## Profile over $s$

The regpoly t-value test reports $t(s)$ for every
$s \in \{2, 3, \ldots, s_{\max}\}$ (the case $s = 1$ is trivially $0$ for any
non-degenerate net). The aggregate figure-of-merit is

$$
\mathrm{se} \;=\; \sum_{s = 2}^{s_{\max}} t(s).
$$

Lower is better; $\mathrm{se} = 0$ iff every queried dimension achieved a
$(0, m, s)$-net. The search-rejection contract mirrors equidistribution:
bail with `verified = False` if any per-dimension $t(s) > \delta(s)$ or if
running $\mathrm{se} > \mathrm{maxTSum}$.

## Algorithm (Pirsic–Schmid primal, prefix-sharing DFS + branch-and-bound)

Build the master matrix $M$ once via `GaussMatrix::prepare` with
$\mathrm{indice\_max} = s_{\max}$ and $L = m$. $M$ has $k_g$ rows (here
$k_g = m$ for a single $\mathtt{DigitalNet}$) and $s_{\max} \cdot m$ columns;
block $j$ (1-based) occupies columns $[(j - 1) m,\; j m)$.

The standard t-value-via-$M$ correspondence: with the `prepare` layout,
the submatrix of $M$ restricted to columns
$\bigcup_j [(j - 1) m,\, (j - 1) m + \ell_j)$ has rank equal to the rank of
the stacked top-$\ell_j$-rows-of-$C_j$ matrix (they are transposes of one
another, and rank is transpose-invariant).

For each $s$ from $2$ to $s_{\max}$ and each target $d$ from $m$ down to
$0$, the algorithm searches for some composition $(\ell_1, \ldots, \ell_s)$
of $d$ (with $\ell_j \in [0, m]$) that brings the elimination rank to $d$.
The search is a depth-first walk through the composition tree, with two
cuts that turn the naïve $O(\binom{d+s-1}{s-1}\cdot m\,k_g)$ enumeration
into a fraction of that work:

1. **Prefix sharing.** At DFS depth $j$, the working matrix $M^{(j)}$
   carries the elimination of blocks $0, \ldots, j-1$ at the chosen
   $(\ell_0, \ldots, \ell_{j-1})$. Every child node at depth $j+1$
   clones $M^{(j)}$ and only eliminates the $\ell_{j+1}$ new columns of
   block $j$. The dominant $O(d\,k_g)$ leaf work the naïve version
   used to redo per composition is paid **once** per tree-prefix.

2. **Rank-deficit branch-and-bound.** Each added column raises the
   F$_2$ rank by at most 1, so the best-case rank reachable from a
   subtree rooted at $\bigl(M^{(j)}, \mathrm{rank}(M^{(j)}), r\bigr)$ —
   where $r$ is the amount still to allocate — is
   $\mathrm{rank}(M^{(j)}) + \min(r, k_g - \mathrm{rank}(M^{(j)}))$. If
   that undershoots $d$, the entire subtree is pruned.

   The lex bound

   $$
   \ell_j^{\min} \;=\; \max\bigl(0,\; r - (s - 1 - j) m\bigr),
   \qquad
   \ell_j^{\max} \;=\; \min(r, m),
   $$

   is still applied at every depth to keep the enumeration within
   admissible compositions.

The sibling values $\ell_j = \ell_j^{\min}, \ell_j^{\min} + 1, \ldots,
\ell_j^{\max}$ are enumerated in **increasing order**, so consecutive
siblings differ by exactly one column (one added). That is the
Gray-code-style adjacency the literature names. We do **not** exploit
it for in-place rewind — that would require a Gauss tableau that
supports column removal — and re-clone the parent state per sibling
instead. With prefix sharing the clone cost is dominated by the saved
prefix-elimination work, so the implementation stays simple without
giving up the asymptotic gain.

The first $d$ that admits some full-rank composition fixes
$d_{\max}(s) = d$, and $t(s) = m - d_{\max}(s)$.

Reference implementation:
[`packages/regpoly-cpp/src/analyses/tvalue_runner.cpp`](../api/cpp.html)
(`run_tvalue_profile_schmid` calls `any_composition_full_rank_dfs`,
which folds in both cuts above; the per-block elimination is a thin
wrapper, `eliminate_columns`, around `GaussMatrix::find_pivot`,
`swap_rows`, and `eliminate_column`).

## Niederreiter–Pirsic dual method (stub)

`run_tvalue_profile_dual` is **registered by name but not yet implemented**.
It throws `std::runtime_error` when called. The dual method (Niederreiter &
Pirsic 2001) reformulates the t-value via the minimum NRT weight of the
dual net's parity-check code and may outperform the primal method for
moderate $(m, s)$; it is the natural follow-up when the v1 caps
($m \le 30$, $s_{\max} \le 32$) become binding.

## v1 caps

- $m \le 30$ (Niederreiter polynomial arithmetic in `uint64_t`).
- $m \le 64$ ($\mathtt{DigitalNet}$ base class).
- $s_{\max} \le 32$ (the t-value kernel's branch-and-bound is calibrated
  for this regime).

Beyond these, the dual method or alternative algorithms are required.

## Cross-references

- [equidistribution-spec.md](equidistribution-spec.md) — the
  equidistribution test that runs *unchanged* on a digital net.
- [notation.md](notation.md) — shared notation.
- [SobolNet generator page](../generators/SobolNet.md)
- [NiederreiterF2Gen generator page](../generators/NiederreiterF2Gen.md)

## References

- H. Niederreiter, "Point sets and sequences with small discrepancy,"
  *Monatsh. Math.* **104**, 273–337, 1987.
- H. Niederreiter, "Low-discrepancy and low-dispersion sequences,"
  *J. Number Theory* **30**, 51–70, 1988.
- W. Ch. Schmid, "Improvements and extensions of the 'Salzburg Tables'
  by using irreducible polynomials," in *Monte Carlo and Quasi-Monte
  Carlo Methods 1998*, Springer, 1999.
- G. Pirsic & W. Ch. Schmid, "Calculation of the quality parameter of
  digital nets and application to their construction," *J. Complexity*
  **17**, 827–839, 2001.
- H. Niederreiter & G. Pirsic, "Duality for digital nets and its
  applications," *Acta Arith.* **97**, 173–182, 2001.
- J. Dick & F. Pillichshammer, *Digital Nets and Sequences: Discrepancy
  Theory and Quasi-Monte Carlo Integration*, Cambridge Univ. Press, 2010.
  §4.2 (t-value definitions), §8.1 (Niederreiter construction).
