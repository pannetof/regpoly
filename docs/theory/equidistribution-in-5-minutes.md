# Equidistribution in 5 minutes

A working intuition for what REGPOLY measures, written for people who
want to *use* a generator without first reading Berlekamp–Massey.
For the design-of-record algorithm see [Equidistribution](equidistribution-spec.md).

## The question

> If I take consecutive outputs from this generator and stack them into
> $t$-dimensional points, are those points spread evenly across the
> $t$-dimensional unit cube?

For a Monte Carlo simulation that draws $t$ random numbers per sample
(option pricing, particle transport, etc.), the answer determines
whether your samples cover the integration domain or cluster on a
lattice — which determines whether your estimator is unbiased.

## The metric

For a $t$-dimensional resolution $\ell$ — meaning, take the leading
$\ell$ bits of each output word and pack $t$ consecutive such words
into a single $\ell t$-bit point — count how many of the $2^{\ell t}$
possible points the generator actually visits over its full period, and
how many times each.

A generator is **equidistributed at $(t, \ell)$** if every $\ell t$-bit
point appears exactly the same number of times. The largest $t$ for
which this holds at resolution $\ell$ is the **dimension of
equidistribution** $d(\ell)$.

You want $d(\ell)$ to be as large as possible. An obvious upper bound is
the state size $k$: you cannot have more independent $\ell$-bit
coordinates than the state has bits, so

$$
d(\ell) \le \lfloor k / \ell \rfloor .
$$

The **dimension gap** at resolution $\ell$ is the gap between that
upper bound and what the generator actually achieves:

$$
\delta(\ell) = \lfloor k / \ell \rfloor - d(\ell) \ge 0 .
$$

A gap of zero means the generator is as good as it can possibly be at
that resolution. The **total dimension defect** is the sum across
resolutions:

$$
\text{SE} = \sum_{\ell=1}^{L} \delta(\ell) ,
$$

where $L$ is the output word width. SE is the single scalar REGPOLY's
search loops use to rank candidate parameter sets.

## What "ME" means

A generator with $\delta(\ell) = 0$ for every $\ell$ — every dimension
gap is zero — is **maximally equidistributed (ME)**. ME at $L = 32$
means the generator is doing the best a 32-bit word size lets it do.
ME is the ceiling; "$\text{SE} = 3$" is a small distance from that
ceiling.

In REGPOLY's search loops:

- A *target* is usually expressed as `mse: N` — search until you find
  a parameter set with $\text{SE} \le N$.
- A *test* prints the full per-$\ell$ gap profile so you can see *where*
  the gaps land — gaps at small $\ell$ (low-resolution coordinates) are
  more harmful than gaps at large $\ell$ for most applications.

## Where the numbers come from

REGPOLY computes $d(\ell)$ via a matricial method:

1. Build the $\ell t \times k$ matrix that maps the generator's state
   to the first $t$ output words truncated to $\ell$ leading bits each.
2. Reduce it via Gaussian elimination over $\mathbb{F}_2$.
3. The rank of that matrix at the largest $t$ where it equals $\ell t$
   is $d(\ell)$.

For full-period generators this is straightforward. For generators
whose characteristic polynomial is reducible (no full period
guaranteed), REGPOLY first factors the characteristic polynomial,
projects onto the largest invariant subspace, and applies the same
matricial computation there — that is what the [equidistribution
spec](equidistribution-spec.md) documents in detail.

## Running it yourself

Python:

```python
from regpoly.library import Catalog
from regpoly.core.generator import Generator
from regpoly.analyses.pis import analyze_single_generator

cat = Catalog("docs/library"); cat.load()
_, params = cat.generator("mt19937")
gen = Generator.create(params.family, L=32, **params.params)
res = analyze_single_generator(gen)
print("SE =", res["se"])           # total defect across all resolutions
print("gaps[:10] =", res["gaps"][:10])
```

C++ — see [C++ usage](../usage/cpp.md) for the equivalent runnable
snippet with `regpoly::core::run_matricial_equidistribution`.

## See also

- [F₂-linear generators](f2-linear-generators.md) — the matrices behind
  $A$, $B$, $T$ and what "linear over $\mathbb{F}_2$" really means.
- [Equidistribution spec](equidistribution-spec.md) — full algorithm,
  including the reducible-$\chi$ case.
- [Search format](search_format.md) — the YAML schema where you
  express a target `mse` (the SE upper bound your search will accept).
