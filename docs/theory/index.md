# Theory

Mathematical and algorithmic background for F₂-linear generators and the analyses REGPOLY performs.

- **[Notation](notation.md)** — canonical math symbols ($k$, $L$, $\ell$, $t$, $d(\ell)$, $\delta(\ell)$, …) used across every page. Read first when picking up a new section.
- **[Equidistribution in 5 minutes](equidistribution-in-5-minutes.md)** — beginner on-ramp: what $d(\ell)$, gaps, SE, and "ME" mean, why the metric matters, no Berlekamp-Massey.
- **[F₂-linear generators](f2-linear-generators.md)** — general framework: state, recurrence, output, characteristic polynomial, full-period condition.
- **[Equidistribution (design of record)](equidistribution-spec.md)** — matricial equidistribution computation on F₂-linear generators that may not be full-period.
- **[Antithetic check](antithetic-check.md)** — algorithm for testing local antitheticity of a linear RNG point set.
- **[Lattice method (Couture–L'Ecuyer)](lattice_method_couture_lecuyer.md)** — dual lattice basis + Lenstra reduction.
- **[Lattice method (Harase)](lattice_method_harase.md)** — Harase–Matsumoto–Saito primal reduction.
- **[Tempering optimization](tempering_optimization.md)** — incremental basis caching for tempering parameter optimization.
- **[LCP irreducibility](lcp-irreducibility.md)** — linear complexity profile / characteristic polynomial irreducibility.
- **[Search format](search_format.md)** — YAML config schema for search drivers.
- **[Tempering optimizer notes](optimizer_temper.md)** — additional notes on the tempering optimizer.

> **Note**: `f2-linear-generators.md` is the umbrella reference; the others go into specific algorithmic detail.

```{toctree}
:hidden:

notation
equidistribution-in-5-minutes
f2-linear-generators
equidistribution-spec
antithetic-check
lattice_method_couture_lecuyer
lattice_method_harase
tempering_optimization
lcp-irreducibility
search_format
optimizer_temper
```
