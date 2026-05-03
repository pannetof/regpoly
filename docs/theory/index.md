# Theory

Mathematical and algorithmic background for F₂-linear generators and the analyses REGPOLY performs.

- **[F₂-linear generators](f2-linear-generators.md)** — general framework: state, recurrence, output, characteristic polynomial, full-period condition.
- **[Equidistribution](equidistribution-spec.md)** — design of record for matricial equidistribution computation on F₂-linear generators that may not be full-period.
- **[Antithetic check](antithetic-check.md)** — algorithm for testing local antitheticity of a linear RNG point set.
- **[Lattice method (Couture–L'Ecuyer)](lattice_method_couture_lecuyer.md)** — dual lattice basis + Lenstra reduction.
- **[Lattice method (Harase)](lattice_method_harase.md)** — Harase–Matsumoto–Saito primal reduction.
- **[Tempering optimization](tempering_optimization.md)** — incremental basis caching for tempering parameter optimization.
- **[LCP irreducibility](LCP Irreducibility.md)** — linear complexity profile / characteristic polynomial irreducibility.
- **[Search format](search_format.md)** — YAML config schema for search drivers.
- **[Tempering optimizer notes](optimzer_temper.md)** — additional notes on the tempering optimizer.

> **Note**: `f2-linear-generators.md` is the umbrella reference; the others go into specific algorithmic detail.
