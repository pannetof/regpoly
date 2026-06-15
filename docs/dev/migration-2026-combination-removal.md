# Migration — dropping Python `Combination` (2026)

> **Follow-up:** later in May 2026, a second migration retired the
> `_cpp.CombinedF2LinearSource` class and the `Generator` C++ base class
> in favour of the `ITestable` / `F2LinearSource` / `Recurrence`
> hierarchy (see
> [migration-2026-hierarchy-refactor.md](migration-2026-hierarchy-refactor.md)).
> The runtime XOR-combined type is now `_cpp.CombinedF2LinearSource`
> (composition, not inheritance). Below the original combination-
> removal text — historically accurate at the time, but the
> `CombinedGenerator` name has since been swept to
> `CombinedF2LinearSource`.

In May 2026 the regpoly monorepo retired the Python `Combination` and
`Component` classes and renamed the C++ search-loop iterator
`regpoly::core::Combination` → `regpoly::core::ComboEnumerator`. The
public Python API around generators is now built on two C++ types
only (names listed here are the post-hierarchy-refactor names):

- `_cpp.CombinedF2LinearSource` — an `ITestable` composition wrapper
  that XORs J components (optionally with per-component tempering).
  This is the *runtime* object every `*Test.run(...)` consumes.
- `_cpp.ComboEnumerator` — the stateful search-loop iterator over
  the cartesian product of J component pools. Used only by `Seek`
  / `TemperingSearch`.

The previous three-class soup (`_cpp.CombinedF2LinearSource`, Python
`Combination`, C++ `Combination`) is gone; nothing in the codebase
constructs the Python `Combination` middleman any more.

## Per-symbol migration

| Removed / renamed                                              | Replacement                                                           |
|----------------------------------------------------------------|------------------------------------------------------------------------|
| `regpoly.Combination`                                          | `regpoly.make_combined(*gens, trans=…, Lmax=…)` (returns a `_cpp.CombinedF2LinearSource`) — or pass a bare `Generator` directly to `*.run(...)` for J=1 / no tempering. |
| `regpoly.Component`                                            | (no replacement — the slot abstraction moved into the C++ enumerator.) |
| `regpoly.Combinaison`                                          | The French alias was deleted with `Combination`.                       |
| `Combination.single(gen)`                                      | `gen` directly — pass to `*.run(...)` unchanged.                       |
| `Combination(J=N, Lmax=L)` + `add_gen` + `reset`               | `make_combined(g1, g2, …, Lmax=L)` for J ≥ 1 analyses, or `regpoly.io.combo_builder.build_cpp_enumerator(gen_pools, tempering_pools, Lmax)` for *search* (pool-of-many) flows. |
| `Combination.CreateFromFiles([[g]], Lmax, [chain])` (pool-of-one) | `make_combined(g, trans=chain, Lmax=Lmax)`.                            |
| `Combination.CreateFromFiles(pool_lists, Lmax, chains)` (multi-pool) | `regpoly.io.combo_builder.build_cpp_enumerator(pool_lists, chains, Lmax)`. |
| `_cpp.Combination`                                             | `_cpp.ComboEnumerator`.                                                |
| `*.run(C: Combination)`                                        | `*.run(gen: Generator \| CombinedF2LinearSource)`. The shim that accepted the legacy Python `Combination` is gone. |

## Test runner changes

`AbstractTest.run(gen_or_C)` now accepts:

- A Python `Generator` wrapper (`Generator.create(...)`).
- A bare `_cpp.Recurrence` / `_cpp.F2LinearSource` (any concrete C++
  generator, including a `_cpp.CombinedF2LinearSource` from
  `make_combined`).
- A duck-typed C++ source object with `.k()` / `.L()` / `.get_output()`.

Anything else raises `TypeError` with a hint pointing at
`Generator(cpp_gen)` / `Generator.create(...)` / `make_combined(...)`.

## Search-loop changes

`Seek.from_yaml(...)` / `TemperingSearch.run()` build the C++
`ComboEnumerator` directly from per-slot generator/tempering pools.
The old `Seek._comb = ...` private-attribute injection has been
replaced by a public method:

```python
seek = Seek()
seek.set_enumerator(
    cpp_enum=cpp_enum,         # _cpp.ComboEnumerator
    gen_pools=gen_pools,        # list[list[Generator]]
    tempering_pools=temp_pools, # list[list[Transformation]]
    Lmax=Lmax,
    tests=tests,
    nbtries=nbtries,
)
seek.run()
```

The C++ side (`Seek::run` in `regpoly-cpp`) owns the per-combo
iteration; Python keeps display + result serialisation only.

## C++ rename

The header and source files were renamed:

- `packages/regpoly-cpp/src/include/core/combination.h` → `combo_enumerator.h`
- `packages/regpoly-cpp/src/core/combination.cpp`        → `combo_enumerator.cpp`

The function `build_combined_from_combination(const Combination&)`
became `build_combined_from_enumerator(const ComboEnumerator&)`.

C++ consumers of `regpoly::core::Combination` (search-loop callbacks
in `seek_search.h` / `tempering_search.h`, the YAML config builder,
`regpoly-cli`) now take `regpoly::core::ComboEnumerator&` — a
mechanical rename with the same surface (`J()`, `Lmax()`, `k_g()`,
`L()`, `reset()`, `next()`, `at(j)`, `component(j)`, `exhausted()`).

## Examples

### Single generator (most common new shape)

```python
from regpoly import Generator
from regpoly.analyses.equidistribution_test import EquidistributionTest

gen = Generator.create("MTGen", L=32, w=32, r=624, m=397, p=31, a=0x9908B0DF)
res = EquidistributionTest(L=32, delta=[10**9]*33, mse=10**9).run(gen)
print(res.se, res.verified)
```

### Combined generator (J ≥ 2 or with tempering)

```python
from regpoly import Generator, make_combined
from regpoly.analyses.equidistribution_test import EquidistributionTest

g1 = Generator.create("TauswortheGen", L=32, k=31, nb_terms=3,
                      poly=[0, 6, 31], s=18, quicktaus=True)
g2 = Generator.create("TauswortheGen", L=32, k=29, nb_terms=3,
                      poly=[0, 2, 29], s=2, quicktaus=True)
combined = make_combined(g1, g2, Lmax=32)
res = EquidistributionTest(L=32, delta=[10**9]*33, mse=10**9).run(combined)
```

### Search loop (was `Combination.CreateFromFiles` with multi-candidate pools)

```python
from regpoly.io.combo_builder import build_cpp_enumerator

cpp_enum = build_cpp_enumerator(gen_pools, tempering_pools, Lmax=32)
# Pass `cpp_enum` directly to `_cpp.run_seek_search(...)` or use
# `Seek.set_enumerator(...)` to wire it into the regpoly search driver.
```

## Why

- `Combination` did double duty as *analysis input* and *search
  config* — two unrelated responsibilities sharing a class name.
- C++ `regpoly::core::Combination` and Python `Combination` were
  separately useful but collided in name in user-facing docs (the
  Python one was a passive container, the C++ one a stateful
  iterator).
- The C++ kernels have always accepted a bare source reference (now
  `const ITestable&` / `const F2LinearSource&`) — the Python
  middleman added a `CombinedF2LinearSource` wrap on every test call
  that was redundant when `J == 1` with no tempering.

The post-2026 API gets the kernels' shape parity on the Python side
and removes the name collision in one cut.
