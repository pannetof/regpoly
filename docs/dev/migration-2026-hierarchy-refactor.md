# Migration — `F2LinearSource` / `Recurrence` hierarchy refactor (2026-05)

Between commits `47c24e6` (v2.0.0) and `46a8815` (2026-05) the C++
generator hierarchy was reorganised to separate three concerns that
the original `Generator` class conflated:

- **Tempering ownership** — chains now live on each source
  intrinsically, not on an external "slot".
- **Recurrence vs digital net** — `Recurrence` (state-evolution PRNGs)
  and `DigitalNet` (point-set sources) are now siblings under
  `F2LinearSource`, instead of `DigitalNet` borrowing the `Generator`
  surface.
- **Single source vs combination** — combinations are composition
  (`CombinedF2LinearSource`), not inheritance, so a J-component combo
  no longer pretends to be a single recurrence.

The architectural target:

```
ITestable                              ← kernel-facing pure interface
   ▲
   ├─── F2LinearSource (abstract)      ← single source; owns a TemperingChain
   │       ▲
   │       ├─── Recurrence (abstract)  ← was `Generator`; state-evolution PRNGs
   │       │       ├── MTGen, WELLGen, SFMTGen, DSFMTGen, MTGPGen,
   │       │       ├── MELGGen, RMT64Gen, TinyMT32Gen, TauswortheGen,
   │       │       ├── PolyLCGGen, F2wBaseGen, TGFSRGen, XoroshiroGen,
   │       │       └── XoshiroGen, MarsaXorshiftGen, CellularAutomataGen
   │       │
   │       └─── DigitalNet (abstract)
   │               ├── SobolNet, NiederreiterF2Net
   │
   └─── CombinedF2LinearSource         ← composition: holds [F2LinearSource₁ … F2LinearSourceₙ]
```

`CombinedF2LinearSource` is the only combined wrapper. It is
**composition**: implements `ITestable` directly, holds J
`Recurrence` components, and exposes `components()` for kernels that
need per-component state evolution. The lattice / matricial kernels
walk the J components themselves — there is no fake "combined
Recurrence" state pretending the wrapper is a single recurrence. The
legacy `CombinedGenerator` class was retired on 2026-05-23 along with
the brief intermediate "promote to Recurrence" variant; the only
combined wrapper that survives is composition-shaped.
`build_combined_from_enumerator` rejects digital-net pools at build
time (a digital net has no state in the recurrence sense).

The `Generator` C++ name is gone entirely (no `using Generator =
Recurrence;` alias). The `_cpp.Generator` Python attribute is also
gone. Code that still spells `Generator` won't compile / import.

## New types

| Type | Header | Role |
|---|---|---|
| `ITestable` | `core/i_testable.h` | Pure interface every analysis kernel accepts. `k`, `L`, `name`, `display_str`, `init`, `next`, `get_output`, `copy`, `sources`, `default_test_method`. |
| `F2LinearSource` | `core/f2_linear_source.h` | Abstract base for a single F_2-linear bit source. Owns a `TemperingChain` by value. Implements `ITestable`. |
| `TemperingChain` | `core/tempering_chain.h` | First-class type wrapping `std::vector<std::unique_ptr<Transformation>>` with `apply`, `randomize_params`, deep-copy. Replaces the `using ComponentTempering = …;` typedef. |
| `CombinedF2LinearSource` | `generators/combined_f2_linear_source.h` | Composition wrapper: implements `ITestable` directly, holds `vector<unique_ptr<Recurrence>>`, XORs each component's already-tempered output. Exposes `components()` so lattice / matricial kernels walk per-component state. Absorbed and replaced the legacy `CombinedGenerator` on 2026-05-23. |
| `Recurrence` | `core/generator.h` (file kept) | Renamed from `Generator`. Recurrence-driven PRNG abstract base. Sibling of `DigitalNet` under `F2LinearSource`. |
| `F2LinearSourcePool` | `core/combo_enumerator.h` | Renamed from `Component`. Per-slot pool of candidate F_2-linear sources used by `ComboEnumerator`. Storage widened to `vector<unique_ptr<F2LinearSource>>` so digital nets can enter the pool. |

## Per-symbol rename table

### C++ class renames

| Removed / renamed | Replacement |
|---|---|
| `regpoly::core::Generator` | `regpoly::core::Recurrence`. No alias — Phase 5.7 (2026-05) dropped the transitional `using Generator = Recurrence;`. Code that still spells `Generator` won't compile. |
| `regpoly::core::Component` | `regpoly::core::F2LinearSourcePool`. |
| `Component::GenPool` typedef | `F2LinearSourcePool::SourcePool` typedef (now `vector<unique_ptr<F2LinearSource>>`). |

### C++ method renames (`F2LinearSourcePool`)

| Removed / renamed | Replacement |
|---|---|
| `Component::add_gen(const Generator&)` | `F2LinearSourcePool::add_source(const F2LinearSource&)` — accepts any F2LinearSource subclass, including `DigitalNet`. |
| `Component::gen_at(int)` | `F2LinearSourcePool::source_at(int)` — returns `F2LinearSource&`. |
| `Component::active_gen()` | `F2LinearSourcePool::active_source()` — returns `F2LinearSource&`. |
| `Component::nb_gen()` | `F2LinearSourcePool::nb_sources()` |
| `Component::current_gen()` | `F2LinearSourcePool::current_index()` |
| `Component::set_current_gen(int)` | `F2LinearSourcePool::set_current_index(int)` |
| `ComboEnumerator::component(int)` | `ComboEnumerator::pool(int)` — returns `F2LinearSourcePool&`. |
| `ComboEnumerator::at(int)` | Unchanged name. Return type widened to `F2LinearSource&`. |

### C++ inheritance changes

| Was | Now |
|---|---|
| `class MTGen : public Generator` (and every other recurrence-driven family) | `class MTGen : public Recurrence` (alias-equivalent during the transition). |
| `class DigitalNet : public Generator` | `class DigitalNet : public F2LinearSource`. Lost the recurrence-specific surface (`state`, `char_poly`, `simd_*`). Public `state()` accessor was added on `DigitalNet` itself. |
| `class CombinedGenerator : public Generator` | Deleted (2026-05-23). Replaced by `class CombinedF2LinearSource : public ITestable` — composition, not inheritance. Holds J `Recurrence` components; exposes `components()` so the lattice / matricial kernels walk per-component state directly (no "combined recurrence" fiction). `combined.h` and `combined.cpp` removed; consumers `#include "combined_f2_linear_source.h"`. |

### Removed / hidden virtuals on `Recurrence`

| Removed | Replacement |
|---|---|
| `virtual std::vector<std::vector<Transformation*>> tempering_chains() const` | Chain ownership moved onto each `F2LinearSource::tempering_` intrinsically in Phase 2 step 6a. Read with `gen.tempering()` per-source. |
| External `trans` list arguments on every lattice kernel signature | Dropped in Phase 5. Chains live on each source; kernels read tempered bits via `gen.get_output()`. |

### Python module attributes

| Existing | New (preferred) |
|---|---|
| `_cpp.Generator` | `_cpp.Recurrence`. Phase 5.7 dropped the transitional `_cpp.Generator` alias — `getattr(_cpp, "Generator")` raises `AttributeError`. |
| `_cpp.Generateur` (legacy French alias) | Also gone (Phase 5.7). |
| `_cpp.Component` | `_cpp.F2LinearSourcePool`. Phase 5.8 dropped the transitional alias — `getattr(_cpp, "Component")` raises `AttributeError`. The legacy method aliases (`add_gen`, `gen_at`, `active_gen`, `nb_gen`, `current_gen`, `set_current_gen`, `ComboEnumerator.component`) are also gone. |
| `_cpp.CombinedGenerator` | Deleted (2026-05-23). Replaced by `_cpp.CombinedF2LinearSource`, which inherits `_cpp.ITestable` (composition) and accepts only Recurrence components. `getattr(_cpp, "CombinedGenerator")` raises `AttributeError`. |
| `_cpp.F2LinearSource` (new) | Marker class exposing `k`, `L`, `name`, `display_str`, `init`, `next`, `get_output`, `copy`, `default_test_method`. |
| `_cpp.ITestable` (new) | Bare marker class. The kernel-facing interface. |

### Kernel API widening (signature changes)

| Kernel | Was | Now |
|---|---|---|
| `run_matricial_equidistribution` | `(const Generator& gen, …)` | `(const ITestable& gen, …)` (Phase 3.1). |
| `run_collision_free` | `(const Generator& gen, …)` | `(const ITestable& gen, …)`. |
| `run_tvalue_profile_schmid` / `_dual` | `(const Generator& gen, …)` | `(const ITestable& gen, …)`. |
| `run_tuplets` | `(const Generator& gen, …)` | `(const ITestable& gen, …)`. |
| `GaussMatrix::prepare` | Took an external per-component `trans` list. | Takes `vector<const F2LinearSource*>` only; chains are read intrinsically via `get_output()`. |
| `test_me_lat` / `test_me_harase` / `test_me_notprimitive` / `test_me_notprimitive_simd` / `compute_kv` (lattice) | `(const Generator& gen, …)` | `(const ITestable& gen, …)` (Phase 5.4). Adapter unpacks `gen.sources()`, dynamic-casts each to `Recurrence*`, throws if any is a non-Recurrence (DigitalNet). Inner vector-form kernels stay typed on `vector<Generator*>`. |
| `HaraseRankCache` / `TemperingOptimizerCache` ctors | `(const Generator&, int, int)` | `(const ITestable&, int, int)` (Phase 5.4). |
| `EquidistributionMethod::run` | `(const Generator&, …)` | `(const ITestable&, …)` (Phase 5.4). All registered methods (`matricial`, `lattice`, `harase`, `notprimitive`, `simd_notprimitive`, `nothing`) widened in lockstep. |
| `SearchPredicate::run` (and subclasses `EquidistributionPredicate`, `CollisionFreePredicate`, `TupletsPredicate`) | `(SeekIterResult&, const Generator&, int, int)` | `(SeekIterResult&, const ITestable&, int, int)` (Phase 5.4). |
| `build_combined_from_enumerator` | `→ unique_ptr<CombinedGenerator>` | `→ unique_ptr<ITestable>` (returns a `CombinedF2LinearSource`). Throws `std::invalid_argument` if any active pool source isn't a Recurrence — combined XOR requires Recurrence components only (Phase 5.5 dispatch + 2026-05-23 CombinedGenerator retirement). |

## Result-struct unification

`MatricialEquidResult` and `MeLatResult` collapsed into a single
`EquidistributionResult { std::vector<int> ecart; int se; bool verified; }`.
Every lattice kernel constructs the unified type; lattice methods set
`verified = true` by construction (the primitive precondition makes
the result unconditionally verified).

## Migration recipes

### Source-side C++ — registering a new recurrence family

```cpp
class MyGen : public regpoly::core::Recurrence { … };
```

The `clone_recurrence()` non-virtual forwarder gives you a
`unique_ptr<Recurrence>` clone with the correct static type via a
`static_cast` + debug-only `dynamic_cast` assert.

### Source-side C++ — registering a digital net family

```cpp
class MyNet : public regpoly::core::DigitalNet {
    // Override raw_output(), copy(), name(), display_str(),
    // default_test_method(). The F2LinearSource layer handles
    // tempering composition. DigitalNet supplies state() / j_ /
    // s_max_ / next() / init() defaults.
};
```

### Building a combined source

Every component must be a `Recurrence` — `CombinedF2LinearSource` is
the Recurrence-typed XOR wrapper, and a digital net has no state in
the recurrence sense.

```cpp
std::vector<std::unique_ptr<Recurrence>> comps;
comps.push_back(my_recurrence_a->clone_recurrence());
comps.push_back(my_recurrence_b->clone_recurrence());
CombinedF2LinearSource cs(std::move(comps), /*Lmax=*/32);
```

Pass `cs` to any `ITestable`-accepting kernel: `run_tvalue_profile`,
`run_matricial_equidistribution`, `run_tuplets`, `run_collision_free`.
The lattice / matricial kernels accept it too — `CombinedF2LinearSource`
inherits `Recurrence`.

### Search-side C++ — setting up a pool

```cpp
ComboEnumerator comb(2, /*Lmax=*/32);
comb.pool(0).add_source(*my_lfsr_a);
comb.pool(0).add_source(*my_lfsr_b);
comb.pool(1).copy_pool_from(comb.pool(0));   // shared-pool C(n, k)
if (comb.reset()) {
    do {
        auto combined = build_combined_from_enumerator(comb);
        // combined is a unique_ptr<ITestable> aliasing a
        // CombinedF2LinearSource. Throws invalid_argument if any
        // slot's active source is not a Recurrence.
    } while (comb.next());
}
```

Digital nets **may** be added to a pool for single-source iteration
use cases (the `F2LinearSourcePool` storage was widened to accept
them), but `build_combined_from_enumerator` itself requires every
active source to be a Recurrence and throws `std::invalid_argument`
otherwise. A digital net has no state in the recurrence sense and
cannot participate in a combined XOR build.

### Python — Generator wrappers and `make_combined`

```python
from regpoly import Generator, make_combined
from regpoly.analyses.tvalue_test import TValueTest

g_mt_a = Generator.create("MTGen", L=32, w=32, r=624, m=397, p=31, a=0x9908B0DF)
g_mt_b = Generator.create("MTGen", L=32, w=32, r=624, m=397, p=31, a=0x9908B0DF)

# All-Recurrence combo → CombinedF2LinearSource.
combined = make_combined(g_mt_a, g_mt_b, Lmax=32)

TValueTest(s_max=4, max_t_sum=100).run(combined)
```

`make_combined` requires every component to be a Recurrence. Passing
a digital net raises `TypeError`.

## Why

- The original `Generator` was the abstract base for *both*
  recurrence-driven PRNGs and digital nets. Both share an interface
  (`init` / `next` / `get_output`) but they mean different things —
  for a PRNG, `next()` advances state via an F_2-linear recurrence;
  for a digital net, it advances the coordinate index `j`. Bundling
  them under one name made the recurrence-specific surface
  (`state`, `char_poly`, `transition_matrix`, `simd_*`) leak into the
  digital-net family.
- Tempering chains were stored half on the `Generator`'s
  `tempering_chains()` accessor and half on the search-side `Component`
  slot. The chain is conceptually a property of the source, not of
  the slot.
- `CombinedGenerator` IS-A `Generator` — inheritance — conflated
  "single F_2-linear source" with "XOR of J sources". The new
  `CombinedF2LinearSource` (which absorbed `CombinedGenerator`'s
  role) keeps the inheritance shape under the architecturally-
  correct name; the lattice / matricial kernels still need it for
  Krylov-based χ-recovery.

The post-refactor architecture has a single kernel-facing interface
(`ITestable`), a clean single-source/composition split, and
intrinsic tempering ownership. Existing source code keeps compiling
through the type aliases; new code can adopt the architecturally
honest names.

## See also

- [`docs/dev/migration-2026-combination-removal.md`](migration-2026-combination-removal.md) — the prior refactor that
  retired the Python `Combination` / `Component` classes.
- [`docs/dev/architecture.md`](architecture.md) — overall system architecture.
