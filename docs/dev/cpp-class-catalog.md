# C++ class & interface catalog

A snapshot of every public type defined in `packages/regpoly-cpp/src/include/`,
grouped by role. The **Rename plan** column lists any pending or
deferred rename per the v2.0 / hierarchy-refactor plans; "—" means
none planned. Compound `*Gen` family classes that match published
literature are intentionally kept.

Generated 2026-05-22 (commit `ad4a47f`); last refreshed 2026-05-24 (commit `22e8f80`).

## Core interfaces / abstract bases

| Type | Header | Kind | Inherits | Role | Rename plan |
|---|---|---|---|---|---|
| `ITestable` | `core/i_testable.h` | abstract class | — | Pure interface every analysis kernel consumes (`k`, `L`, `name`, `init`, `next`, `get_output`, `copy`, `sources`, `default_test_method`). | — |
| `F2LinearSource` | `core/f2_linear_source.h` | abstract class | `ITestable` | Single F_2-linear bit source; owns a `TemperingChain` by value; supplies `clone_source()` typed-forwarder. | — |
| `Recurrence` | `core/generator.h` | abstract class | `F2LinearSource` | Recurrence-driven PRNG (state evolution under an F_2-linear map). Supplies `state()`, `char_poly()`, `simd_*`, `clone_recurrence()`. Renamed from `Generator` in Phase 5.1. | — |
| `DigitalNet` | `core/digital_net.h` | abstract class | `F2LinearSource` | Digital net (point-set source over s_max dimensions; `next()` advances coordinate index `j`). Re-anchored to `F2LinearSource` in Phase 3. | — |
| `Transformation` | `core/transformation.h` | abstract class | — | Single tempering step (bitvector → bitvector). | — |
| `TemperingChain` | `core/tempering_chain.h` | concrete class | — | Owns `vector<unique_ptr<Transformation>>`. First-class type since Phase 1. | — |

## Recurrence families (concrete PRNG classes)

Every class below inherits `Recurrence`. The `*Gen` suffix marks
"the canonical class for this PRNG family" and is kept by design
(established literature names; renaming propagates to every notebook
and paper-cited parameter set).

| Type | Header | Family | Rename plan |
|---|---|---|---|
| `MTGen` | `generators/mt.h` | Mersenne Twister (MT19937 et al.) | — |
| `WELLGen` | `generators/well.h` | WELL family (well1024, well19937a, …) | — |
| `SFMTGen` | `generators/sfmt.h` | SIMD-oriented Fast MT | — |
| `DSFMTGen` | `generators/dsfmt.h` | Double-precision SFMT | — |
| `MTGPGen` | `generators/mtgp.h` | MT for graphics processors | — |
| `MELGGen` | `generators/melg.h` | Maximally-Equidistributed Long-period | — |
| `RMT64Gen` | `generators/rmt64.h` | Reflective MT64 | — |
| `TinyMT32Gen` | `generators/tinymt32.h` | Tiny MT (32-bit) | — |
| `TauswortheGen` | `generators/tausworthe.h` | Tausworthe LFSR | — |
| `PolyLCGGen` | `generators/polylcg.h` | Polynomial LCG | — |
| `TGFSRGen` | `generators/tgfsr.h` | Twisted GFSR | — |
| `XoroshiroGen` | `generators/xoroshiro.h` | Xoroshiro family | — |
| `XoshiroGen` | `generators/xoshiro.h` | Xoshiro family | — |
| `MarsaXorshiftGen` | `generators/marsaxorshift.h` | Marsaglia Xorshift (subsumes legacy XorShift128) | — |
| `CellularAutomataGen` | `generators/cellular_automata.h` | Cellular-automaton recurrence | — |
| `F2wBaseGen` | `generators/f2w_base.h` | Abstract base for F_{2^w} recurrences | — |
| `F2wLFSRGen` | `generators/f2w_lfsr.h` | F_{2^w}-LFSR | — |
| `F2wPolyLCGGen` | `generators/f2w_polylcg.h` | F_{2^w} polynomial LCG | — |

## Digital nets (concrete classes)

| Type | Header | Description | Rename plan |
|---|---|---|---|
| `SobolNet` | `generators/sobol.h` | Sobol net via embedded Joe–Kuo 2003 table. | — |
| `NiederreiterF2Net` | `generators/niederreiter_f2.h` | Niederreiter F_2 net (irreducibles generated on the fly). Renamed from `NiederreiterF2Gen` in Phase 4. | — |
| `SobolDirNumbers` | `generators/sobol.h` | Plain struct: Sobol direction-number triple. | — |

## Combined wrappers (J-component XOR)

| Type | Header | Kind | Description | Rename plan |
|---|---|---|---|---|
| `CombinedF2LinearSource` | `generators/combined_f2_linear_source.h` | `class : ITestable` | The combined XOR wrapper — **composition**, not inheritance. Holds J `F2LinearSource` components (Recurrence PRNGs **and/or** DigitalNet sources); implements `ITestable` directly. Exposes `components() → vector<F2LinearSource*>` for kernels that consume any source, and `recurrence_components("kernel_name") → vector<Recurrence*>` for kernels that genuinely need recurrence-state evolution — the latter throws `std::invalid_argument` if any component is a DigitalNet. Each component carries its own intrinsic tempering chain; the wrapper has no top-level chain. | — |

## Search-loop iteration

| Type | Header | Description | Rename plan |
|---|---|---|---|
| `F2LinearSourcePool` | `core/combo_enumerator.h` | Per-slot pool of candidate F_2-linear sources + a per-slot tempering-chain template. Storage: `vector<unique_ptr<F2LinearSource>>` (widened in Phase 5.2-storage so digital nets can enter). Renamed from `Component` in Phase 5.2. | — |
| `ComboEnumerator` | `core/combo_enumerator.h` | Stateful iterator over the cartesian product of J `F2LinearSourcePool`s with identity-uniqueness + shared-pool C(n,k) selection. | — |
| `ParameterSpaceEnumerator` | `core/parameter_space_enumerator.h` | Iterator over the parameter space of a single PRNG family. Renamed from `GenEnumerator` (Tier 3 follow-up); header file renamed in lockstep. | — |
| `ParamBag` | `core/params.h` | Heterogeneous string-keyed parameter bag (ints, bools, strings, int_vecs, uint_vecs). Renamed from `Params` (Tier 3 follow-up). | — |
| `ParamSpec` | `core/params.h` | Schema for a single parameter (name + type + optional range). | — |
| `GeneratorRegistry` | `core/generator_registry.h` | Factory registry keyed by PRNG family name. | — |

## Transformations (tempering)

| Type | Header | Inherits | Description | Rename plan |
|---|---|---|---|---|
| `TemperMKTrans` | `transforms/temper_mk.h` | `Transformation` | Mask/shift tempering step. | — |
| `PermutationTrans` | `transforms/permutation.h` | `Transformation` | Bit-permutation tempering step. | — |
| `LaggedTemperingTrans` | `transforms/lag_mask.h` | `Transformation` | Lagged-bit-mask tempering. Renamed from `LaggedTempering` in Phase 4. | — |
| `TransformationRegistry` | `core/transformation_registry.h` | — | Factory registry keyed by transformation type. | — |

## Equidistribution methods (polymorphic dispatch)

| Type | Header | Kind | Description | Rename plan |
|---|---|---|---|---|
| `EquidistributionMethod` | `analyses/equidistribution_method.h` | abstract class | Polymorphic method interface; concrete subclasses wrap a lattice / matricial kernel. | — |
| `MatricialMethod` / `LatticeMethod` / `HaraseMethod` / `NotPrimitiveMethod` / `SimdNotPrimitiveMethod` / `NothingMethod` | `analyses/equidistribution_method.cpp` | concrete subclasses | One wrapper per registered method name. | — |
| `EquidistributionMethodRegistry` | `analyses/equidistribution_method.h` | concrete class | Name-string registry for the above. Renamed from `MethodRegistry` in Phase 4. | — |
| `EquidistributionMethodResult` | `analyses/equidistribution_method.h` | struct | Uniform result: per-resolution gap, cumulative SE, verified flag. | — |

## Lattice / dual-lattice

| Type | Header | Description | Rename plan |
|---|---|---|---|
| `HaraseRankCache` | `lattice/me_harase.h` | Per-`v` PIS-reduction cache for the Harase lattice method. Renamed from `PISCache` in Phase 4 (clearer than the algorithm-jargon original). | — |
| `TemperingOptimizerCache` | `lattice/temper_optimizer.h` | Per-component cache for the tempering optimiser. Renamed from `TemperOptCache` in Phase 4. | — |
| `DualLatticeBase` | `lattice/dual_lattice.h` | Lenstra dual-lattice reduction. | — |
| `DualLatticeBasisRow` | `lattice/dual_lattice.h` | One row of a `DualLatticeBase`. Renamed from `PolVect` in Phase 4. | — |

## Algebra primitives

| Type | Header | Description | Rename plan |
|---|---|---|---|
| `BitVect` | `algebra/bitvect.h` | Packed bit-vector with F_2 arithmetic (XOR, shifts, get/set bit). Domain-canonical. | — |
| `GaussMatrix` | `algebra/gauss.h` | F_2 Gaussian-elimination matrix. The matricial-equidistribution and t-value kernels operate on it. | — |

## Search drivers

| Type | Header | Description | Rename plan |
|---|---|---|---|
| `SearchPredicate` | `search/test.h` | Abstract predicate run on each combo by the search loop. Renamed from `Test` in Phase 4. | — |
| `EquidistributionPredicate` | `search/test.h` | Wraps an `EquidistributionMethod`. Renamed from `EquidistributionTestRunner` in Phase 4. | — |
| `CollisionFreePredicate` | `search/test.h` | Wraps `run_collision_free`. Renamed from `CollisionFreeTestRunner` in Phase 4. | — |
| `TupletsPredicate` | `search/test.h` | Wraps `run_tuplets`. Renamed from `TupletsTestRunner` in Phase 4. | — |

## YAML config

| Type | Header | Description | Rename plan |
|---|---|---|---|
| `SeekConfig` | `yaml_config/seek_config.h` | Parsed YAML config for a search run. | — |
| `SeekBuild` | `yaml_config/seek_config.h` | Materialised search (enumerator + predicates) built from a `SeekConfig`. Renamed from `BuiltSearch` in Phase 4. | — |
| `ComponentSpec` | `yaml_config/seek_config.h` | Per-slot YAML spec (source family + params + tempering chain). | — |
| `SeekTestSpec` | `yaml_config/seek_config.h` | YAML spec for one search predicate. | — |
| `SeekIterResult` | `search/search_types.h` | Per-iteration result row (ME, CF, tuplets fields). | — |
| `SeekResult` | `search/search_types.h` | Aggregated search result. | — |
| `SearchProgress` | `search/search_types.h` | Generic `{tries, elapsed_seconds}` snapshot. Emitted by `primitive_search` and `tempering_search`. | — |
| `SeekProgress` | `search/seek_search.h` | Seek-loop snapshot adding `{nbgen, nb_select, nb_me}` counters that don't apply to the other drivers. | — (kept separate; investigated 2026-05-23 — shapes differ materially) |
| `PrimitiveSearchConfig` | `search/primitive_search.h` | Config for the primitive-polynomial sub-search. | — |

## Library catalog

| Type | Header | Description | Rename plan |
|---|---|---|---|
| `Catalog` | `library/catalog.h` | Top-level library of curated generators + papers. | — |
| `CatalogGenerator` | `library/catalog.h` | A registered generator entry in the catalog. | — |
| `CatalogComponent` | `library/catalog.h` | A component of a catalog generator entry. Renamed from `library::Component` in Phase 4 (avoided collision with `regpoly::core::Component` / now `F2LinearSourcePool`). | — |
| `Paper` | `library/catalog.h` | Citation / metadata for a referenced paper. | — |
| `Author` | `library/catalog.h` | Author info attached to a `Paper`. | — |
| `TemperingStep` | `library/catalog.h` | One declared tempering step inside a `CatalogGenerator`. | — |
| `ParamValue` | `library/catalog.h` | Tagged-union value type for catalog params. | — |

## Tempering search / optimiser

| Type | Header | Description | Rename plan |
|---|---|---|---|
| `TemperingSearchConfig` | `search/tempering_search.h` | Config for the per-component tempering search. | — |
| `TemperingTryResult` | `search/tempering_search.h` | One per-try outcome. Renamed from `TemperingTryOutcome` in Phase 4. | — |
| `TemperingSearchResult` | `search/tempering_search.h` | Final result of a tempering search. | — |
| `TemperingOptimizerConfig` | `lattice/temper_optimizer.h` | Config for the tempering optimiser. | — |
| `TemperingOptimizerResult` | `lattice/temper_optimizer.h` | Per-iteration optimiser result. Renamed from `TemperingOptResult` in Phase 4. | — |
| `TemperParamLocator` | `lattice/temper_optimizer.h` | Locates mutable tempering parameters inside a chain. | — |

## Result / value types

| Type | Header | Description | Rename plan |
|---|---|---|---|
| `EquidistributionResult` | `analyses/equidistribution_runner.h` | Unified result: per-resolution ecart, cumulative SE, verified flag. Unifies the former `MatricialEquidResult` + `MeLatResult` (Phase 4). | — |
| `CollisionFreeResult` | `analyses/equidistribution_runner.h` | Per-dimension rank gap + cumulative SE-CF. | — |
| `TValueResult` | `analyses/tvalue_runner.h` | T-value profile per dimension. | — |
| `TupletsResult` | `analyses/tuplets_runner.h` | Tuplets-uniformity result. Renamed from `TupletsRunResult` in Phase 4. | — |
| `RandomParamResult` | `core/factory.h` | Result of a random-param sample (params + the constructed source). | — |
| `TestedGenerator` | `library/tested.h` | A generator plus its accumulated test results. | — |

## Internal helpers (impl-only; not part of the public API)

| Type | Where | Description | Rename plan |
|---|---|---|---|
| `UnpackedGenerator` | (anonymous namespace in `analyses/*.cpp`) | Helper struct holding `sources()` + per-source `k`. | — |
| `PhiPick` | `lattice/me_notprimitive.cpp` (anonymous) | Result of `select_phi`: selected χ factor + degree + primitivity-certified flag. | — |
| `LatticeVec`, `SimdLinVec`, `FactorEntry` | various `lattice/*.cpp` | Internal vector / factor structs. | — |
| `HaraseRankCache::Impl` | `lattice/me_harase.cpp` | Pimpl-style impl. | — |
| `Slot` | `analyses/equidistribution_method.cpp` (anonymous) | Registry slot. | — |
| `MarsaXorshiftEnumeratorBase` + 5 concrete derived classes | `generators/marsaxorshift.cpp` | Internal per-type parameter enumerators (`Type1`, `Type2`, `Type3`, `Type4`, `Type100`). | — |
| `TausworthePolyEnumerator`, `XoroshiroEnumerator`, `XoshiroEnumerator` | various `generators/*.cpp` | Internal per-family parameter enumerators. | — |

## Renames already completed

| Phase | Was | Now | Commit |
|---|---|---|---|
| 4 | `library::Component` | `CatalogComponent` | bundled |
| 4 | `NiederreiterF2Gen` | `NiederreiterF2Net` | bundled |
| 4 | `LaggedTempering` | `LaggedTemperingTrans` | bundled |
| 4 | `TemperingTryOutcome` | `TemperingTryResult` | bundled |
| 4 | `MatricialEquidResult` + `MeLatResult` | `EquidistributionResult` | bundled |
| 4 | `TupletsRunResult` | `TupletsResult` | bundled |
| 4 | `PISCache` | `HaraseRankCache` | bundled |
| 4 | `TemperOptCache` | `TemperingOptimizerCache` | bundled |
| 4 | `TemperingOptResult` | `TemperingOptimizerResult` | bundled |
| 4 | `PolVect` | `DualLatticeBasisRow` | bundled |
| 4 | `MethodRegistry` | `EquidistributionMethodRegistry` | bundled |
| 4 | `BuiltSearch` | `SeekBuild` | bundled |
| 4 | `Test` | `SearchPredicate` | bundled |
| 4 | `EquidistributionTestRunner` / `CollisionFreeTestRunner` / `TupletsTestRunner` | `*Predicate` | bundled |
| 5.1 | `Generator` | `Recurrence` | `d8311d9` |
| 5.2 | `Component` | `F2LinearSourcePool` | `6932a69` |
| 5.7 | (drop) `using Generator = Recurrence;` | — | `b36edff` |
| 5.8 | (drop) `_cpp.Component` alias + method aliases | — | `ad4a47f` |
| Tier 3 | `Params` | `ParamBag` | `f93efb5` |
| Tier 3 | `GenEnumerator` | `ParameterSpaceEnumerator` | `3eb77c1` |
| Tier 3 | `make_gen_enumerator` / `build_gen_enumerator` | `make_parameter_space_enumerator` / `build_parameter_space_enumerator` | `3eb77c1` |
| Tier 3 | header `core/gen_enumerator.{h,cpp}` | `core/parameter_space_enumerator.{h,cpp}` | `3eb77c1` |
| cleanup | `CombinedGenerator::ComponentTempering` typedef | replaced by `TemperingChain` at the API surface | `5db9b6a` |
| retirement | `CombinedGenerator` (Recurrence subclass holding J Recurrences) + `combined.{h,cpp}` | first folded into `CombinedF2LinearSource` (`9a49130`), then `CombinedF2LinearSource` was promoted to **composition**: it now inherits `ITestable` directly (no Recurrence inheritance) and the lattice / matricial kernels walk per-component state instead of treating the wrapper as a single Recurrence. `me_notprimitive.cpp`'s `step_once` / `recover_char_poly` / `output_phases` are composition-aware free functions; the combined wrapper exposes `components()` for kernels to walk. | `9a49130` + follow-up |
| widen | `CombinedF2LinearSource` components | `vector<unique_ptr<Recurrence>>` → `vector<unique_ptr<F2LinearSource>>`. DigitalNets now compose freely. Lattice/Harase/notprimitive/SIMD kernels enforce the Recurrence requirement via `cs.recurrence_components("kernel_name")`, which throws `std::invalid_argument` on a DigitalNet. | `22e8f80` |
| rename | `run_matricial_equidistribution` | `test_me_matricial`. Now every equidistribution-method kernel (`test_me_{lat,harase,notprimitive,notprimitive_simd,matricial}`) shares the `test_me_*` prefix. | `c87b9ce` (sweep) / `22e8f80` |
| migrate | kernel signatures: `vector<Recurrence*>` → `const CombinedF2LinearSource&` | Every lattice/Harase/notprimitive/SIMD/cache entry-point now takes `const CombinedF2LinearSource& cs` and internally unpacks via `cs.recurrence_components("kernel_name")`. Legacy list-form bindings, `_gen` ITestable overloads, and the `unpack_recurrences` adapter (`single_gen_adapters.cpp`) all deleted. Python `_to_cpp_gen` always returns a `CombinedF2LinearSource` (wrapping bare Generators in a J=1 wrapper). | `22e8f80` |

## Pending / deferred renames

| Item | Status | Notes |
|---|---|---|
| `SeekProgress` vs `SearchProgress` | **Resolved 2026-05-23 — keep separate.** | Shapes differ materially: `SearchProgress{tries, elapsed_seconds}` (generic, used by `primitive_search` + `tempering_search`); `SeekProgress{nbgen, nb_select, nb_me, elapsed_seconds}` (seek-loop-specific selection counters). Merging would either bloat `SearchProgress` with always-zero seek fields or introduce inheritance for no gain. |

## See also

- [`migration-2026-hierarchy-refactor.md`](migration-2026-hierarchy-refactor.md) — per-symbol migration table for the v2 hierarchy refactor.
- [`migration-2026-combination-removal.md`](migration-2026-combination-removal.md) — the prior `Combination` / `Component` Python rename.
- [`architecture.md`](architecture.md) — overall system architecture.
