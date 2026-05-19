# The Python / C++ bridge

REGPOLY is a C++ library with a thin Python wrapper. Most users see only
the Python side, but every algorithm runs in C++ — the wrapper exists to
make the same routines pleasant to drive from a notebook, plot results
in pandas, and share fixtures with the FastAPI web app. This page
explains where the seam is, why it is shaped the way it is, and how the
two layers cross-reference each other in the documentation.

If you only want to *use* REGPOLY, skip to
[Python usage](python.md) or [C++ usage](cpp.md). Read this page when
you want to know which side owns what, or when you are reading a
docstring that says *"thin wrapper around `regpoly::core::X`"* and
wondering whether the wrapper does anything you should care about.

## The shape of the split

| Concern                                          | Lives in C++ | Lives in Python |
|--------------------------------------------------|:------------:|:---------------:|
| GF(2) algebra: `BitVect`, `BitMatrix`, BM, Gauss | ✅           |                 |
| Generator families: MT, WELL, MELG, SFMT, …      | ✅           |                 |
| Tempering transformations                         | ✅           |                 |
| Combined-generator iterator                       | ✅           |                 |
| Equidistribution / collision-free / tuplets       | ✅           |                 |
| Catalog (paper-centric YAML) reader               | ✅           |                 |
| Search drivers (Seek, PrimitiveSearch, …)         | ✅           |                 |
| YAML config parsing                               | ✅           |                 |
| `regpoly-cli` standalone binary                   | ✅           |                 |
| `regpoly` CLI argv handling + ergonomic output    |              | ✅              |
| Result decoration (DataFrames, plots)             |              | ✅              |
| FastAPI web shell + task queue                    |              | ✅              |
| Pre-v2 `.dat` parameter reader (`regpoly-legacy`) |              | ✅              |

Two layered import contracts (enforced by `import-linter` in CI) keep
dependencies one-way:

```text
regpoly-web  →  regpoly  →  regpoly-cpp        (catalog / YAML / web path)
regpoly-legacy → regpoly  →  regpoly-cpp        (pre-v2 .dat path)
```

`regpoly-web` may not import `_regpoly_cpp` directly; it goes through
`regpoly`. `regpoly-legacy` and `regpoly-web` never import each other.

## What the Python wrapper actually does

When a docstring says **"thin wrapper around `regpoly::core::X`"**, the
wrapper is doing one or more of these:

1. **Translating Python types into C++ inputs.** A `pathlib.Path` becomes
   `std::string`; a Python `list[int]` becomes `std::vector<int>`; a
   missing optional argument is filled by a randomiser declared in the
   C++ `ParamSpec`.
2. **Packaging C++ outputs into idiomatic Python.** A
   `MatricialEquidResult` struct becomes the flat dict
   `{"se", "gaps", "elapsed", "char_poly_int", "hamming_weight", "k", "L"}`
   that {func}`regpoly.analyses.pis.analyze_single_generator` returns,
   and that the web `results` table persists.
3. **Holding onto a C++ handle.** Every Python wrapper class
   (`Generator`, `Combination`, `Catalog`, …) has a private `_cpp_*`
   attribute pointing at the underlying C++ object. Lifetime is managed
   by pybind11 — when the Python wrapper is garbage-collected, so is
   the C++ handle.

Nothing in the wrapper layer reimplements an algorithm. If you see a
loop over output bits, equidistribution dimensions, or tempering masks
in a Python file, it is **either** glue around a C++ call **or** legacy
code that has not yet been migrated.

## The two-`Generator` confusion

REGPOLY has two classes both called `Generator`:

- {class}`regpoly.core.generator.Generator` — the **runtime** class. A
  constructed, parameterised, seedable bit source you can iterate. This
  is the one search loops and analyses operate on.
- {class}`regpoly.library.Generator` — a **parameter set**, re-exported
  from `regpoly_cpp.CatalogGenerator`. It is a record inside a paper's
  catalog YAML — a name, a family, and a dict of parameters. It is not
  iterable; you build a runtime `Generator` from it.

The relationship:

```python
from regpoly.library import Catalog
from regpoly.core.generator import Generator as RuntimeGenerator

cat = Catalog("docs/library"); cat.load()
_, params = cat.generator("mt19937")        # params is a CatalogGenerator
gen = RuntimeGenerator.create(              # gen is a runtime Generator
    params.family, L=32, **params.params
)
```

The C++ side has the same shape (and the same names, in different
namespaces): `regpoly::core::Generator` is the runtime class,
`regpoly::library::CatalogGenerator` is the parameter record.

## Cross-references in the docs

Every Python wrapper that has a C++ counterpart links to it via a
`See Also` block with the `:cpp:class:` Sphinx role:

```
See Also
--------
:cpp:class:`regpoly::core::Generator` : the C++ implementation this wraps.
```

The reverse — `@see :py:class:` in a C++ Doxygen block — is added on
the C++ side as part of the per-header Doxygen backfill. When a class
has both a Python and a C++ page, clicking either side's "See Also"
hops to the other. The wrapper map:

| Python                                              | C++                                  |
|-----------------------------------------------------|--------------------------------------|
| {class}`regpoly.core.generator.Generator`           | `regpoly::core::Generator`           |
| {class}`regpoly.core.combination.Combination`       | `regpoly::core::Combination`         |
| {class}`regpoly.core.transformation.Transformation` | `regpoly::core::Transformation`      |
| {class}`regpoly.core.component.Component`           | `regpoly::core::Component`           |
| {class}`regpoly.library.Catalog`                    | `regpoly::library::Catalog`          |
| {class}`regpoly.library.Paper`                      | `regpoly::library::Paper`            |
| {class}`regpoly.library.Generator` (CatalogGenerator) | `regpoly::library::CatalogGenerator` |
| {class}`regpoly.search.seek.Seek`                   | `regpoly::core::Seek`                |
| {class}`regpoly.search.tempering_optimizer.TemperingOptimizer` | `regpoly::core::TemperingOptimizer` |

A handful of public Python helpers (e.g. `regpoly.introspection.get_gen_param_specs`)
wrap a C++ free function rather than a class — those use `:cpp:func:`.

## Choosing a layer

| If you want to …                                     | Use                                     |
|------------------------------------------------------|-----------------------------------------|
| run a search from a YAML config on the command line | `regpoly` (Python CLI) or `regpoly-cli` |
| script a search from a notebook                      | Python — {class}`regpoly.search.seek.Seek` |
| embed the search loop in a C++ application           | C++ — `regpoly::core::Seek` (see [C++ usage](cpp.md)) |
| browse a catalog through a web UI                    | `regpoly-web` (FastAPI; Python only)    |
| read a legacy `.dat` parameter file                  | `regpoly-legacy` (Python only)          |
| publish a tested-generator entry into the catalog    | either: `regpoly publish` or `regpoly-cli publish` |

The C++ layer is feature-complete for the catalog / YAML / search /
publish path. Python adds ergonomic decoration around the same calls,
plus the `.dat` reader and the web UI — neither of which has a C++
counterpart.

## See also

- [Python usage](python.md) — programmatic entry points.
- [C++ usage](cpp.md) — pure-C++ build + runnable snippet.
- [Architecture](../dev/architecture.md) — package layout + dependency
  arrows in more detail.
