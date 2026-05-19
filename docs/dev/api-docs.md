# API documentation style guide

Conventions for writing and rendering REGPOLY's API reference. This
page describes how to document a new symbol so it slots into the
existing site without introducing warnings or breaking the Python/C++
cross-reference graph.

The framework is **Sphinx** orchestrating **autodoc + Napoleon** for
Python and **Doxygen + Breathe + Exhale** for C++. The conventions
here are tuned to that stack; they may not transfer wholesale to a
different framework.

## When to add documentation

Every new public symbol ships with documentation in the same commit:

- **Python:** every member of a module's `__all__` (or, when there is
  no `__all__`, every non-underscore module-level symbol).
- **C++:** every class, struct, enum, free function, public method
  declared in a header under `packages/regpoly-cpp/src/include/`.
  Doxygen's `EXTRACT_PRIVATE` is `NO`, so private members are
  implicitly excluded.

Code that doesn't yet have documentation is a backlog item, not a
licence to add more undocumented code.

## Python — NumPy style with Sphinx domain cross-refs

`docs/conf.py` configures Napoleon with `napoleon_numpy_docstring = True`
and `napoleon_google_docstring = False`. Google-style docstrings will
render but with degraded section headers — write **NumPy style** for
every new docstring.

### Required sections

| Section      | When                                           |
|--------------|------------------------------------------------|
| Summary line | Always — one short imperative sentence.        |
| Long description | When the summary line isn't self-contained. |
| `Parameters` | Whenever the function takes arguments.         |
| `Returns`    | Whenever the function returns a non-`None` value. |
| `Raises`     | Whenever the function may raise; one entry per exception type. |

### Recommended sections

| Section      | When                                                   |
|--------------|--------------------------------------------------------|
| `See Also`   | When a C++ counterpart exists (cross-ref it with `:cpp:class:`), or when a closely-related Python symbol is worth pointing at. |
| `Examples`   | On non-trivial public entries. Doctest-shaped. |
| `Notes`      | Math derivations, performance characteristics, references to a paper or theory page. |
| `Attributes` | On dataclass-style wrappers (`Catalog`, `Paper`, …).   |
| `Yields`     | On iterator / generator methods.                       |
| `Warnings`   | Caveats users will hit.                                |
| `Warns`      | `warnings.warn(...)` calls the function may emit.      |

### Skeleton

```python
def analyze_single_generator(gen, *, L=32, method="matricial"):
    """Compute the equidistribution profile of one F2-linear generator.

    Thin Python shim over the C++ matricial-equidistribution runner.
    The output dict is the canonical shape downstream code (the web
    analysis worker, the `regpoly` CLI's `show` mode) expects.

    Parameters
    ----------
    gen : regpoly.core.generator.Generator
        A constructed generator instance.
    L : int, optional
        Output resolution in bits (default 32).
    method : {'matricial', 'lattice'}, optional
        Equidistribution method.

    Returns
    -------
    dict
        Keys ``'se'``, ``'gaps'``, ``'k'``, ``'L'``, ``'elapsed'``,
        ``'char_poly_int'``, ``'hamming_weight'``.

    Raises
    ------
    ValueError
        If ``L`` is outside ``[1, gen.w]``.
    RuntimeError
        If the C++ driver fails to converge.

    See Also
    --------
    :cpp:func:`regpoly::core::run_matricial_equidistribution` :
        The C++ implementation this wraps.

    Examples
    --------
    >>> from regpoly.library import Catalog
    >>> from regpoly.core.generator import Generator
    >>> from regpoly.analyses.pis import analyze_single_generator
    >>> cat = Catalog("docs/library")
    >>> cat.load()  # doctest: +SKIP
    >>> _, params = cat.generator("mt19937")  # doctest: +SKIP
    >>> gen = Generator.create(params.family, L=32, **params.params)  # doctest: +SKIP
    >>> res = analyze_single_generator(gen)  # doctest: +SKIP
    >>> res["se"]  # doctest: +SKIP
    0
    """
```

### Cross-reference roles

In docstring prose and `See Also` entries:

| Target                          | Role                                       |
|---------------------------------|--------------------------------------------|
| Python class                    | `` :class:`regpoly.x.Y` ``                 |
| Python method                   | `` :meth:`regpoly.x.Y.method` ``           |
| Python function                 | `` :func:`regpoly.x.foo` ``                |
| Python module                   | `` :mod:`regpoly.x` ``                     |
| C++ class                       | `` :cpp:class:`regpoly::core::Generator` `` |
| C++ free function               | `` :cpp:func:`regpoly::core::run_matricial_equidistribution` `` |
| C++ struct member               | `` :cpp:member:`regpoly::core::SeekIterResult::me_se` `` |
| External (NumPy, stdlib)        | The `intersphinx_mapping` in `conf.py`
covers `python`, `numpy`, `pathlib`; just write the type name. |

Always pick the shortest role that resolves uniquely (`:class:` over
`:py:class:` — both work, the shorter form keeps docstrings legible).

### Two-`Generator` confusion

Two classes share the name `Generator`. Disambiguate **explicitly** in
every wrapper-class docstring:

- :class:`regpoly.core.generator.Generator` — the **runtime** class.
  A constructed, parameterised, seedable bit source.
- :class:`regpoly.library.Generator` — a **parameter set**
  (re-exported from `regpoly_cpp.CatalogGenerator`). A record inside
  a paper's catalog YAML.

When in doubt, write the fully-qualified name. See
[The Python / C++ bridge](../usage/python-cpp-bridge.md) for the
wider explanation.

## C++ — Doxygen with Sphinx domain cross-refs

`docs/Doxyfile` runs Doxygen over `packages/regpoly-cpp/src/include/`,
emits XML, and feeds it to Breathe + Exhale. The conventions below
match the `library/catalog.h` and `core/generator.h` exemplars
verbatim.

### File-level block

Every header opens with an `@file` block immediately after the SPDX
licence line:

```cpp
/**
 * @file <filename>.h
 * @brief <one-line summary>.
 * @defgroup <area> <Title>           // FIRST header in the subdirectory only
 *
 * <Optional 1-3 paragraph prose: what this header defines, how it
 * fits the larger picture, cross-header invariants.>
 *
 * @ingroup <area>
 */
```

- The **first header** alphabetically in each subdirectory carries
  the `@defgroup <area> <Title>` directive. Subsequent headers in the
  same subdirectory carry `@ingroup <area>` only.
- Current groups (one per `src/include/<subdir>/`):
  `algebra` · `analyses` · `core` · `generators` · `lattice` ·
  `library` · `search` · `transforms` · `yaml_config`.
- Headers that live in a non-`regpoly::core` namespace (e.g.
  `core/random_samplers.h` in `regpoly::random`,
  `algebra/primitivity.h` partially in `regpoly::internal`) still
  tag `@ingroup core` — the group is a *subdirectory* organising
  device, the Exhale "Namespace Hierarchy" handles the namespace
  view separately.

### Per-symbol block

Every public class / struct / enum / free function / non-trivial
method:

```cpp
/**
 * @brief One-line summary. Imperative mood. Ends with a period.
 *
 * <Long description. Required on every class/struct: preserve any
 * existing prose-level invariants from previous // comments — do
 * NOT strip those during a mechanical pass.>
 *
 * @param name1  What it is, including unit / range / precondition.
 * @param name2  ...
 * @return       What is returned. SINGULAR @return — Doxygen accepts
 *               @returns as an alias but the project standard is
 *               @return. Drop entirely for void.
 * @throws std::runtime_error  When this is thrown.
 * @pre   Precondition the caller must establish.
 * @post  Postcondition this method guarantees.
 *
 * @code{.cpp}
 *   // REQUIRED on every user-constructible class and every
 *   // user-callable free function. Keep examples short (3-6 lines)
 *   // and illustrative — they are NOT compile-checked.
 * @endcode
 *
 * @see :cpp:class:`regpoly::core::Other`  — related C++ symbol.
 * @see :py:class:`regpoly.core.x.Other`   — Python counterpart.
 *
 * @ingroup <area>
 */
```

Mandatory:
- `@brief` on every documented symbol.
- `@param <name>` for every parameter (each on its own line, matching
  the signature exactly — a `@param` for a non-existent parameter
  emits a Sphinx warning).
- `@return` for every non-void function. Singular form.
- `@throws <Type>` for every exception that may propagate out.
- `@ingroup <area>` on every **class / struct only**. **Do NOT** add
  `@ingroup` to free functions, enums, type aliases, or macros —
  doing so breaks Breathe's `doxygenfunction` / `doxygenenum`
  resolution (the symbol gets pulled out of namespace XML into
  group XML, where Breathe can't find it). The file-level
  `@ingroup` plus the per-class `@ingroup` is enough to populate
  the group tree pages.
- **Class-level long description** preserving existing prose
  invariants. A mechanical pass that only writes `@brief`s and strips
  the prose is rejected at review.
- `@code{.cpp}` block on every class a downstream user would
  construct directly — every `*Gen` family, `Combination`,
  `Catalog`, `TemperingOptimizer`, etc.

Trivial accessors (one-line getters):

```cpp
/** @brief State size in bits. */
int k() const { return k_; }
```

### Doxygen prose pitfalls

Symbols and patterns that confuse Doxygen / Pygments / docutils:

| Symptom                                     | Fix                                  |
|---------------------------------------------|--------------------------------------|
| Indented multi-line block reinterpreted as docutils block-quote / definition-list | Wrap in `@verbatim` / `@endverbatim`. |
| `**` / `++` / `*` in prose parsed as Markdown emphasis | Backtick them or wrap in `@verbatim`. |
| `<` / `>` in template syntax inside prose   | Write in a `@code` block or escape as `\<`. |
| `/*foo=*/` inline-C-comment inside `@code{.cpp}` confuses Pygments | Use line comments `// foo` instead. |
| Doxygen "nested comment" error              | Look for `/*` inside an existing `/** */` block — usually a `lattice/*.h` glob mention in prose; reword to "the `lattice/` headers". |

### Cross-reference roles inside Doxygen blocks

Breathe re-parses Doxygen output through Sphinx, so Sphinx domain
roles work inside `@see` / `@brief` prose:

```cpp
/// @see :py:class:`regpoly.core.generator.Generator`  — Python wrapper.
/// @see :cpp:class:`regpoly::core::Combination`       — sibling C++ type.
```

Use `:py:func:` for free functions, `:py:meth:` for methods,
`:py:mod:` for modules. The full mapping from C++ to Python lives in
[The Python / C++ bridge](../usage/python-cpp-bridge.md#cross-references-in-the-docs).

## Building + previewing the docs

```bash
# One-time setup:
uv sync --group docs
sudo apt install -y doxygen graphviz

# One-shot build:
uv run --group docs sphinx-build -b html docs site

# Live-reload preview:
uv run --group docs sphinx-autobuild docs site

# Strict (warnings-as-errors — the CI gate):
uv run --group docs sphinx-build -W -b html docs site
```

When `sphinx-build -W` lands as the CI gate, every new warning blocks
merge — write documentation with that in mind.

### Known-benign warnings

The build currently emits ~62 warnings that are **not** regressions:

- **Exhale "Duplicate C++ declaration"** on nested structs (every
  `*Gen::MiEntry`, `Catalog::PapersFilter`, etc.). Exhale registers
  these both as members of the parent class and as standalone pages.
  Suppressed via `cpp.duplicate_declaration` in `conf.py`, but the
  underlying warnings still print.
- **docutils "line-length-limit"** on the Exhale-generated
  `class_view_hierarchy.rst.include` (one warning, one long line).

Any other warning category in a fresh build is a regression to fix
before submitting.

## Reviewing a documentation change

Mechanical checklist for reviewers (and for AI-assisted backfills):

- [ ] Every new public symbol carries a docstring / Doxygen block.
- [ ] Section style is NumPy (Python) / Doxygen `@param @return @throws @ingroup` (C++).
- [ ] `@return` is singular, not `@returns`.
- [ ] `@ingroup <area>` is on every class/struct, **off** every free
      function and enum.
- [ ] Class-level prose preserves invariants from the previous
      `//` comment (if any) — no mechanical strip.
- [ ] `@code{.cpp}` example present on every user-constructible class
      and every user-callable free function.
- [ ] Cross-references use the Sphinx domain roles (`:cpp:class:`,
      `:py:class:`, etc.) — not Markdown links.
- [ ] `sphinx-build` warning count is at-or-below baseline.
- [ ] `sphinx-build -b doctest` (when wired) passes.

## See also

- [The Python / C++ bridge](../usage/python-cpp-bridge.md) — the
  layered split + the full wrapper-to-C++ cross-reference table.
- [Building](building.md) — install, run tests, build docs.
- [Contributing](contributing.md) — branching, PR conventions, CI
  lanes.
- `packages/regpoly-cpp/src/include/library/catalog.h` — the gold
  standard exemplar for C++ Doxygen blocks.
- `packages/regpoly/src/regpoly/library/__init__.py` — the gold
  standard exemplar for Python NumPy-style docstrings.
