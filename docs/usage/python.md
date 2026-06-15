# Using REGPOLY from Python

The `regpoly` package is a thin Python wrapper around the C++ core (`regpoly_cpp`). Use it when you want full search-loop ergonomics — DataFrames, plots, notebook integration — on top of the same algorithms the standalone CLI runs.

## Install

```bash
git clone https://github.com/pannetof/regpoly_monorepo
cd regpoly_monorepo
uv sync
```

`uv sync` installs the workspace in editable mode and builds the C++ extension via `scikit-build-core`. After it completes, the `regpoly` and `regpoly-web` console scripts are on the venv's `PATH`.

## Run a search from a YAML config

```bash
uv run regpoly shared/yaml/equidist/mt19937.yaml
```

Test config templates live under [`shared/yaml/`](https://github.com/pannetof/regpoly/tree/master/shared/yaml). The search loop iterates a `_cpp.ComboEnumerator` over generator pools, runs the configured analysis (`equidistribution`, `collision_free`, `tuplets`), and writes a tested-generator YAML for any candidate that passes the acceptance predicate.

## Use the library programmatically

```python
from regpoly import make_combined
from regpoly.core.combination_build import build_combinaison_inputs
from regpoly.library import Catalog

cat = Catalog("docs/library")
cat.load()
_, gen = cat.generator("mt19937")

gen_lists, temperings = build_combinaison_inputs(gen.components, gen.Lmax)
# Each pool is a singleton; flatten to per-slot active gens.
gens = [pool[0] for pool in gen_lists]
comb = make_combined(*gens, trans=temperings, Lmax=gen.Lmax)
```

`comb` is a `_cpp.CombinedF2LinearSource` ready to feed any test. From here:

- Run an equidistribution analysis on the first component:

    ```python
    from regpoly.analyses.pis import analyze_single_generator
    res = analyze_single_generator(comb[0])
    print(res["se"], res["gaps"])
    ```

- Iterate a Seek search:

    ```python
    from regpoly.search.seek import Seek
    seek = Seek.from_yaml("shared/yaml/equidist/mt19937.yaml")
    seek.run()
    ```

## Layer boundaries

The `regpoly` package may import `regpoly_cpp`, but the web app must go through `regpoly` only — never `_cpp` directly. `regpoly.introspection` re-exports the small set of catalog/parameter helpers the web layer needs. This is enforced by `import-linter` contracts declared in the workspace `pyproject.toml`.

## See also

- [The Python / C++ bridge](python-cpp-bridge.md) — what each layer
  owns, the wrapper-to-C++ map, the *two-`Generator` confusion*
  (`regpoly.library.Generator` is a parameter set;
  `regpoly.core.generator.Generator` is the runtime class).
- [C++ usage](cpp.md) — for users who do not want a Python runtime.
- [Web UI](web.md) — for browser-based search and result browsing.
- [Notebooks](../notebooks/index.md) — per-family demos and equidistribution exploration.
- [Architecture](../dev/architecture.md) — for a deeper look at the package layout.
