# API documentation style guide

```{note}
Placeholder — the full style guide (NumPy-style Python docstrings,
catalog.h Doxygen precedent, cross-reference roles, `@ingroup`
discipline) lands as part of Step 10 of the Sphinx framework
migration. Until then, the conventions in active use are:

* **Python:** NumPy-style sections (`Parameters`, `Returns`, `Raises`,
  `See Also`, `Examples`). Cross-reference C++ counterparts with the
  `:cpp:class:` / `:cpp:func:` roles inside a `See Also` block. See
  any of the six rewritten modules (`regpoly.library`,
  `regpoly.core.generator`, `regpoly.analyses.pis`,
  `regpoly.introspection`, `regpoly_legacy.reader`,
  `regpoly_legacy.seek_factory`) for examples.
* **C++:** see `packages/regpoly-cpp/src/include/library/catalog.h`
  for the per-symbol Doxygen block template — required tags are
  `@brief`, `@param`, `@return`, `@throws`, `@ingroup`; required prose
  is a class-level long description preserving invariants.
```
