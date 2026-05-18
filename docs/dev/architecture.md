# Architecture

> **Status**: Phase 0 placeholder. Phase 6 will lift content from the active redesign plan at `~/.claude/plans/i-want-to-have-rippling-starlight.md`.

REGPOLY is a `uv` workspace of four packages with strict one-way dependencies:

```
regpoly-web      →  regpoly  →  regpoly-cpp
   FastAPI           Python      C++ core
regpoly-legacy   →  regpoly  →  regpoly-cpp
   pure-Python       (same arrow — `regpoly-legacy` is an optional add-on)
```

`import-linter` enforces both layered arrows in CI. The C++ core carries every algorithm, the YAML search schema, the catalog reader/writer, and a standalone `regpoly-cli` binary so a user working exclusively in C++ has feature parity with a Python user for the *catalog-driven* path. The historical pre-v2 `.dat` parameter format is **not** supported by the C++ core; reading it requires the optional pure-Python `regpoly-legacy` add-on, which parses `.dat` files in Python and constructs `Generator` objects via the existing C++ factory.
