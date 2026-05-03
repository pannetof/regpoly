# Architecture

> **Status**: Phase 0 placeholder. Phase 6 will lift content from the active redesign plan at `~/.claude/plans/i-want-to-have-rippling-starlight.md`.

REGPOLY is a `uv` workspace of three packages with strict one-way dependencies:

```
regpoly-web  →  regpoly  →  regpoly-cpp
   FastAPI       Python      C++ core
```

`import-linter` enforces this in CI. The C++ core carries every algorithm, every YAML schema, the catalog reader/writer, the legacy `.dat` reader, and a standalone `regpoly-cli` binary so a user working exclusively in C++ has full feature parity with a Python user.
