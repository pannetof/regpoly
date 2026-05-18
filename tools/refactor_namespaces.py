#!/usr/bin/env python3
"""Wrap regpoly-cpp public symbols into the `regpoly::*` namespace tree.

Phase A — rename existing sub-namespaces:
    regpoly_internal    -> regpoly::internal
    regpoly_random      -> regpoly::random
    regpoly_catalog     -> regpoly::library
    regpoly_yaml_config -> regpoly::yaml_config

Phase B — wrap every other public header (and its matching .cpp impl) in
    namespace regpoly::core { ... }

Phase C — in every other .cpp file that consumes regpoly-cpp symbols,
    add `using namespace regpoly::core;` (+ other using-directives as needed)
    after the last #include line.

The script is idempotent. Re-running it does no harm.
"""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CPP_PKG = ROOT / "packages" / "regpoly-cpp"
SRC_INCLUDE = CPP_PKG / "src" / "include"
SRC = CPP_PKG / "src"
TESTS = CPP_PKG / "tests"

# Headers that already live in a regpoly_X namespace — rename to regpoly::X.
EXISTING_NS_RENAMES = {
    "regpoly_internal":    "regpoly::internal",
    "regpoly_random":      "regpoly::random",
    "regpoly_catalog":     "regpoly::library",
    "regpoly_yaml_config": "regpoly::yaml_config",
}

# Which header path -> which target namespace.
HEADER_TO_NS: dict[Path, str] = {}
for h in SRC_INCLUDE.rglob("*.h"):
    text = h.read_text()
    matched = None
    for old, new in EXISTING_NS_RENAMES.items():
        if re.search(rf"^\s*namespace\s+{re.escape(old)}\s*\{{", text, re.M):
            matched = new
            break
    HEADER_TO_NS[h] = matched if matched is not None else "regpoly::core"


def rename_existing_namespaces() -> None:
    """Phase A: pure string-rewrite across every .h, .cpp, .py."""
    for ext in ("*.h", "*.cpp", "*.hpp"):
        for f in (ROOT / "packages").rglob(ext):
            if "/build/" in str(f) or "/third_party/" in str(f):
                continue
            text = orig = f.read_text()
            for old, new in EXISTING_NS_RENAMES.items():
                # `namespace regpoly_X {` -> `namespace regpoly::X {`
                text = re.sub(
                    rf"\bnamespace\s+{re.escape(old)}\b",
                    f"namespace {new}",
                    text,
                )
                # `regpoly_X::` -> `regpoly::X::`
                text = re.sub(
                    rf"\b{re.escape(old)}::",
                    f"{new}::",
                    text,
                )
                # `// namespace regpoly_X` (closing brace comments) -> renamed
                text = re.sub(
                    rf"//\s*namespace\s+{re.escape(old)}\b",
                    f"// namespace {new}",
                    text,
                )
            if text != orig:
                f.write_text(text)


# ── Phase B: wrap headers + matching .cpp in `namespace regpoly::core { ... }`

# Skip wrapping these — they already have a renamed namespace via Phase A.
def _has_target_namespace(text: str, target: str) -> bool:
    """Already wrapped in the target namespace?"""
    # `namespace regpoly::core` literal, or `namespace regpoly { namespace core`
    return bool(re.search(rf"^\s*namespace\s+{re.escape(target)}\s*\{{", text, re.M))


def wrap_in_namespace(file_path: Path, ns: str) -> bool:
    """Wrap a file's body in `namespace <ns> { ... }`. Returns True if changed."""
    text = file_path.read_text()
    if _has_target_namespace(text, ns):
        return False

    lines = text.splitlines(keepends=True)

    # Find insertion point: after the last #include / #pragma / preamble
    # comment / blank line, before the first code declaration.
    insert_at = 0
    in_preamble = True
    for i, line in enumerate(lines):
        s = line.lstrip()
        if (s.startswith("//")
                or s.startswith("/*")
                or s.startswith("*")
                or s.startswith("#")
                or s.strip() == ""):
            if in_preamble:
                insert_at = i + 1
        else:
            in_preamble = False
            break

    # Special-case: pybind11 forward decls like
    #   namespace pybind11 { class module_; }
    # Skip over them when looking for our insertion point.
    while insert_at < len(lines):
        s = lines[insert_at].lstrip()
        if s.startswith("namespace pybind11"):
            insert_at += 1
            while insert_at < len(lines) and lines[insert_at].strip() == "":
                insert_at += 1
            continue
        break

    open_block = f"namespace {ns} {{\n\n"
    close_block = f"\n}}  // namespace {ns}\n"

    # Insert open block at insert_at.
    new_lines = lines[:insert_at] + [open_block] + lines[insert_at:]

    # Append close block at EOF (after stripping trailing whitespace).
    while new_lines and new_lines[-1].strip() == "":
        new_lines.pop()
    new_lines.append(close_block)

    file_path.write_text("".join(new_lines))
    return True


def wrap_headers_and_impls() -> int:
    """Phase B. Returns count of files modified."""
    changed = 0
    for header in sorted(HEADER_TO_NS):
        ns = HEADER_TO_NS[header]
        if ns != "regpoly::core":
            continue  # Phase A already handled these
        if wrap_in_namespace(header, ns):
            changed += 1

        # Find matching .cpp file. The headers live at src/include/<area>/<name>.h;
        # impls at src/<area>/<name>.cpp.
        area = header.parent.name
        name = header.stem
        impl = SRC / area / f"{name}.cpp"
        if impl.exists() and wrap_in_namespace(impl, ns):
            changed += 1
    return changed


# ── Phase C: add using-namespace directives to consumer .cpp files

# Map include-file basename → set of namespaces it introduces.
# Anything not in this map is presumed to introduce `regpoly::core`.
# Headers that split content across namespaces list both.
NS_BY_INCLUDE: dict[str, set[str]] = {
    "primitivity.h":     {"regpoly::core", "regpoly::internal"},
    "random_samplers.h": {"regpoly::core", "regpoly::random"},
    "catalog.h":         {"regpoly::core", "regpoly::library"},
    "seek_config.h":     {"regpoly::core", "regpoly::yaml_config"},
}


def _ns_for(basename: str) -> set[str]:
    """Return the set of namespaces introduced by including `basename`."""
    return NS_BY_INCLUDE.get(basename, {"regpoly::core"})


REGPOLY_AREAS = {
    "algebra", "analyses", "core", "generators", "lattice",
    "library", "search", "transforms", "yaml_config",
}


def _namespaces_introduced_by_includes(text: str) -> set[str]:
    """Scan #include lines and return the set of regpoly:: namespaces brought in.

    Only counts includes that look like ours (quoted, with one of our area
    prefixes OR a bare header basename matching one of our public headers).
    System includes like <gtest/gtest.h> or <yaml-cpp/yaml.h> are ignored.
    """
    ns: set[str] = set()
    for line in text.splitlines():
        s = line.lstrip()
        if not s.startswith("#include"):
            continue
        # Capture both "..." and <...> paths.
        m = re.search(r'#include\s+["<]([^">]+)[">]', s)
        if not m:
            continue
        path = m.group(1)
        parts = path.split("/")
        basename = parts[-1]
        # Match `<regpoly/foo.h>`: installed-style; basename is what matters.
        if parts[0] == "regpoly":
            if basename in NS_BY_INCLUDE:
                ns |= NS_BY_INCLUDE[basename]
            elif any((SRC_INCLUDE / a / basename).exists() for a in REGPOLY_AREAS):
                ns.add("regpoly::core")
            continue
        # Match `"area/foo.h"`.
        if len(parts) >= 2 and parts[0] in REGPOLY_AREAS:
            ns |= _ns_for(basename)
            continue
        # Match bare quoted basename `"foo.h"` only if it's one of ours.
        # System headers like <gtest/gtest.h> are excluded here.
        if s.startswith("#include \""):
            if basename in NS_BY_INCLUDE:
                ns |= NS_BY_INCLUDE[basename]
            elif any((SRC_INCLUDE / a / basename).exists() for a in REGPOLY_AREAS):
                ns.add("regpoly::core")
    return ns


def add_using_directives(file_path: Path) -> bool:
    """Insert include-aware `using namespace ...;` directives after the last
    include line. Skips if any of them is already present.
    """
    text = file_path.read_text()
    namespaces = _namespaces_introduced_by_includes(text)
    if not namespaces:
        return False
    # Idempotent: ensure every required directive is present. If all are
    # already there, skip.
    missing = [
        ns for ns in namespaces
        if f"using namespace {ns};" not in text
    ]
    if not missing:
        return False
    namespaces = set(missing)

    lines = text.splitlines(keepends=True)
    last_include = -1
    for i, line in enumerate(lines):
        if line.lstrip().startswith("#include"):
            last_include = i
    if last_include < 0:
        return False  # nothing to attach to

    block = "\n" + "".join(
        f"using namespace {ns};\n" for ns in sorted(namespaces)
    ) + "\n"
    new_lines = lines[:last_include + 1] + [block] + lines[last_include + 1:]
    file_path.write_text("".join(new_lines))
    return True


def add_using_to_consumers() -> int:
    """Phase C. Add using-directives to EVERY .cpp file under src/ and tests/.

    The directives are harmless in files that don't reference the namespace,
    and necessary in any file that references symbols from another area
    (e.g. `regpoly::library::Catalog::load()` calls into `regpoly::core`).
    Even impls of `regpoly::core::*` benefit from `using namespace
    regpoly::library;` etc. if they cross-reference.
    """
    changed = 0
    # Skip impls whose body is wrapped INSIDE `namespace regpoly::core { ... }`
    # by Phase B — for those, `using namespace regpoly::core;` would be a
    # redundancy (and gcc accepts it as no-op). But we still WANT
    # `using namespace regpoly::library;` etc. inside them. So: always add.
    for cpp in list((CPP_PKG / "src").rglob("*.cpp")) + list((CPP_PKG / "tests").rglob("*.cpp")):
        if "/build/" in str(cpp) or "/third_party/" in str(cpp):
            continue
        if add_using_directives(cpp):
            changed += 1
    return changed


def main() -> None:
    print("Phase A — renaming existing namespaces …")
    rename_existing_namespaces()

    # Refresh HEADER_TO_NS since Phase A renamed some files in place.
    HEADER_TO_NS.clear()
    for h in SRC_INCLUDE.rglob("*.h"):
        text = h.read_text()
        matched = None
        for new in ("regpoly::internal", "regpoly::random",
                    "regpoly::library", "regpoly::yaml_config"):
            if re.search(rf"^\s*namespace\s+{re.escape(new)}\s*\{{", text, re.M):
                matched = new
                break
        HEADER_TO_NS[h] = matched if matched is not None else "regpoly::core"

    print("Phase B — wrapping headers + impls into `regpoly::core` …")
    n_b = wrap_headers_and_impls()
    print(f"  modified {n_b} files")

    print("Phase C — adding using-namespace to consumer .cpp files …")
    n_c = add_using_to_consumers()
    print(f"  modified {n_c} files")


if __name__ == "__main__":
    main()
