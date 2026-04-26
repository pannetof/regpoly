#!/usr/bin/env python3
"""Materialise per-family `*_mttoolbox_reference.py` modules.

Reads catalog entries from ``docs/library/*_params.yaml``, invokes the
matching MTToolBox `calc_equidist`-style binary per entry, parses
``k(v)/d(v)`` lines, keeps only ``v = 1..MAX_V`` (default 5), and writes
``tests/<family>_mttoolbox_reference.py`` modules consumed by
``tests/test_mttoolbox_crosscheck.py``.

Run manually after editing a catalog or after rebuilding a binary:

    cd MinimalCode/cpp_regpoly
    python tests/_gen_mttoolbox_reference.py            # all families
    python tests/_gen_mttoolbox_reference.py sfmt mt    # subset

The script is committed alongside the catalogs so the materialisation
is reproducible.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_LIB = REPO_ROOT / "docs" / "library"
TESTS_DIR = REPO_ROOT / "tests"
MTTOOLBOX_SAMPLES = Path(
    "/home/frpan/projets/claude_projects/regpoly/MTToolBox/samples"
)

MAX_V = 5  # cross-check granularity (per user spec)
DEFAULT_TIMEOUT_S = int(os.environ.get("MTTOOLBOX_REF_TIMEOUT_S", "1800"))
# Skip mexps > this in the default pass; pass --slow to include them.
SLOW_MEXP = 20000

# Two output formats observed:
#   sfmtdc / dSFMTdc / RMT:  "k(NN) = NNN  d(NN) = NNN"
#   MTDC mt1 / XORSHIFT-*:   "k(NN):NNN  d(NN):NNN"
LINE_RE = re.compile(
    r"k\(\s*(\d+)\)\s*[=:]\s*(\d+)\s+d\(\s*\1\)\s*[=:]\s*(\d+)"
)

# Map family-yaml file → output-module name and Python dict identifier.
# Output modules end in `_d5.py` to keep them distinct from any pre-existing
# full-d(1..32) reference modules (e.g. sfmt_mttoolbox_reference.py).
CATALOG_TO_REF = {
    "sfmt_params.yaml": ("sfmt_mttoolbox_d5.py", "SFMT_MTTOOLBOX_D5"),
    "mt_params.yaml": ("mt_mttoolbox_d5.py", "MT_MTTOOLBOX_D5"),
    "dsfmt_params.yaml": ("dsfmt_mttoolbox_d5.py", "DSFMT_MTTOOLBOX_D5"),
    "xorshift_params.yaml": ("xorshift_mttoolbox_d5.py", "XORSHIFT_MTTOOLBOX_D5"),
    "tinymt_params.yaml": ("tinymt_mttoolbox_d5.py", "TINYMT_MTTOOLBOX_D5"),
    "rmt_params.yaml": ("rmt_mttoolbox_d5.py", "RMT_MTTOOLBOX_D5"),
}


def _binary_path(catalog: dict, entry: dict) -> Path:
    """Resolve the MTToolBox binary, allowing per-entry override."""
    rel = entry.get("mttoolbox_binary") or catalog.get("mttoolbox_binary")
    if rel is None:
        raise RuntimeError(
            f"no mttoolbox_binary in catalog or entry {entry.get('id')!r}"
        )
    path = MTTOOLBOX_SAMPLES / rel
    if not path.exists():
        raise FileNotFoundError(
            f"binary {path} not built — see tests/_mttoolbox_build.md"
        )
    return path


def _invoke(binary: Path, argv: list[str]) -> str:
    """Run binary with optional argv (from yaml) and -v for verbose output."""
    cmd = [str(binary), "-v", *argv] if argv else [str(binary)]
    # mt1 / xorshift-* don't accept flags. Try -v first, fall back without.
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=DEFAULT_TIMEOUT_S,
            check=False,
        )
        out = result.stdout
        if LINE_RE.search(out):
            return out
        # Fall through and retry without -v.
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"timeout running {' '.join(cmd)}")

    cmd = [str(binary), *argv] if argv else [str(binary)]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=DEFAULT_TIMEOUT_S,
        check=False,
    )
    return result.stdout


def _parse_dv(text: str, max_v: int) -> list[int]:
    """Extract d(v) for v=1..max_v from the FIRST k(v)/d(v) block.

    Some calc_equidist binaries print multiple blocks (e.g. 32-bit /
    64-bit / 128-bit accuracy for SFMT).  We treat the first contiguous
    sequence starting at v=1 as the canonical 32-bit-accuracy block
    (matches existing tests/sfmt_mttoolbox_reference.py convention).
    """
    found: dict[int, int] = {}
    last_v = 0
    for m in LINE_RE.finditer(text):
        v = int(m.group(1))
        d = int(m.group(3))
        # Stop reading once v drops back to 1 after we've moved past v=1
        # (= start of the next accuracy block).
        if v <= last_v and 1 in found:
            break
        if 1 <= v <= max_v:
            found[v] = d
        last_v = v
    missing = [v for v in range(1, max_v + 1) if v not in found]
    if missing:
        raise RuntimeError(
            f"missing d(v) for v={missing}; raw output:\n{text[:500]}"
        )
    return [found[v] for v in range(1, max_v + 1)]


def _emit_module(out_path: Path, dict_name: str, family: str,
                 entries: list[tuple[str, int, list[int]]]) -> None:
    """Write a Python module with one dict mapping id → d(1..MAX_V)."""
    lines = [
        '"""',
        f"{family} d(1..{MAX_V}) reference values from MTToolBox",
        f"calc_equidist (auto-generated by tests/_gen_mttoolbox_reference.py).",
        "",
        "Edit docs/library/*_params.yaml and re-run the generator instead",
        "of editing this file directly.",
        '"""',
        "",
        f"{dict_name} = {{",
    ]
    for gen_id, mexp, dv in entries:
        formatted = ", ".join(str(d) for d in dv)
        lines.append(
            f"    {gen_id!r}: [{formatted}],   # mexp={mexp}"
        )
    lines.append("}")
    lines.append("")
    out_path.write_text("\n".join(lines))


def process_catalog(catalog_file: Path, include_slow: bool,
                    refresh: bool = False) -> None:
    """Process one *_params.yaml file → matching reference module."""
    with catalog_file.open() as f:
        catalog = yaml.safe_load(f)
    family = catalog["family"]
    if catalog.get("deferred"):
        print(f"[skip] {catalog_file.name}: {catalog.get('deferred_reason')}")
        return

    out_module, dict_name = CATALOG_TO_REF[catalog_file.name]
    out_path = TESTS_DIR / out_module

    # Preserve previously-captured entries when re-running (so a fast
    # pass doesn't drop entries the slow pass already populated).  Look
    # up each entry's mexp from the catalog up-front so partial flushes
    # don't write mexp=0 placeholders.
    catalog_mexps: dict[str, int] = {}
    for entry in catalog.get("generators", []):
        catalog_mexps[entry["id"]] = entry["mexp"]
    entries: dict[str, tuple[int, list[int]]] = {}
    if out_path.exists():
        # Best-effort read of the previous module to keep slow-lane entries.
        try:
            ns: dict = {}
            exec(out_path.read_text(), ns)
            prev = ns.get(dict_name) or {}
            for k, v in prev.items():
                # If the catalog lists this id, recover its mexp;
                # otherwise drop the stale entry rather than write mexp=0.
                m = catalog_mexps.get(k)
                if m is not None:
                    entries[k] = (m, v)
        except Exception:
            pass

    def _flush() -> None:
        """Write the current entries dict to disk so partial captures survive
        a later timeout/crash.  Called after every successful capture."""
        if not entries:
            return
        sorted_entries = [
            (gid, mexp if mexp is not None else 0, dv)
            for gid, (mexp, dv) in entries.items()
        ]
        sorted_entries.sort(key=lambda t: (t[1], t[0]))
        _emit_module(out_path, dict_name, family, sorted_entries)

    threshold = catalog.get("slow_threshold_mexp", SLOW_MEXP)
    for entry in catalog.get("generators", []):
        gen_id = entry["id"]
        mexp = entry["mexp"]
        if mexp > threshold and not include_slow:
            if gen_id in entries:
                # Keep the previously captured value, just refresh mexp.
                entries[gen_id] = (mexp, entries[gen_id][1])
                print(f"[keep] {family} {gen_id} (mexp={mexp}, "
                      f"slow lane — pass --slow to refresh)")
            else:
                print(f"[defer] {family} {gen_id} (mexp={mexp}, "
                      f"slow lane — pass --slow to capture)")
            continue
        if gen_id in entries and not refresh:
            # Already captured (loaded from prior file).  Skip re-running
            # so a resumed slow-lane pass doesn't re-pay multi-hour costs
            # for entries already on disk.  Pass --refresh to override.
            print(f"[skip] {family} {gen_id} (mexp={mexp}, already captured)")
            continue
        argv = entry.get("mttoolbox_argv") or []
        binary = _binary_path(catalog, entry)
        print(f"[run] {family} {gen_id} (mexp={mexp}) ...", flush=True)
        try:
            text = _invoke(binary, argv)
            dv = _parse_dv(text, MAX_V)
        except Exception as exc:
            print(f"[fail] {family} {gen_id}: {exc}", flush=True)
            _flush()  # persist anything captured before this failure
            continue
        print(f"      d(1..{MAX_V}) = {dv}", flush=True)
        entries[gen_id] = (mexp, dv)
        _flush()  # persist after every success — survives later timeouts

    if not entries:
        print(f"[skip] {family}: no generator entries")
        return

    print(f"[ok] wrote {out_path}", flush=True)


def main(argv: list[str]) -> int:
    args = list(argv[1:])
    include_slow = False
    refresh = False
    if "--slow" in args:
        include_slow = True
        args.remove("--slow")
    if "--refresh" in args:
        refresh = True
        args.remove("--refresh")
    selected = set(args)
    for catalog_file in sorted(DOCS_LIB.glob("*_params.yaml")):
        prefix = catalog_file.stem.replace("_params", "")
        if selected and prefix not in selected:
            continue
        process_catalog(catalog_file, include_slow=include_slow, refresh=refresh)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
