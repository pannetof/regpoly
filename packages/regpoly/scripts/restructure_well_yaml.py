#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""One-shot WELL YAML restructure: per-slot positional → T0..T7 map.

Converts the legacy
    matrices:
      - {type: 3, pi: [-25, 0, 0], pu: [0, 0, 0]}
      - {type: 1, pi: [0, 0, 0], pu: [0, 0, 0]}
      - ...
form to the paper-aligned
    matrices:
      T0: {M: 3, t: -25}
      T1: {M: 1}
      ...
form.

Mirrors the per-Mi positional decoding done by `legacy_reader.cpp` —
M0/M1 take no args; M2/M3 use pi[0] as t; M4 uses pu[0] as a; M5 uses
(pi[0], pu[0]) as (t, b); M6 uses (pi[0], pu[0..2]) decoded back to
paper (q, t, s, a) via the same single-bit / single-zero-bit pattern
recognition as the C++ helpers `decode_d_s_mask` / `decode_test_mask`.

Drops `mat_types_version: 2` if present (the structural shape is now
the discriminant — no version key needed).

Usage:
    uv run --with ruamel.yaml python packages/regpoly/scripts/restructure_well_yaml.py

Targets the two committed WELL YAMLs by default. Pass paths to override:
    ... restructure_well_yaml.py path/to/file.yaml [more.yaml ...]
"""

from __future__ import annotations

import sys
from pathlib import Path

from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap

DEFAULT_TARGETS = [
    Path("shared/yaml/generators/WELLRNG/p31_r624_w32.yaml"),
    Path("shared/yaml/generators/WELLRNG/p15_r1391_w32.yaml"),
]

M32 = 0xFFFFFFFF


def decode_d_s_mask(mask: int, ctx: str) -> int:
    """Return s such that mask == ~(1 << (31 - s)). Raises if mask is
    not a single-zero-bit pattern."""
    inv = (~int(mask)) & M32
    if inv == 0 or inv & (inv - 1):
        raise ValueError(
            f"{ctx}: cannot decode WELL paper M6 d_s mask 0x{int(mask):08x}: "
            f"expected exactly one zero bit, got mask with {bin(inv).count('1')} "
            f"flipped bits."
        )
    lsb = (inv & -inv).bit_length() - 1
    return 31 - lsb


def decode_test_mask(mask: int, ctx: str) -> int:
    """Return t such that mask == (1 << (31 - t))."""
    m = int(mask) & M32
    if m == 0 or m & (m - 1):
        raise ValueError(
            f"{ctx}: cannot decode WELL paper M6 test mask 0x{int(mask):08x}: "
            f"expected exactly one set bit."
        )
    lsb = (m & -m).bit_length() - 1
    return 31 - lsb


def slot_to_paper_form(slot_idx: int, entry: dict, ctx: str) -> dict:
    """Convert one legacy {type, pi, pu} slot to the paper-form dict."""
    Mi = int(entry["type"])
    pi = list(entry.get("pi", [0, 0, 0]))
    pu = list(entry.get("pu", [0, 0, 0]))
    while len(pi) < 3:
        pi.append(0)
    while len(pu) < 3:
        pu.append(0)
    out: dict = {"M": Mi}
    slot_ctx = f"{ctx} slot T{slot_idx}"
    if Mi in (0, 1):
        pass  # no args
    elif Mi in (2, 3):
        out["t"] = int(pi[0])
    elif Mi == 4:
        out["a"] = int(pu[0]) & M32
    elif Mi == 5:
        out["t"] = int(pi[0])
        out["b"] = int(pu[0]) & M32
    elif Mi == 6:
        out["q"] = int(pi[0])
        out["t"] = decode_test_mask(pu[2], slot_ctx)
        out["s"] = decode_d_s_mask(pu[1], slot_ctx)
        out["a"] = int(pu[0]) & M32
    else:
        raise ValueError(f"{slot_ctx}: unknown M-class {Mi}")
    return out


def restructure_file(path: Path, yaml: YAML) -> None:
    if not path.exists():
        print(f"  {path}: missing, skipping", file=sys.stderr)
        return
    data = yaml.load(path.read_text())

    # Detect already-restructured (matrices is a CommentedMap, not list).
    sample = data.get("generators", [{}])
    if sample:
        first = sample[0] if isinstance(sample, list) else None
        if first and isinstance(first.get("matrices"), dict):
            print(f"  {path}: already in T0..T7 map form, skipping")
            return

    # Drop the old mat_types_version marker if present (shape is the
    # discriminant now).
    if "mat_types_version" in data:
        del data["mat_types_version"]

    n_gens = 0
    for gen in data.get("generators", []):
        legacy = gen.get("matrices")
        if not isinstance(legacy, list):
            continue
        new_matrices = CommentedMap()
        for j, slot in enumerate(legacy):
            if not isinstance(slot, dict) or "type" not in slot:
                raise ValueError(
                    f"{path}: generator #{n_gens} slot #{j} is not a "
                    f"legacy {{type, pi, pu}} dict"
                )
            new_matrices[f"T{j}"] = slot_to_paper_form(
                j, slot, f"{path} gen#{n_gens}"
            )
        gen["matrices"] = new_matrices
        n_gens += 1

    yaml.dump(data, path.open("w"))
    print(f"  {path}: restructured {n_gens} generators")


def main() -> None:
    yaml = YAML()
    yaml.preserve_quotes = True
    yaml.indent(mapping=2, sequence=4, offset=2)
    yaml.width = 4096

    repo_root = Path.cwd()
    targets = (
        [Path(p) for p in sys.argv[1:]]
        if len(sys.argv) > 1
        else [repo_root / t for t in DEFAULT_TARGETS]
    )
    for t in targets:
        print(f"Restructuring {t.relative_to(repo_root) if t.is_relative_to(repo_root) else t}")
        restructure_file(t, yaml)


if __name__ == "__main__":
    main()
