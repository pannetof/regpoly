"""
search_primitive.py — Search for full-period generators.

Reads a YAML configuration specifying the generator family, structural
parameters, and search constraints.  Iterates by generating random
search parameters, creating the generator, and testing
is_full_period().

Results are stored in::

    <generators_dir>/<Family>/<structural_params>.yaml

The output file uses the standard generator file format (family, common,
generators) so it can be passed directly to the equidistribution search.
New generators are appended to the existing file.  Duplicates are
skipped.  A .partial file protects against data loss on interruption.
"""

from __future__ import annotations

import json
import os
import sys
import time

import yaml

from regpoly.generateur import Generateur, resolve_family
from regpoly.parametric import generate_random
import regpoly._regpoly_cpp as _cpp


class PrimitiveSearch:
    """
    Search for generators whose characteristic polynomial is primitive.

    Create from a YAML config file, then call run().
    """

    def __init__(
        self,
        family: str,
        L: int,
        structural_params: dict,
        fixed_params: dict,
        max_tries: int | None,
        max_seconds: float | None,
        generators_dir: str,
    ) -> None:
        self.family = family
        self.L = L
        self.structural_params = structural_params
        self.fixed_params = fixed_params
        self.max_tries = max_tries
        self.max_seconds = max_seconds
        self.generators_dir = generators_dir
        self.output_file = self._build_output_path()

    @classmethod
    def from_yaml(cls, config_file: str) -> "PrimitiveSearch":
        with open(config_file) as f:
            config = yaml.safe_load(f)

        search = config["search"]
        family = search["family"]
        L = search["L"]
        limit = search.get("limit", {})
        max_tries = limit.get("max_tries")
        max_seconds = limit.get("max_seconds")
        generators_dir = search.get("generators_dir", "yaml/generators")

        structural_params = config.get("structural_params",
                                       config.get("structural", {}))
        fixed_params = config.get("fixed_params",
                                   config.get("search_params", {}))

        return cls(
            family=family,
            L=L,
            structural_params=structural_params,
            fixed_params=fixed_params,
            max_tries=max_tries,
            max_seconds=max_seconds,
            generators_dir=generators_dir,
        )

    def run(self) -> list[dict]:
        family = self.family
        L = self.L

        # Recover partial results from a previous interrupted run
        partial_file = self.output_file + ".partial"
        if os.path.exists(partial_file):
            print(f"Found partial results from an interrupted run: "
                  f"{partial_file}")
            new_entries = _read_partial(partial_file)
            self._merge_into_output(new_entries)
            os.remove(partial_file)
            print(f"Recovered {len(new_entries)} generators into "
                  f"{self.output_file}")
            print()

        # Build the fixed params: structural + fixed_params with explicit values
        fixed = dict(self.structural_params)
        for key, val in self.fixed_params.items():
            if val is not None:
                fixed[key] = val

        # Identify which params will be randomized
        resolved = resolve_family(family, fixed)
        specs = _cpp.get_gen_param_specs(resolved)
        randomized_names = [
            s["name"] for s in specs
            if s["name"] not in fixed
            and not s["structural"]
            and not s["has_default"]
            and s["rand_type"]
            and s["rand_type"] != "none"
        ]

        # Load existing generators to deduplicate
        existing = self._load_existing_generators()

        # Header — create one throwaway generator to determine k
        try:
            probe_params = self._randomize(fixed, specs)
            probe_gen = Generateur.create(family, L, **probe_params)
            k_str = str(probe_gen.k)
        except Exception:
            k_str = "?"
        print(f"Searching for full-period {family} generators")
        print(f"  L = {L},  k = {k_str}")
        fixed_search = {k: v for k, v in self.fixed_params.items()
                        if v is not None}
        if fixed_search:
            print(f"  Fixed search params: {fixed_search}")
        print(f"  Randomized: {randomized_names}")
        if self.max_tries:
            print(f"  Max tries: {self.max_tries}")
        if self.max_seconds:
            print(f"  Max seconds: {self.max_seconds}")
        print(f"  Output: {self.output_file}")
        if existing:
            print(f"  Existing generators in file: {len(existing)}")
        print()
        sys.stdout.flush()

        os.makedirs(os.path.dirname(self.output_file), exist_ok=True)
        partial_fd = open(partial_file, "w")

        results = []
        tries = 0
        found = 0
        t_start = time.time()
        use_progress = self.max_tries is not None and sys.stderr.isatty()

        try:
            while True:
                if self.max_tries and tries >= self.max_tries:
                    break
                if self.max_seconds and (time.time() - t_start) >= self.max_seconds:
                    break

                tries += 1
                full_params = self._randomize(fixed, specs)

                try:
                    gen = Generateur.create(family, L, **full_params)
                except Exception as e:
                    if not use_progress:
                        print(f"  [try {tries}] Error: {e}", file=sys.stderr)
                    continue

                if gen.is_full_period():
                    # Extract only the search params (not structural)
                    entry = {k: _yaml_safe(v)
                             for k, v in full_params.items()
                             if k not in self.structural_params}

                    # Deduplicate
                    entry_key = _entry_key(entry)
                    if entry_key in existing:
                        continue

                    existing.add(entry_key)
                    results.append(entry)
                    found += 1

                    partial_fd.write(json.dumps(entry) + "\n")
                    partial_fd.flush()

                    if use_progress:
                        _clear_progress()
                    elapsed = time.time() - t_start
                    print(f"  [{tries:>8d}] Found #{found:>4d}  "
                          f"k={gen.k}  ({elapsed:.1f}s)")
                    sys.stdout.flush()

                if use_progress and tries % 100 == 0:
                    _show_progress(tries, self.max_tries, found,
                                   time.time() - t_start)

        except KeyboardInterrupt:
            if use_progress:
                _clear_progress()
            print(f"\n  Interrupted after {tries} tries.")

        partial_fd.close()
        elapsed = time.time() - t_start
        if use_progress:
            _clear_progress()

        # Merge new results into the output file
        if results:
            self._merge_into_output(results)
        if os.path.exists(partial_file):
            os.remove(partial_file)

        total = len(existing)
        print()
        print(f"Search complete: {tries} tries, "
              f"{found} new found, {total} total in file, {elapsed:.1f}s")
        print(f"Results: {self.output_file}")

        return results

    # -- Output path -------------------------------------------------------

    def _build_output_path(self) -> str:
        """Compute results/<Family>/<structural_params>.yaml."""
        resolved = resolve_family(self.family, self.structural_params)
        parts = []
        for key, val in sorted(self.structural_params.items()):
            parts.append(f"{key}{val}")
        filename = "_".join(parts) + ".yaml" if parts else "default.yaml"
        return os.path.join(self.generators_dir, resolved, filename)

    # -- File I/O ----------------------------------------------------------

    def _load_existing_generators(self) -> set[str]:
        """Load existing generator entries as a set of keys for dedup."""
        if not os.path.exists(self.output_file):
            return set()
        with open(self.output_file) as f:
            data = yaml.safe_load(f)
        if not data or "generators" not in data:
            return set()
        return {_entry_key(g) for g in data["generators"]}

    def _merge_into_output(self, new_entries: list[dict]) -> None:
        """Append new generator entries to the output file."""
        os.makedirs(os.path.dirname(self.output_file), exist_ok=True)

        if os.path.exists(self.output_file):
            with open(self.output_file) as f:
                data = yaml.safe_load(f)
        else:
            data = None

        if not data or "generators" not in data:
            resolved = resolve_family(self.family, self.structural_params)
            data = {
                "family": resolved,
                "common": dict(self.structural_params),
                "generators": [],
            }

        data["generators"].extend(new_entries)

        with open(self.output_file, "w") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)

    # -- Helpers -----------------------------------------------------------

    @staticmethod
    def _randomize(fixed: dict, specs: list[dict]) -> dict:
        """Fill in missing non-structural params with random values."""
        full = dict(fixed)
        for spec in specs:
            name = spec["name"]
            if name in full or spec["has_default"] or spec["structural"]:
                continue
            rt = spec["rand_type"]
            if not rt or rt == "none":
                continue
            full[name] = generate_random(spec, full)
        return full



# ═══════════════════════════════════════════════════════════════════════════
# Partial file helpers
# ═══════════════════════════════════════════════════════════════════════════

def _read_partial(path: str) -> list[dict]:
    entries = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                entries.append(json.loads(line))
    return entries


def _entry_key(entry: dict) -> str:
    """Stable string key for deduplication."""
    return json.dumps(entry, sort_keys=True)


# ═══════════════════════════════════════════════════════════════════════════
# Progress bar
# ═══════════════════════════════════════════════════════════════════════════

_BAR_WIDTH = 30


def _show_progress(tries: int, max_tries: int, found: int,
                   elapsed: float) -> None:
    pct = tries / max_tries
    filled = int(_BAR_WIDTH * pct)
    bar = "#" * filled + "-" * (_BAR_WIDTH - filled)
    rate = tries / elapsed if elapsed > 0 else 0
    sys.stderr.write(
        f"\r  [{bar}] {pct:5.1%}  "
        f"{tries}/{max_tries}  "
        f"found: {found}  "
        f"({rate:.0f}/s)"
    )
    sys.stderr.flush()


def _clear_progress() -> None:
    sys.stderr.write("\r" + " " * 80 + "\r")
    sys.stderr.flush()


# ═══════════════════════════════════════════════════════════════════════════
# YAML helpers
# ═══════════════════════════════════════════════════════════════════════════

def _yaml_safe(v):
    """Convert numpy/large ints to plain Python ints for YAML serialization."""
    if isinstance(v, int) and not isinstance(v, bool):
        return int(v)
    if isinstance(v, list):
        return [_yaml_safe(x) for x in v]
    return v
