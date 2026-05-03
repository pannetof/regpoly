"""
tested_generator.py — Save and load tested generator YAML files.

A tested generator file describes one generator (possibly multi-component)
with its tempering chain and the results of the test that validated it.

File format (single component)::

    generator:
      family: TGFSR
      L: 32
      w: 32
      r: 3
      m: 1
      a: 2571403067

    tempering:
      - type: tempMK
        w: 32
        eta: 7
        mu: 15
        b: 2636928640
        c: 2611425280

    results:
      equidistribution:
        method: lattice
        se: 42
        ecart: [0, 0, 0, 2, ...]

File format (multi-component)::

    components:
      - generator:
          family: TGFSR
          L: 32
          w: 32
          r: 3
          m: 1
          a: 2571403067
        tempering:
          - type: tempMK
            ...
      - generator:
          family: TGFSR
          L: 31
          w: 31
          r: 5
          m: 2
          a: 3456789012

    results:
      equidistribution:
        se: 12
        ecart: [0, 0, 0, 0, ...]

Directory structure::

    yaml/testedgenerators/<Family>/<structural>.<testname>.<serial>.yaml
"""

from __future__ import annotations

import os
import re

import yaml

from regpoly.core.generator import Generator
from regpoly.core.transformation import Transformation
from regpoly.core.combination import Combination


def structural_filename(gen: Generator) -> str:
    """Build the structural-params part of the filename (e.g. 'r3_w32')."""
    sp = gen.structural_params()
    if not sp:
        return f"k{gen.k}"
    return "_".join(f"{k}{v}" for k, v in sorted(sp.items()))


def save_tested_generator(
    output_dir: str,
    test_name: str,
    comb: Combination,
    results: dict,
) -> str:
    """
    Save a tested generator to a YAML file.

    Parameters
    ----------
    output_dir : base directory (e.g. "yaml/testedgenerators")
    test_name  : e.g. "equidist", "collision_free"
    comb       : the Combination that was tested
    results    : dict of test results to include in the file

    Returns the path of the created file.
    """
    # Determine family and structural params from the first component
    gen0 = comb[0]
    family = gen0.type_name or "Unknown"
    struct_str = structural_filename(gen0)

    # Build directory and find next serial number
    family_dir = os.path.join(output_dir, family)
    os.makedirs(family_dir, exist_ok=True)
    serial = _next_serial(family_dir, struct_str, test_name)
    filename = f"{struct_str}.{test_name}.{serial:04d}.yaml"
    filepath = os.path.join(family_dir, filename)

    # Build the YAML content
    if comb.J == 1:
        data = _single_component_data(comb, 0)
    else:
        data = _multi_component_data(comb)

    data["results"] = results

    with open(filepath, "w") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)

    return filepath


def load_tested_generator(filepath: str) -> tuple[Combination, dict]:
    """
    Load a tested generator YAML file.

    Returns (Combination, results_dict).
    The Combination is ready to use (generator created, tempering applied).
    """
    with open(filepath) as f:
        data = yaml.safe_load(f)

    results = data.get("results", {})

    if "components" in data:
        comb = _load_multi(data)
    else:
        comb = _load_single(data)

    return comb, results


# ═══════════════════════════════════════════════════════════════════════════
# Save helpers
# ═══════════════════════════════════════════════════════════════════════════

def _single_component_data(comb: Combination, j: int) -> dict:
    gen = comb[j]
    comp = comb.components[j]

    gen_data = {"family": gen.type_name, "L": gen.L}
    if gen.params:
        gen_data.update({k: _yaml_safe(v)
                         for k, v in gen.params.items()})

    data: dict = {"generator": gen_data}

    if comp.trans:
        data["tempering"] = [_trans_to_dict(t) for t in comp.trans]

    return data


def _multi_component_data(comb: Combination) -> dict:
    components = []
    for j in range(comb.J):
        gen = comb[j]
        comp = comb.components[j]

        gen_data = {"family": gen.type_name, "L": gen.L}
        if gen.params:
            gen_data.update({k: _yaml_safe(v)
                             for k, v in gen.params.items()})

        entry: dict = {"generator": gen_data}
        if comp.trans:
            entry["tempering"] = [_trans_to_dict(t) for t in comp.trans]
        components.append(entry)

    return {"components": components}


def _trans_to_dict(t: Transformation) -> dict:
    d = {"type": t._type_name}
    d.update({k: _yaml_safe(v) for k, v in t._params.items()})
    return d


# ═══════════════════════════════════════════════════════════════════════════
# Load helpers
# ═══════════════════════════════════════════════════════════════════════════

def _load_single(data: dict) -> Combination:
    gen_data = data["generator"]
    family = gen_data["family"]
    L = gen_data["L"]
    params = {k: v for k, v in gen_data.items() if k not in ("family", "L")}

    gen = Generator.create(family, L, **params)

    comb = Combination(J=1, Lmax=L)
    comb.components[0].add_gen(gen)

    for t_data in data.get("tempering", []):
        trans_type = t_data["type"]
        t_params = {k: v for k, v in t_data.items() if k != "type"}
        t = Transformation.create(trans_type, **t_params)
        comb.components[0].add_trans(t)

    comb.reset()
    return comb


def _load_multi(data: dict) -> Combination:
    comp_list = data["components"]
    J = len(comp_list)
    Lmax = 0

    # First pass: find Lmax
    for entry in comp_list:
        L = entry["generator"]["L"]
        if L > Lmax:
            Lmax = L

    comb = Combination(J=J, Lmax=Lmax)

    for j, entry in enumerate(comp_list):
        gen_data = entry["generator"]
        family = gen_data["family"]
        L = gen_data["L"]
        params = {k: v for k, v in gen_data.items()
                  if k not in ("family", "L")}

        gen = Generator.create(family, L, **params)
        comb.components[j].add_gen(gen)

        for t_data in entry.get("tempering", []):
            trans_type = t_data["type"]
            t_params = {k: v for k, v in t_data.items() if k != "type"}
            t = Transformation.create(trans_type, **t_params)
            comb.components[j].add_trans(t)

    comb.reset()
    return comb


# ═══════════════════════════════════════════════════════════════════════════
# Utilities
# ═══════════════════════════════════════════════════════════════════════════

def _next_serial(directory: str, struct_str: str, test_name: str) -> int:
    """Find the next available serial number for the given prefix."""
    pattern = re.compile(
        rf'^{re.escape(struct_str)}\.{re.escape(test_name)}\.(\d+)\.yaml$'
    )
    max_serial = 0
    if os.path.isdir(directory):
        for name in os.listdir(directory):
            m = pattern.match(name)
            if m:
                max_serial = max(max_serial, int(m.group(1)))
    return max_serial + 1


def _yaml_safe(v):
    if isinstance(v, int) and not isinstance(v, bool):
        return int(v)
    if isinstance(v, list):
        return [_yaml_safe(x) for x in v]
    return v
