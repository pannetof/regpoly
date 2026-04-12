"""
seek.py — Seek class: runs the combined-generator search.

Supports two input formats:
  - YAML: single search config file (search.*.yaml)
  - Legacy: nb_comp + test_file + gen_data_files (C-compatible)
"""

from __future__ import annotations

import os
import random
import socket
import sys
import time

from regpoly.generateur import Generateur
from regpoly.transformation import Transformation
from regpoly.combinaison import Combinaison
from regpoly.legacy_reader import LegacyReader
from regpoly.analyses.equidistribution_test import (
    EquidistributionTest,
    METHOD_MATRICIAL,
    METHOD_DUALLATTICE,
    METHOD_NOTHING,
)
from regpoly.analyses.collision_free_test import CollisionFreeTest
from regpoly.analyses.tuplets_test import TupletsTest
from regpoly.analyses.tuplets_results import _MAX_TYPE

_SEP      = "\n\n" + "+" * 104
_EQ66     = "=" * 66
_EQ66dash = "-" * 66


class Seek:
    """
    Search for combined generators with good equidistribution properties.

    Two ways to create:
        Seek.from_yaml("search.config.yaml")
        Seek.from_legacy(nb_comp, test_file, gen_data_files)
    """

    def __init__(self) -> None:
        self._comb: Combinaison = None
        self._tests: list = []
        self._nbtries: int = 1

    # -- Constructors ------------------------------------------------------

    @classmethod
    def from_yaml(cls, config_file: str) -> "Seek":
        """Build a Seek from a YAML search config file."""
        import yaml

        s = cls()
        base_dir = os.path.dirname(os.path.abspath(config_file))

        with open(config_file) as f:
            config = yaml.safe_load(f)

        search = config.get("search", {})
        seeds = search.get("seed", [-1, 0])
        Lmax = search.get("Lmax", 32)
        s._nbtries = search.get("nbtries", 1)

        # Seed RNG
        seed1, seed2 = seeds
        Seed1, Seed2, seed = _compute_seeds(seed1, seed2)
        random.seed(seed)

        # Components
        components_cfg = config.get("components", [])
        nb_comp = len(components_cfg)
        gen_lists = []
        temperings = []
        prev_gen_list = None

        for comp_cfg in components_cfg:
            gen_cfg = comp_cfg.get("generators", {})
            if gen_cfg == "same":
                gen_lists.append(prev_gen_list)
            elif isinstance(gen_cfg, dict) and "file" in gen_cfg:
                path = _resolve_path(gen_cfg["file"], base_dir)
                gen_list = Generateur.from_yaml(path, Lmax)
                gen_lists.append(gen_list)
                prev_gen_list = gen_list
            elif isinstance(gen_cfg, dict) and "legacy_file" in gen_cfg:
                path = _resolve_path(gen_cfg["legacy_file"], base_dir)
                gen_list = LegacyReader.read_generators(path, Lmax)
                gen_lists.append(gen_list)
                prev_gen_list = gen_list
            elif isinstance(gen_cfg, dict) and "family" in gen_cfg:
                gen_list = _build_inline_generators(gen_cfg, Lmax)
                gen_lists.append(gen_list)
                prev_gen_list = gen_list
            else:
                raise ValueError(f"Invalid generators spec: {gen_cfg}")

            trans_cfg = comp_cfg.get("tempering", [])
            trans_list = _parse_tempering(trans_cfg)
            temperings.append(trans_list)

        # Tests
        tests_cfg = config.get("tests", [])
        if isinstance(tests_cfg, dict) and "file" in tests_cfg:
            from regpoly.analyses.test_base import AbstractTest
            path = _resolve_path(tests_cfg["file"], base_dir)
            s._tests = AbstractTest.from_yaml(path, Lmax)
        else:
            s._tests = _parse_tests(tests_cfg, Lmax)

        # Build Combinaison
        has_tempering = any(len(tl) > 0 for tl in temperings)
        if not has_tempering:
            s._nbtries = 1

        print(_format_header(nb_comp, Seed1, Seed2, temperings))

        s._comb = Combinaison(J=nb_comp, Lmax=Lmax)
        for j, (gen_list, trans_list) in enumerate(zip(gen_lists, temperings)):
            for t in trans_list:
                s._comb.components[j].add_trans(t)
            if gen_list is not None:
                if j > 0 and gen_list is gen_lists[j - 1]:
                    s._comb.components[j].gens = s._comb.components[j - 1].gens
                else:
                    for gen in gen_list:
                        s._comb.components[j].add_gen(gen)

        print(_format_search_summary(s._tests, s._nbtries, has_tempering, gen_lists))
        sys.stdout.flush()

        return s

    @classmethod
    def from_legacy(cls, nb_comp: int, test_file: str, gen_data_files: list) -> "Seek":
        """Build a Seek from legacy C-format files."""
        s = cls()
        s._comb, tests, s._nbtries = _read_legacy(nb_comp, test_file, gen_data_files)
        s._tests = tests
        return s

    # -- Run ---------------------------------------------------------------

    def run(self) -> None:
        """Execute the search, printing results to stdout."""
        comb = self._comb
        nbtries = self._nbtries
        tests = self._tests

        if not comb.reset():
            print("First combined generator not found or invalid Combinaison")
            return

        nbgen = nbsel = nbME = nbCF = 0
        no_try = 1
        t_start = time.time()

        while True:
            comb.update_trans()
            nbgen += 1

            # Run tests in order
            passed = True
            me_results = None
            tup_results = None

            for test in tests:
                if isinstance(test, EquidistributionTest):
                    me_results = test.run(comb)
                    if not me_results.is_presque_me():
                        passed = False
                        break

                elif isinstance(test, TupletsTest):
                    tup_results = test.run(comb)
                    if not tup_results.is_ok():
                        passed = False
                        break

                elif isinstance(test, CollisionFreeTest):
                    test.run(comb, me_results=me_results)

            if passed and me_results is not None:
                print(_format_current_comb(comb))
                if me_results.is_me():
                    msg = me_results.display()
                    if msg:
                        print(msg)
                    nbME += 1
                else:
                    print("\n  Dimension gaps for every resolution", end="")
                    table_str, _ = me_results.display_table(comb, 'l')
                    print(table_str)
                if tup_results is not None:
                    msg = tup_results.display()
                    if msg:
                        print(msg)
                nbsel += 1
                print(_SEP)

            sys.stdout.flush()

            if no_try < nbtries:
                no_try += 1
            else:
                try:
                    next(comb)
                except StopIteration:
                    break
                no_try = 1

        elapsed = time.time() - t_start
        print(_format_summary(nbgen, nbME, nbCF, nbsel, elapsed))


# ═══════════════════════════════════════════════════════════════════════════
# Private helpers
# ═══════════════════════════════════════════════════════════════════════════

def _compute_seeds(seed1, seed2):
    if seed1 == -1:
        seed = int(time.time())
        Seed1 = float(seed)
        Seed2 = Seed1 + 111119903.0
    else:
        Seed1 = float(seed1)
        Seed2 = float(seed2)
        seed = seed1 * (2**32) + seed2
    return Seed1, Seed2, seed


def _resolve_path(path: str, base_dir: str) -> str:
    if os.path.isabs(path):
        return path
    return os.path.join(base_dir, path)


def _build_inline_generators(gen_cfg: dict, L: int) -> list:
    family = gen_cfg["family"]
    family_params = {k: v for k, v in gen_cfg.items()
                     if k not in ("family", "common", "generators")}
    common = {**family_params, **gen_cfg.get("common", {})}
    generators = []
    for entry in gen_cfg["generators"]:
        params = {**common, **entry}
        generators.append(Generateur.create(family, params, L))
    return generators


def _parse_tempering(trans_cfg) -> list:
    if not isinstance(trans_cfg, list):
        return []
    trans_list = []
    for entry in trans_cfg:
        params = {}
        random_specs = {}
        for k, v in entry.items():
            if k == "type":
                continue
            if isinstance(v, dict) and "random" in v:
                random_specs[k] = v
                params[k] = 0
            else:
                params[k] = v
        w_original = params.get("w", 0)
        t = Transformation(entry["type"], params, w_original)
        t._random_specs = random_specs
        trans_list.append(t)
    return trans_list


def _parse_tests(tests_cfg: list, Lmax: int) -> list:
    tests = []
    for test_cfg in tests_cfg:
        test_type = test_cfg["type"]
        if test_type == "equidistribution":
            tests.append(EquidistributionTest._from_params(test_cfg, Lmax))
        elif test_type == "collision_free":
            tests.append(CollisionFreeTest._from_params(test_cfg, Lmax))
        elif test_type == "tuplets":
            tests.append(TupletsTest._from_params(test_cfg, Lmax))
        else:
            raise ValueError(f"Unknown test type: {test_type}")
    return tests


def _read_legacy(nb_comp, test_file, gen_data_files):
    """Read legacy C-format files and return (comb, tests, nbtries)."""
    with open(test_file, encoding="latin-1") as f:
        lines = [ln.split('#')[0].strip() for ln in f]
    lines = [ln for ln in lines if ln]
    i = 0

    parts = lines[i].split(); i += 1
    seed1 = int(parts[0])
    seed2 = int(parts[1]) if seed1 != -1 else 0
    seeds = [seed1, seed2]

    Lmax = int(lines[i]); i += 1

    temperings = []
    mk_opt = False
    has_tempering = False
    for _ in range(nb_comp):
        parts = lines[i].split(); i += 1
        deftrans = int(parts[0])
        if deftrans:
            trans_list, mk_j = LegacyReader.read_transformations(parts[1])
            temperings.append(trans_list)
            mk_opt |= mk_j
            has_tempering = True
        else:
            temperings.append([])

    nbtries = int(lines[i]); i += 1
    if not has_tempering:
        nbtries = 1

    parts = lines[i].split(); i += 1
    mse_val = int(parts[0])
    if mse_val == -1:
        meverif = False
        mse = sys.maxsize
        method = METHOD_NOTHING
    else:
        meverif = True
        mse = mse_val
        method_str = parts[1].lower()
        if method_str == "matrix":
            method = METHOD_MATRICIAL
        elif method_str == "lattice":
            method = METHOD_DUALLATTICE
        else:
            raise ValueError(f"Unknown method '{method_str}' in {test_file}")

    delta = [10_000_000] * (Lmax + 1)
    delta_lines = []
    nb_lines = int(lines[i]); i += 1
    for _ in range(nb_lines):
        parts = lines[i].split(); i += 1
        jmin, jmax, ec = int(parts[0]), int(parts[1]), int(parts[2])
        delta_lines.append((jmin, jmax, ec))
        for j in range(jmin, jmax + 1):
            delta[j] = ec

    parts = lines[i].split(); i += 1
    tupverif = int(parts[0])
    if tupverif:
        d = int(parts[1])
        s_list = [0] + [int(parts[2 + j]) for j in range(d)]
        mDD = float(parts[2 + d])
        tuplets_test = TupletsTest(tupletsverif=True, d=d, s=s_list, mDD=mDD, testtype=_MAX_TYPE)
    else:
        tuplets_test = TupletsTest(tupletsverif=False)

    gen_lists = []
    old_gen_list = None
    for path in gen_data_files:
        if path.lower() == "same":
            gen_lists.append(old_gen_list)
        else:
            old_gen_list = LegacyReader.read_generators(path, Lmax)
            gen_lists.append(old_gen_list)

    # Seed RNG
    Seed1, Seed2, seed = _compute_seeds(seed1, seed2)
    random.seed(seed)

    print(_format_header(nb_comp, Seed1, Seed2, temperings))

    comb = Combinaison.CreateFromFiles(gen_lists, Lmax, temperings)

    if has_tempering:
        print(f"Number of tries per combined generator : {nbtries}")
    if meverif:
        print(f"Upperbound for the sum of dimension gaps for resolutions in psi_12 : {mse}")
    if nb_lines > 0:
        print("delta for particular bits:")
        for jmin, jmax, ec in delta_lines:
            print(f"   >>Bits from {jmin} to {jmax} delta_i={ec}")
    if tupverif:
        print("Verification of DELTA( " +
              ", ".join(str(s_list[j]) for j in range(1, d)) +
              f"{s_list[d]})")
    print("=" * 68)
    for j, gen_list in enumerate(gen_lists):
        print(f"- Component {j + 1}: {gen_list[0].name()} ")
    sys.stdout.flush()

    me_test = EquidistributionTest(L=Lmax, delta=delta, mse=mse, meverif=meverif, method=method)
    tests = [me_test, tuplets_test]
    return comb, tests, nbtries


# ═══════════════════════════════════════════════════════════════════════════
# Display helpers
# ═══════════════════════════════════════════════════════════════════════════

def _format_header(nb_comp, Seed1, Seed2, temperings) -> str:
    lines = []
    lines.append("=" * 68)
    lines.append("SUMMARY OF THE SEARCH PARAMETERS\n")
    lines.append(f"Computer : {socket.gethostname()}\n")
    lines.append(f"Seed of RNG for tempering parameters = ( {Seed1:12.0f}, {Seed2:12.0f} )\n")
    lines.append("1 component:" if nb_comp == 1 else f"{nb_comp} components:")
    for trans_list in temperings:
        if trans_list:
            lines.append("  Tempering transformations:")
            for t in trans_list:
                lines.append(f"   * {t.name}")
    return "\n".join(lines)


def _format_search_summary(tests, nbtries, has_tempering, gen_lists) -> str:
    lines = []
    if has_tempering:
        lines.append(f"Number of tries per combined generator : {nbtries}")
    for test in tests:
        if isinstance(test, EquidistributionTest):
            lines.append(f"Upperbound for the sum of dimension gaps for resolutions in psi_12 : {test.mse}")
    lines.append("=" * 68)
    for j, gen_list in enumerate(gen_lists):
        if gen_list is not None:
            lines.append(f"- Component {j + 1}: {gen_list[0].name()} ")
    return "\n".join(lines)


def _format_current_comb(comb) -> str:
    lines = []
    lines.append(_EQ66)
    lines.append(f"Number of points   : 2^({comb.k_g})")
    lines.append("")
    for j, comp in enumerate(comb.components):
        gen = comb[j]
        poly_bv = gen.char_poly()
        hw = bin(poly_bv._val).count('1') + 1
        lines.append(f"hammingweigth = {hw}")
        lines.append(f"{gen.name()}:")
        disp = gen._cpp_gen.display_str()
        for line in disp.split('\n'):
            if line.lstrip().startswith('w=') or line.lstrip().startswith(' w='):
                lines.append(f"{line}  hamingweight poly = {hw}")
            else:
                lines.append(line)
        comp_str = comp.display()
        if comp_str:
            lines.append(comp_str)
        if j < comb.J - 1:
            lines.append(_EQ66dash)
        else:
            lines.append(_EQ66)
    return "\n".join(lines)


def _format_summary(nbgen, nbME, nbCF, nbsel, elapsed) -> str:
    lines = []
    lines.append("\n===========================")
    lines.append(f"   Total   =  {nbgen:10d}  ")
    lines.append("")
    lines.append(f"     ME    =  {nbME:10d}  ")
    lines.append(f"   CF-ME   =  {nbCF:10d}  ")
    lines.append(f"  retained =  {nbsel:10d}  ")
    lines.append("---------------------------")
    lines.append(f" CPU (sec) =   {elapsed:5.2f}       ")
    lines.append("===========================")
    return "\n".join(lines)
