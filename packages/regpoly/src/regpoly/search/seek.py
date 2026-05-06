# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""
seek.py — Seek class: runs the equidistribution search.

Supports two input formats:
  - YAML: single config file (equidist.*.yaml)
  - Legacy: nb_comp + test_file + gen_data_files (C-compatible)
"""

from __future__ import annotations

import os
import random
import socket
import sys
import time

import regpoly._regpoly_cpp as _cpp
from regpoly.analyses.collision_free_results import CollisionFreeResults
from regpoly.analyses.collision_free_test import CollisionFreeTest
from regpoly.analyses.equidistribution_results import EquidistributionResults
from regpoly.analyses.equidistribution_test import (
    METHOD_DUALLATTICE,
    METHOD_HARASE,
    METHOD_MATRICIAL,
    METHOD_NOTHING,
    METHOD_NOTPRIMITIVE,
    METHOD_SIMD_NOTPRIMITIVE,
    EquidistributionTest,
)
from regpoly.analyses.tuplets_results import _MAX_TYPE, TupletsResults
from regpoly.analyses.tuplets_test import TupletsTest
from regpoly.core.combination import Combination
from regpoly.core.generator import Generator
from regpoly.core.transformation import Transformation
from regpoly.io.legacy_reader import LegacyReader
from regpoly.io.tested_generator import save_tested_generator

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
        self._comb: Combination = None
        self._tests: list = []
        self._nbtries: int = 1
        self._output_dir: str | None = None

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
        s._output_dir = search.get("output_dir", "yaml/testedgenerators")

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
                gen_list = Generator.from_yaml(path, Lmax)
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
            from regpoly.analyses.abstract_test import AbstractTest
            path = _resolve_path(tests_cfg["file"], base_dir)
            s._tests = AbstractTest.from_yaml(path, Lmax)
        else:
            s._tests = _parse_tests(tests_cfg, Lmax)

        # Build Combination
        has_tempering = any(len(tl) > 0 for tl in temperings)
        if not has_tempering:
            s._nbtries = 1

        print(_format_header(nb_comp, Seed1, Seed2, temperings))

        s._comb = Combination(J=nb_comp, Lmax=Lmax)
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
        """Execute the search, printing results to stdout. Phase 2.4b:
        the per-combo iteration loop now lives in C++; Python keeps
        ownership of YAML parsing, result-object synthesis for display,
        and persistence."""
        py_comb = self._comb
        nbtries = self._nbtries
        tests = self._tests

        if not py_comb.reset():
            print("First combined generator not found or invalid Combination")
            return

        # Mirror each test object so the on_iter callback can rebuild
        # EquidistributionResults / TupletsResults instances that the
        # existing display path consumes.
        eq_test = next(
            (t for t in tests if isinstance(t, EquidistributionTest)), None)
        tup_test = next(
            (t for t in tests if isinstance(t, TupletsTest)), None)
        cf_test = next(
            (t for t in tests if isinstance(t, CollisionFreeTest)), None)

        # Build the C++ Combination + the typed test specs for the driver.
        cpp_comb = _build_cpp_comb_from_python(py_comb)
        test_specs = _build_test_specs(tests)

        state = {"first_done": False, "nbsel": 0, "nbME": 0, "nbCF": 0}

        def on_prep(_cpp_c, is_retry):
            # Lockstep: the C++ comb starts at the same position as
            # py_comb (built by _build_cpp_comb_from_python). On
            # subsequent fresh-combo iterations (not retries), advance
            # py_comb to match the C++ side so the display path sees
            # the right state.
            if not state["first_done"]:
                state["first_done"] = True
                return
            if not is_retry:
                try:
                    next(py_comb)
                except StopIteration:
                    pass

        def on_iter(_cpp_c, iter_result):
            me_results = _synth_me_results(iter_result, eq_test, py_comb)
            tup_results = _synth_tup_results(iter_result, tup_test)
            cf_results = _synth_cf_results(iter_result, cf_test)

            print(_format_current_comb(py_comb))

            if me_results is not None:
                if me_results.is_me():
                    msg = me_results.display()
                    if msg:
                        print(msg)
                    state["nbME"] += 1
                else:
                    print("\n  Dimension gaps for every resolution", end="")
                    table_str, _ = me_results.display_table(py_comb, "l")
                    print(table_str)
            if tup_results is not None:
                msg = tup_results.display()
                if msg:
                    print(msg)
            if cf_results is not None and cf_results.verified:
                msg = cf_results.display()
                if msg:
                    print(msg)

            state["nbsel"] += 1

            if self._output_dir:
                test_results = _build_results_dict(me_results, tup_results)
                path = save_tested_generator(
                    self._output_dir, "equidist", py_comb, test_results)
                print(f"  Saved: {path}")

            print(_SEP)
            sys.stdout.flush()

        cpp_result = _cpp.run_seek_search(
            cpp_comb, test_specs, nbtries,
            1,  # progress_interval (unused since on_progress is None)
            on_prep=on_prep, on_iter=on_iter, on_progress=None)

        print(_format_summary(
            cpp_result.nbgen, state["nbME"], state["nbCF"],
            state["nbsel"], cpp_result.elapsed_seconds))


# ═══════════════════════════════════════════════════════════════════════════
# Private helpers
# ═══════════════════════════════════════════════════════════════════════════

_INT_MAX_C = 2**31 - 1


def _build_cpp_comb_from_python(py_comb: Combination):
    """Construct a C++ Combination mirroring the Python one's pool +
    transformation structure. Detects shared pools by Python list
    object identity and replays them via share_pool_with on the C++
    side (preserving C(n,k) selection semantics)."""
    cpp_comb = _cpp.Combination(py_comb.J, py_comb.Lmax)
    pool_owner: dict[int, int] = {}
    for j, py_comp in enumerate(py_comb.components):
        cpp_comp = cpp_comb.component(j)
        pool_key = id(py_comp.gens)
        if pool_key in pool_owner:
            cpp_comp.share_pool_with(cpp_comb.component(pool_owner[pool_key]))
        else:
            pool_owner[pool_key] = j
            for py_gen in py_comp.gens:
                cpp_comp.add_gen(py_gen._cpp_gen)
        for py_trans in py_comp.trans:
            cpp_comp.add_trans(py_trans._cpp_trans)
    if not cpp_comb.reset():
        raise RuntimeError(
            "_build_cpp_comb_from_python: C++ combination has no valid combo "
            "(empty pool?)")
    return cpp_comb


# Map from the Python METHOD_* integer constants to the canonical
# string names known to C++ MethodRegistry. The C++ side is the source
# of truth — this map only translates Python's integer enum (kept for
# backward compatibility with notebooks/tests) into those strings.
_EQ_METHOD_TO_NAME = {
    METHOD_MATRICIAL:         "matricial",
    METHOD_DUALLATTICE:       "lattice",
    METHOD_HARASE:            "harase",
    METHOD_NOTPRIMITIVE:      "notprimitive",
    METHOD_SIMD_NOTPRIMITIVE: "simd_notprimitive",
    METHOD_NOTHING:           "nothing",
}


def _build_test_specs(tests: list) -> list:
    """Convert each Python test instance to a SeekTestSpec for the C++
    driver. Order is preserved — the driver runs them in order and
    short-circuits on equidistribution / tuplets failure.

    Equidistribution variants are dispatched via the canonical
    SeekTestKind.Equidistribution + SeekTestSpec.method_name string;
    the C++ MethodRegistry resolves the string at run() time. No more
    parallel enum mapping between Python and C++.
    """
    out = []
    for t in tests:
        spec = _cpp.SeekTestSpec()
        if isinstance(t, EquidistributionTest):
            spec.kind = _cpp.SeekTestKind.Equidistribution
            spec.method_name = _EQ_METHOD_TO_NAME.get(t.method, "matricial")
            spec.eq_L_max_test = t.L
            spec.eq_delta = [min(d, _INT_MAX_C) for d in t.delta]
            spec.eq_mse = min(t.mse, _INT_MAX_C)
        elif isinstance(t, CollisionFreeTest):
            spec.kind = _cpp.SeekTestKind.CollisionFree
        elif isinstance(t, TupletsTest):
            spec.kind = _cpp.SeekTestKind.Tuplets
            spec.tup_d = t.d
            spec.tup_h = list(t.s) if t.s else [0]
            spec.tup_threshold = float(t.mDD)
            spec.tup_testtype = int(t.testtype)
        else:
            raise TypeError(f"Unknown test type: {type(t).__name__}")
        out.append(spec)
    return out


def _synth_me_results(iter_result, eq_test, py_comb):
    """Build an EquidistributionResults object the existing display
    code expects, from the SeekIterResult fields."""
    if not iter_result.me_ran or eq_test is None:
        return None
    psi12 = list(_cpp.compute_psi12(py_comb.k_g, eq_test.L))
    return EquidistributionResults(
        L=eq_test.L,
        ecart=list(iter_result.me_ecart),
        psi12=psi12,
        se=iter_result.me_se,
        verified=iter_result.me_verified,
        mse=eq_test.mse,
        meverif=eq_test.meverif,
        delta=eq_test.delta,
    )


def _synth_tup_results(iter_result, tup_test):
    """Best-effort TupletsResults reconstruction. The C++ runner
    returns the firstpart_/secondpart_ aggregates; the per-tuple
    arrays (gap, DELTA, pourcentage, tuph) are not currently surfaced
    on the SeekIterResult — so for display purposes we keep them
    empty. The display() method handles empty arrays gracefully (it
    only prints the firstpart/secondpart summary)."""
    if not iter_result.tup_ran or tup_test is None:
        return None
    return TupletsResults(
        tupletsverif=True,
        tupd=tup_test.d,
        tuph=list(tup_test.s) if tup_test.s else [],
        gap=[],
        DELTA=[],
        pourcentage=[],
        firstpart_max=iter_result.tup_firstpart_max,
        firstpart_sum=iter_result.tup_firstpart_sum,
        secondpart_max=iter_result.tup_secondpart_max,
        secondpart_sum=iter_result.tup_secondpart_sum,
        treshold=tup_test.mDD,
        testtype=tup_test.testtype,
        verified=True,
    )


def _synth_cf_results(iter_result, cf_test):
    if not iter_result.cf_ran or cf_test is None:
        return None
    return CollisionFreeResults(
        ecart_cf=list(iter_result.cf_ecart_cf),
        secf=iter_result.cf_secf,
        verified=iter_result.cf_verified,
        msecf=cf_test.msecf,
    )


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
        generators.append(Generator.create(family, L, **params))
    return generators


def _parse_tempering(trans_cfg) -> list:
    if not isinstance(trans_cfg, list):
        return []
    trans_list = []
    for entry in trans_cfg:
        params = {}
        for k, v in entry.items():
            if k == "type":
                continue
            if isinstance(v, dict) and "random" in v:
                continue  # omit — will be randomized by fill_params
            else:
                params[k] = v
        t = Transformation.create(entry["type"], **params)
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

    comb = Combination.CreateFromFiles(gen_lists, Lmax, temperings)

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


def _build_results_dict(me_results, tup_results=None) -> dict:
    """Build a results dict from test results for saving."""
    results = {}
    if me_results is not None and me_results.verified:
        eq = {"se": me_results.se}
        if me_results.is_me():
            eq["status"] = "ME"
        # Store non-zero ecart values compactly
        ecart = {}
        for l in range(1, me_results.L + 1):
            if me_results.ecart[l] != 0:
                ecart[l] = me_results.ecart[l]
        if ecart:
            eq["ecart"] = ecart
        results["equidistribution"] = eq

    if tup_results is not None and tup_results.verified:
        tup = {
            "firstpart_max": tup_results.firstpart_max,
            "firstpart_sum": tup_results.firstpart_sum,
            "secondpart_max": tup_results.secondpart_max,
            "secondpart_sum": tup_results.secondpart_sum,
        }
        results["tuplets"] = tup

    return results


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
