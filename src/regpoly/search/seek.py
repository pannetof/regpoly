# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Full-stack equidistribution search driver.

The :class:`regpoly.search.seek.Seek` class glues together YAML
parsing, generator/tempering pool construction, the test driver, and
the result writer for the canonical ``equidist.*.yaml`` flow. The
public entry points are:

- :meth:`regpoly.search.seek.Seek.from_yaml` —
  build a `Seek` from a YAML config.
- :meth:`regpoly.search.seek.Seek.run` — execute the search loop.

The pre-v2 `legacy_file:` YAML key (which embedded `.dat` parameter
pools directly) is no longer supported; author configs using the
inline ``family:`` / ``file:`` source forms instead.

Cross-layer note: the inner search loop lives in
:cpp:class:`regpoly::core::Seek` (header: ``src/include/search/seek.h``);
this module owns configuration ingestion and result output.
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
from regpoly.core.generator import Generator
from regpoly.core.transformation import Transformation
from regpoly.io.combo_builder import build_cpp_enumerator
from regpoly.io.tested_generator import save_tested_generator

_SEP      = "\n\n" + "+" * 104
_EQ66     = "=" * 66
_EQ66dash = "-" * 66


class Seek:
    """Search for combined generators with good equidistribution properties.

    Construct via ``Seek.from_yaml("search.config.yaml")``. The pre-v2
    positional-arg constructor is no longer supported.

    See Also
    --------
    :cpp:class:`regpoly::core::Seek` : the C++ search-loop counterpart.
    """

    def __init__(self) -> None:
        # Search state: per-slot pools + the C++ enumerator (built by
        # `from_yaml` or injected via `set_enumerator`). The Python
        # `ComboEnumerator` middleman is gone.
        self._gen_pools: list[list[Generator]] = []
        self._tempering_pools: list[list[Transformation]] = []
        self._cpp_enum = None      # `_cpp.ComboEnumerator`
        self._Lmax: int = 32
        self._tests: list = []
        self._nbtries: int = 1
        self._output_dir: str | None = None

    def set_enumerator(
        self,
        cpp_enum,
        gen_pools: list[list[Generator]],
        tempering_pools: list[list[Transformation]],
        Lmax: int,
        tests: list,
        nbtries: int,
    ) -> None:
        """Inject a pre-built enumerator + the Python-side pool state.

        Public replacement for the previous "write to private
        `_comb`/`_tests`/`_nbtries`" pattern. The Python
        `gen_pools` / `tempering_pools` are retained for display and
        result serialisation (which need `.type_name` / `.params` from
        the Python `Generator` / `Transformation` wrappers).
        """
        self._cpp_enum = cpp_enum
        self._gen_pools = gen_pools
        self._tempering_pools = tempering_pools
        self._Lmax = Lmax
        self._tests = tests
        self._nbtries = nbtries

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
                raise ValueError(
                    "legacy_file: (the pre-v2 .dat parameter source) is no "
                    "longer supported. Re-author this config using inline "
                    "'family:' / 'file:' source forms."
                )
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

        # Build the C++ enumerator directly from the pools — no Python
        # ComboEnumerator middleman.
        has_tempering = any(len(tl) > 0 for tl in temperings)
        if not has_tempering:
            s._nbtries = 1

        print(_format_header(nb_comp, Seed1, Seed2, temperings))

        s._gen_pools = gen_lists
        s._tempering_pools = temperings
        s._Lmax = Lmax
        s._cpp_enum = build_cpp_enumerator(gen_lists, temperings, Lmax)

        print(_format_search_summary(s._tests, s._nbtries, has_tempering, gen_lists))
        sys.stdout.flush()

        return s

    # -- Run ---------------------------------------------------------------

    def run(self) -> None:
        """Execute the search, printing results to stdout.

        The per-combo iteration loop lives in C++. Python keeps
        ownership of YAML parsing, result-object synthesis for display,
        and persistence. The C++ `_cpp.ComboEnumerator` (renamed to
        `ComboEnumerator` in a follow-up) is the single source of truth
        for the active combo's state.
        """
        cpp_enum = self._cpp_enum
        nbtries = self._nbtries
        tests = self._tests
        gen_pools = self._gen_pools
        tempering_pools = self._tempering_pools

        if cpp_enum is None:
            raise RuntimeError(
                "Seek.run: no enumerator built. Call from_yaml(...) or "
                "set_enumerator(...) first."
            )

        # Mirror each test object so the on_iter callback can rebuild
        # EquidistributionResults / TupletsResults instances that the
        # existing display path consumes.
        eq_test = next(
            (t for t in tests if isinstance(t, EquidistributionTest)), None)
        tup_test = next(
            (t for t in tests if isinstance(t, TupletsTest)), None)
        cf_test = next(
            (t for t in tests if isinstance(t, CollisionFreeTest)), None)

        test_specs = _build_test_specs(tests)

        state = {"nbsel": 0, "nbME": 0, "nbCF": 0}

        def on_prep(_cpp_c, _is_retry):
            # `_cpp_c` is the same `_cpp.ComboEnumerator` we passed in; the
            # C++ driver advances it via reset()/next(). Nothing for the
            # Python side to do here — display/save read from cpp_enum
            # directly on each on_iter.
            return

        def on_iter(_cpp_c, iter_result):
            me_results = _synth_me_results(iter_result, eq_test, cpp_enum)
            tup_results = _synth_tup_results(iter_result, tup_test)
            cf_results = _synth_cf_results(iter_result, cf_test)

            print(_format_current_comb(cpp_enum, gen_pools))

            if me_results is not None:
                if me_results.is_me():
                    msg = me_results.display()
                    if msg:
                        print(msg)
                    state["nbME"] += 1
                else:
                    print("\n  Dimension gaps for every resolution", end="")
                    table_str, _ = me_results.display_table(cpp_enum, "l")
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
                    self._output_dir, "equidist",
                    _current_combo_snapshot(cpp_enum, gen_pools, tempering_pools),
                    test_results,
                )
                print(f"  Saved: {path}")

            print(_SEP)
            sys.stdout.flush()

        cpp_result = _cpp.run_seek_search(
            cpp_enum, test_specs, nbtries,
            1,  # progress_interval (unused since on_progress is None)
            on_prep=on_prep, on_iter=on_iter, on_progress=None)

        print(_format_summary(
            cpp_result.nbgen, state["nbME"], state["nbCF"],
            state["nbsel"], cpp_result.elapsed_seconds))


# ═══════════════════════════════════════════════════════════════════════════
# Private helpers
# ═══════════════════════════════════════════════════════════════════════════

_INT_MAX_C = 2**31 - 1


def _active_python_gens(
    cpp_enum, gen_pools: list[list[Generator]],
) -> list[Generator]:
    """Return the Python `Generator` wrappers currently active in each
    slot of `cpp_enum`. Looks up `gen_pools[j][current_gen]` to find
    the wrapper carrying `.type_name` / `.params` for serialisation."""
    out = []
    for j in range(cpp_enum.J):
        idx = cpp_enum.pool(j).current_index()
        out.append(gen_pools[j][idx])
    return out


def _current_combo_snapshot(
    cpp_enum,
    gen_pools: list[list[Generator]],
    tempering_pools: list[list[Transformation]],
):
    """Build a passive snapshot of the active combo's state for the save
    path. The legacy `save_tested_generator(comb, ...)` expected a
    ComboEnumerator-shaped object; we satisfy the same surface (J,
    components, k_g, L, __getitem__) without instantiating a Python
    ComboEnumerator."""
    return _ComboSnapshot(cpp_enum, gen_pools, tempering_pools)


class _ComboSnapshot:
    """ComboEnumerator-API-shaped snapshot over a `_cpp.ComboEnumerator`."""

    def __init__(
        self,
        cpp_enum,
        gen_pools: list[list[Generator]],
        tempering_pools: list[list[Transformation]],
    ) -> None:
        self._cpp = cpp_enum
        self._active = _active_python_gens(cpp_enum, gen_pools)
        # Each "component" exposes a `.trans` attribute mirroring the
        # Python Component shape that save_tested_generator reads from.
        self.components = [
            _ComponentSnapshot(tempering_pools[j])
            for j in range(cpp_enum.J)
        ]

    @property
    def J(self) -> int:
        return self._cpp.J

    @property
    def L(self) -> int:
        return self._cpp.L

    @property
    def k_g(self) -> int:
        return self._cpp.k_g

    def __getitem__(self, j: int) -> Generator:
        return self._active[j]


class _ComponentSnapshot:
    """Component-API-shaped snapshot over a `tempering_pool[j]`."""

    def __init__(self, trans_list: list[Transformation]) -> None:
        self.trans = trans_list


# Map from the Python METHOD_* integer constants to the canonical
# string names known to C++ EquidistributionMethodRegistry. The C++ side is the source
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
    the C++ EquidistributionMethodRegistry resolves the string at run() time. No more
    parallel enum mapping between Python and C++.

    A ``TupletsTest(tupletsverif=False)`` is a no-op marker (the legacy
    seek-factory unconditionally adds one so the test list has a stable
    shape) and gets filtered out here. Letting it reach the C++ runner
    causes a segfault: with ``d=0`` and ``tup_h=[0]``, the runner
    indexes ``tup_h[1]`` (out of bounds) when sizing its output buffer.
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
            if not t.tupletsverif:
                continue  # no-op marker — skip
            spec.kind = _cpp.SeekTestKind.Tuplets
            spec.tup_d = t.d
            spec.tup_h = list(t.s) if t.s else [0]
            spec.tup_threshold = float(t.mDD)
            spec.tup_testtype = int(t.testtype)
        else:
            raise TypeError(f"Unknown test type: {type(t).__name__}")
        out.append(spec)
    return out


def _synth_me_results(iter_result, eq_test, cpp_enum):
    """Build an EquidistributionResults object the existing display
    code expects, from the SeekIterResult fields. Reads `k_g` directly
    from the active `_cpp.ComboEnumerator`."""
    if not iter_result.me_ran or eq_test is None:
        return None
    psi12 = list(_cpp.compute_psi12(cpp_enum.k_g, eq_test.L))
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


def _format_current_comb(cpp_enum, gen_pools: list[list[Generator]]) -> str:
    """Format the active combo for display.

    Reads structure from `cpp_enum` directly (`J`, `k_g`, `component(j).display()`),
    and looks up Python `Generator` wrappers via the current pool indices
    in `gen_pools` for char_poly / name / display_str access.
    """
    lines = []
    lines.append(_EQ66)
    lines.append(f"Number of points   : 2^({cpp_enum.k_g})")
    lines.append("")
    J = cpp_enum.J
    for j in range(J):
        cpp_comp = cpp_enum.pool(j)
        idx = cpp_comp.current_index()
        gen = gen_pools[j][idx]   # Python Generator wrapper
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
        comp_str = cpp_comp.display()
        if comp_str:
            lines.append(comp_str)
        if j < J - 1:
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
