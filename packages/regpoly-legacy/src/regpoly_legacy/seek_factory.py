# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""seek_from_legacy(): build a regpoly.search.seek.Seek from the legacy
positional-argument input form (test_file + per-component .dat files).

The implementation mirrors the body of the now-deleted
``Seek.from_legacy`` / ``_read_legacy`` in ``regpoly.search.seek``. It
reaches into ``Seek``'s private ``_comb`` / ``_tests`` / ``_nbtries``
attributes — acceptable since both modules live in the same workspace
and the previous ``Seek.from_legacy`` classmethod did the same.
"""

from __future__ import annotations

import random
import sys

from regpoly.analyses.equidistribution_test import (
    METHOD_DUALLATTICE,
    METHOD_MATRICIAL,
    METHOD_NOTHING,
    EquidistributionTest,
)
from regpoly.analyses.tuplets_results import _MAX_TYPE
from regpoly.analyses.tuplets_test import TupletsTest
from regpoly.core.combination import Combination
from regpoly.search.seek import (
    Seek,
    _compute_seeds,
    _format_header,
)

from regpoly_legacy.reader import LegacyReader


def _read_legacy(nb_comp, test_file, gen_data_files):
    """Read legacy C-format files and return (comb, tests, nbtries).

    Verbatim port of the now-deleted ``regpoly.search.seek._read_legacy``.
    """
    with open(test_file, encoding="latin-1") as f:
        lines = [ln.split('#')[0].strip() for ln in f]
    lines = [ln for ln in lines if ln]
    i = 0

    parts = lines[i].split(); i += 1
    seed1 = int(parts[0])
    seed2 = int(parts[1]) if seed1 != -1 else 0
    # seed2 captured for parity with the C-format; unused below since
    # _compute_seeds rederives both Seed1/Seed2 from seed1.
    _ = seed2

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
        tuplets_test = TupletsTest(
            tupletsverif=True, d=d, s=s_list, mDD=mDD, testtype=_MAX_TYPE,
        )
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

    Seed1, Seed2, seed = _compute_seeds(seed1, seed2)
    random.seed(seed)

    print(_format_header(nb_comp, Seed1, Seed2, temperings))

    comb = Combination.CreateFromFiles(gen_lists, Lmax, temperings)

    if has_tempering:
        print(f"Number of tries per combined generator : {nbtries}")
    if meverif:
        print(
            f"Upperbound for the sum of dimension gaps for resolutions "
            f"in psi_12 : {mse}"
        )
    if nb_lines > 0:
        print("delta for particular bits:")
        for jmin, jmax, ec in delta_lines:
            print(f"   >>Bits from {jmin} to {jmax} delta_i={ec}")
    if tupverif:
        print(
            "Verification of DELTA( "
            + ", ".join(str(s_list[j]) for j in range(1, d))
            + f"{s_list[d]})"
        )
    print("=" * 68)
    for j, gen_list in enumerate(gen_lists):
        print(f"- Component {j + 1}: {gen_list[0].name()} ")
    sys.stdout.flush()

    me_test = EquidistributionTest(
        L=Lmax, delta=delta, mse=mse, meverif=meverif, method=method,
    )
    tests = [me_test, tuplets_test]
    return comb, tests, nbtries


def seek_from_legacy(
    nb_comp: int,
    test_file: str,
    gen_data_files: list[str],
) -> Seek:
    """Build a ``Seek`` from legacy C-format files.

    Replaces the now-deleted ``Seek.from_legacy`` classmethod.
    """
    s = Seek()
    s._comb, tests, s._nbtries = _read_legacy(nb_comp, test_file, gen_data_files)
    s._tests = tests
    return s
