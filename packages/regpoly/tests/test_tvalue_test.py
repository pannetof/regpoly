# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""
Phase 5 Python tests: TValueTest end-to-end + YAML dispatch.
"""

from __future__ import annotations

import sys
import textwrap

import pytest

import regpoly._regpoly_cpp as _cpp
from regpoly.analyses.abstract_test import AbstractTest
from regpoly.analyses.test_registry import all_specs, get
from regpoly.analyses.tvalue_test import (
    METHOD_NIEDERREITER_PIRSIC,
    METHOD_SCHMID,
    TValueTest,
)
from regpoly.analyses.tvalue_results import TValueResults
from regpoly.core.combination import Combination
from regpoly.core.generator import Generator


# ── Registry ─────────────────────────────────────────────────────────


def test_registry_knows_tvalue():
    spec = get("tvalue")
    assert spec.test_cls is TValueTest
    assert spec.results_cls is TValueResults


def test_registry_lists_canonical_tests():
    names = {spec.name for spec in all_specs()}
    assert {"equidistribution", "collision_free", "tuplets", "tvalue"} <= names


def test_registry_unknown_test_raises():
    with pytest.raises(KeyError):
        get("not_a_real_test")


# ── TValueTest construction ──────────────────────────────────────────


def test_construct_with_defaults():
    t = TValueTest(s_max=4, max_t_sum=10)
    assert t.s_max == 4
    assert t.max_t_sum == 10
    assert t.method == METHOD_SCHMID
    assert len(t.delta) == 5  # 0..s_max inclusive


def test_construct_rejects_smax_below_two():
    with pytest.raises(ValueError):
        TValueTest(s_max=1, max_t_sum=10)


def test_construct_with_string_method():
    t = TValueTest(s_max=4, max_t_sum=10, method="niederreiter_pirsic")
    assert t.method == METHOD_NIEDERREITER_PIRSIC


# ── _from_params (YAML) ──────────────────────────────────────────────


def test_from_params_minimal():
    t = TValueTest._from_params({"s_max": 4}, Lmax=4)
    assert t.s_max == 4
    assert t.max_t_sum == sys.maxsize
    assert t.method == METHOD_SCHMID


def test_from_params_with_delta_rules():
    t = TValueTest._from_params(
        {
            "s_max": 6,
            "max_t_sum": 5,
            "delta": [
                {"from": 2, "to": 4, "max": 0},
                {"from": 5, "to": 6, "max": 2},
            ],
        },
        Lmax=10,
    )
    assert t.delta[2] == 0
    assert t.delta[3] == 0
    assert t.delta[4] == 0
    assert t.delta[5] == 2
    assert t.delta[6] == 2


# ── from_yaml registry-driven dispatch ───────────────────────────────


def test_from_yaml_dispatches_tvalue(tmp_path):
    yaml_text = textwrap.dedent(
        """
        tests:
          - type: tvalue
            s_max: 4
            max_t_sum: 99
        """
    )
    path = tmp_path / "tv.yaml"
    path.write_text(yaml_text)
    tests = AbstractTest.from_yaml(str(path), Lmax=4)
    assert len(tests) == 1
    assert isinstance(tests[0], TValueTest)
    assert tests[0].s_max == 4


def test_from_yaml_unknown_type_raises(tmp_path):
    path = tmp_path / "bad.yaml"
    path.write_text("tests:\n  - type: not_a_real_test\n")
    with pytest.raises(ValueError, match="Unknown test type"):
        AbstractTest.from_yaml(str(path), Lmax=4)


# ── Generator.create routes Sobol / Niederreiter ─────────────────────


def test_generator_create_sobol():
    g = Generator.create("SobolNet", L=2, s_max=2)
    assert g._type_name == "SobolNet"
    assert g.L == 2
    assert g.k == 2


def test_generator_create_niederreiter():
    g = Generator.create("NiederreiterF2Gen", L=2, s_max=2)
    assert g._type_name == "NiederreiterF2Gen"
    assert g.L == 2


def test_generator_create_sobol_alias():
    g = Generator.create("sobol", L=2, s_max=2)
    assert g._type_name == "SobolNet"


# ── End-to-end: TValueTest.run on a SobolNet ─────────────────────────


def _make_combination(family: str, m: int, s_max: int) -> Combination:
    gen = Generator.create(family, L=m, s_max=s_max)
    return Combination.single(gen)


def test_run_on_sobol_m2_s2_gives_zero():
    C = _make_combination("SobolNet", m=2, s_max=2)
    res = TValueTest(s_max=2, max_t_sum=sys.maxsize).run(C)
    assert res.verified
    assert res.tvals[2] == 0
    assert res.se == 0


def test_run_on_niederreiter_m2_s2_gives_zero():
    C = _make_combination("NiederreiterF2Gen", m=2, s_max=2)
    res = TValueTest(s_max=2, max_t_sum=sys.maxsize).run(C)
    assert res.verified
    assert res.tvals[2] == 0


def test_run_short_circuits_on_delta():
    C = _make_combination("SobolNet", m=2, s_max=2)
    res = TValueTest(s_max=2, max_t_sum=sys.maxsize, delta=[0, 0, -1]).run(C)
    assert not res.verified


def test_dual_method_raises_at_runtime():
    C = _make_combination("SobolNet", m=2, s_max=2)
    t = TValueTest(s_max=2, max_t_sum=sys.maxsize, method="niederreiter_pirsic")
    with pytest.raises(RuntimeError, match="not yet implemented"):
        t.run(C)


# ── Architectural claim: existing EquidistributionTest runs on a DigitalNet ─


def test_equidistribution_runs_unchanged_on_sobol():
    from regpoly.analyses.equidistribution_test import EquidistributionTest

    C = _make_combination("SobolNet", m=4, s_max=4)
    eqt = EquidistributionTest(L=4, delta=[sys.maxsize] * 5, mse=sys.maxsize)
    res = eqt.run(C)
    # The architectural claim: this returns a meaningful result
    # without any digital-net-aware code in the equidist machinery.
    assert res.verified
    assert res.se >= 0


def test_equidistribution_runs_unchanged_on_niederreiter():
    from regpoly.analyses.equidistribution_test import EquidistributionTest

    C = _make_combination("NiederreiterF2Gen", m=4, s_max=4)
    eqt = EquidistributionTest(L=4, delta=[sys.maxsize] * 5, mse=sys.maxsize)
    res = eqt.run(C)
    assert res.verified
    assert res.se >= 0


# ── Results serialisation ────────────────────────────────────────────


def test_results_roundtrip():
    r = TValueResults(s_max=3, tvals=[0, 0, 0, 1], se=1,
                      verified=True, delta=[sys.maxsize] * 4,
                      max_t_sum=sys.maxsize, method="schmid")
    d = r.to_dict()
    r2 = TValueResults.from_dict(d)
    assert r2.s_max == 3
    assert r2.tvals == [0, 0, 0, 1]
    assert r2.se == 1
    assert r2.verified


def test_results_display_contains_table():
    r = TValueResults(s_max=3, tvals=[0, 0, 1, 2], se=3,
                      verified=True, delta=[sys.maxsize] * 4,
                      max_t_sum=sys.maxsize, method="schmid")
    out = r.display()
    assert "t-value profile" in out
    assert "se = 3" in out
