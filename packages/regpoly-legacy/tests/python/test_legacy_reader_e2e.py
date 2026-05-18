# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""End-to-end smoke: a few canonical fixtures parse through LegacyReader
and yield Generator / Transformation objects of the expected family."""

from __future__ import annotations

import pytest

from pathlib import Path

from regpoly_legacy.reader import LegacyReader


def fixtures_dir() -> Path:
    return Path(__file__).resolve().parents[2] / "shared" / "legacy_parameters"


@pytest.mark.parametrize(
    "fname,expected_name_substring,expected_count",
    [
        ("MT19937.dat",                "Mersenne Twister",   1),
        ("32poly.dat",                 "Polynomial LCG",    13),
        ("trinomials.dat",             "Tausworthe",         9),
        ("96_1.dt",                    "TGFSR",              5),
        ("genf2w_test.dat",            "F_{2^w}",            2),
        ("marsaxorshift_test.dat",     "Marsaglia Xor-shift", 4),
        ("carry32_624_31_final.dat",   "Carry Generator",    1),
    ],
)
def test_read_generators_yields_correct_family(fname, expected_name_substring, expected_count):
    fdir = fixtures_dir()
    gens = LegacyReader.read_generators(str(fdir / fname), 32)
    assert len(gens) == expected_count, (
        f"expected {expected_count} generators in {fname}, got {len(gens)}"
    )
    actual = gens[0].name()
    assert expected_name_substring.lower() in actual.lower(), (
        f"{fname}: expected name containing '{expected_name_substring}', got '{actual}'"
    )


def test_read_transformations_returns_pair():
    fdir = fixtures_dir()
    trans, mk_opt = LegacyReader.read_transformations(str(fdir / "trans32.dat"))
    assert len(trans) == 2
    assert mk_opt is True   # trans32.dat uses tempMKopt
