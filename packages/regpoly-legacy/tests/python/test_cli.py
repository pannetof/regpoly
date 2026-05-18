# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Subprocess smoke tests for the ``regpoly-legacy`` console script."""

from __future__ import annotations

import subprocess
import sys

from pathlib import Path

import pytest


def fixtures_dir() -> Path:
    return Path(__file__).resolve().parents[2] / "shared" / "legacy_parameters"


def _run(*args):
    return subprocess.run(
        [sys.executable, "-m", "regpoly_legacy.cli", *args],
        capture_output=True, text=True, timeout=60,
    )


def test_help_lists_three_subcommands():
    out = _run("--help")
    assert out.returncode == 0
    for sub in ("info", "trans", "seek"):
        assert sub in out.stdout


def test_info_mt19937():
    fdir = fixtures_dir()
    out = _run("info", str(fdir / "MT19937.dat"))
    assert out.returncode == 0, out.stderr
    assert "MersenneTwister" in out.stdout
    assert "Generators: 1" in out.stdout


def test_trans_returns_two_rows():
    fdir = fixtures_dir()
    out = _run("trans", str(fdir / "trans32.dat"))
    assert out.returncode == 0, out.stderr
    assert "Transformations: 2" in out.stdout
    assert "mk_opt=True" in out.stdout


def test_seek_argument_count_mismatch():
    fdir = fixtures_dir()
    # nb_comp=2 but only 1 gen file supplied
    out = _run("seek", "2", str(fdir / "example1c"), str(fdir / "96_1.dt"))
    assert out.returncode != 0
    assert "Expected 2 generator file(s)" in out.stderr


@pytest.mark.slow
def test_seek_runs_to_completion():
    """A short seek run with the canonical example1c + 96_1.dt fixture pair.
    Marked slow — actually runs the search, so it can take seconds."""
    fdir = fixtures_dir()
    out = _run("seek", "1", str(fdir / "example1c"), str(fdir / "96_1.dt"))
    # The search prints to stdout; we don't pin the wording, just the exit.
    assert out.returncode == 0, out.stderr
