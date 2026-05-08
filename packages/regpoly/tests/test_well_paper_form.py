# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Byte-for-byte equivalence test for WELL constructed via the
paper-aligned `matrices` map.

Two construction paths must produce identical output streams:

1. The legacy reader path: read shared/legacy_parameters/carry32_624_31_final.dat,
   which writes a StructMap into Params via the .dat-format token decoder.
2. The structured-kwarg path: pass `matrices={T0: {M: 3, t: -25}, ...}`
   directly to Generator.create from Python.

Both go through WELLGen::from_params with the same MatrixEntry array.
The byte-for-byte fixture at packages/regpoly-cpp/tests/fixtures/well19937a_1024.bin
was captured pre-rename with the deterministic seed (bit 0 set, all
others zero); after the structured-Params refactor the same fixture
must still match.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from regpoly.core.generator import Generator
from regpoly_cpp._regpoly_cpp import BitVect


FIXTURE = (
    Path(__file__).resolve().parents[2]
    / "regpoly-cpp"
    / "tests"
    / "fixtures"
    / "well19937a_1024.bin"
)

# WELL19937a from paper Table II.
WELL19937A_KWARGS = {
    "w": 32,
    "r": 624,
    "p": 31,
    "m1": 70,
    "m2": 179,
    "m3": 449,
    "matrices": {
        "T0": {"M": 3, "t": -25},
        "T1": {"M": 3, "t":  27},
        "T2": {"M": 2, "t":   9},
        "T3": {"M": 3, "t":   1},
        "T4": {"M": 1},
        "T5": {"M": 3, "t":  -9},
        "T6": {"M": 3, "t": -21},
        "T7": {"M": 3, "t":  21},
    },
}


def _generate_1024_bytes(gen) -> bytes:
    """Pull 1024 32-bit outputs (little-endian) from a deterministically
    seeded generator. Matches the capture procedure in
    packages/regpoly-cpp/tests/test_well_byte_for_byte.cpp."""
    cpp_gen = gen._cpp_gen
    k = cpp_gen.k() if callable(cpp_gen.k) else cpp_gen.k
    seed = BitVect(k)
    seed.set_bit(0, 1)
    cpp_gen.init(seed)

    out = bytearray()
    for _ in range(1024):
        cpp_gen.next()
        bv = cpp_gen.get_output()
        word = 0
        for b in range(min(32, bv.nbits())):
            if bv.get_bit(b):
                word |= 1 << b
        out += int(word).to_bytes(4, "little")
    return bytes(out)


def test_fixture_exists() -> None:
    assert FIXTURE.exists(), (
        f"WELL byte-for-byte fixture missing: {FIXTURE}. "
        "Run the C++ test_well_byte_for_byte once to capture it."
    )


def test_paper_form_matches_byte_fixture() -> None:
    expected = FIXTURE.read_bytes()
    assert len(expected) == 4096

    gen = Generator.create("WELLRNG", 32, **WELL19937A_KWARGS)
    actual = _generate_1024_bytes(gen)

    assert actual == expected, (
        "WELL19937a constructed from the paper-aligned `matrices` map "
        "diverged from the byte-for-byte fixture. The new construction "
        "path produces a different output stream than the legacy reader."
    )


def test_string_M_class_alias_works() -> None:
    """`M: "M3"` must be equivalent to `M: 3`."""
    kwargs = dict(WELL19937A_KWARGS)
    kwargs["matrices"] = {
        slot: {**entry, "M": f"M{entry['M']}"}
        for slot, entry in kwargs["matrices"].items()
    }
    gen = Generator.create("WELLRNG", 32, **kwargs)
    actual = _generate_1024_bytes(gen)
    assert actual == FIXTURE.read_bytes()


def test_legacy_flat_triple_is_rejected() -> None:
    """Old mat_types/mat_pi/mat_pu must error: either at fill_params
    (because `matrices` is now the required spec) or at the binding's
    explicit legacy-key rejection. In both cases the error message
    must point at `matrices` so the user knows where to migrate."""
    bad_kwargs = {
        "w": 32, "r": 624, "p": 31, "m1": 70, "m2": 179, "m3": 449,
        "mat_types": [3, 3, 2, 3, 1, 3, 3, 3],
        "mat_pi": [-25, 0, 0, 27, 0, 0, 9, 0, 0, 1, 0, 0,
                    0, 0, 0, -9, 0, 0, -21, 0, 0, 21, 0, 0],
        "mat_pu": [0] * 24,
    }
    with pytest.raises(Exception, match="matrices"):
        Generator.create("WELLRNG", 32, **bad_kwargs)


def test_legacy_keys_alongside_matrices_are_rejected_at_binding() -> None:
    """Even if `matrices` is provided correctly, stray legacy keys must
    be rejected by the binding with a clear migration pointer."""
    bad_kwargs = {**WELL19937A_KWARGS, "mat_types": [0] * 8}
    with pytest.raises(Exception) as excinfo:
        Generator.create("WELLRNG", 32, **bad_kwargs)
    msg = str(excinfo.value)
    assert "mat_types" in msg
    assert "matrices" in msg


def test_missing_slot_raises() -> None:
    bad = dict(WELL19937A_KWARGS)
    bad["matrices"] = {k: v for k, v in bad["matrices"].items() if k != "T3"}
    with pytest.raises(Exception, match="T3"):
        Generator.create("WELLRNG", 32, **bad)


def test_wrong_arg_raises() -> None:
    """M3 takes `t`, not `q`."""
    bad = dict(WELL19937A_KWARGS)
    bad["matrices"] = {**bad["matrices"], "T0": {"M": 3, "q": 5}}
    with pytest.raises(Exception):
        Generator.create("WELLRNG", 32, **bad)


# ─── Cost-bounded sampler (regpoly.well) ──────────────────────────────


def test_random_matrices_respects_cost_cap() -> None:
    from regpoly import well

    for cap in (1, 4, 12, 32, 64):
        for seed in range(50):
            m = well.random_matrices(32, cap, seed=seed)
            assert well.total_cost(m) <= cap, (
                f"cap={cap} seed={seed} produced cost {well.total_cost(m)}: {m}"
            )


def test_random_matrices_builds_a_valid_generator() -> None:
    """100 cost-capped samples each produce a generator that
    `Generator.create` accepts. This exercises the full bind →
    decode → MatrixEntry validator path."""
    from regpoly import well

    for seed in range(100):
        matrices = well.random_matrices(32, 12, seed=seed)
        kwargs = {
            "w": 32, "r": 16, "p": 0,
            "m1": 13, "m2": 9, "m3": 5,
            "matrices": matrices,
        }
        gen = Generator.create("WELLRNG", 32, **kwargs)
        assert gen is not None
        assert well.total_cost(matrices) <= 12


def test_random_matrices_invalid_cap_raises() -> None:
    from regpoly import well

    with pytest.raises(Exception):
        well.random_matrices(32, 0, seed=1)
    with pytest.raises(Exception):
        well.random_matrices(32, -3, seed=1)


def test_random_matrices_unsupported_w_raises() -> None:
    from regpoly import well

    with pytest.raises(Exception):
        well.random_matrices(64, 12, seed=1)
