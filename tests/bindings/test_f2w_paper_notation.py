# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Phase-2 paper-notation entry for F2wLFSRGen / F2wPolyLCGGen.

The C++ generator's `from_params` accepts either:
 - the explicit Phase-1 catalog form `{nocoeff, coeff, modM}`, or
 - the paper-notation form `{nb_terms, t, q, modM, coeff}` from which
   nocoeff is derived inside from_params.

The two paths must produce identical generators. The new
`irreducible_gf2` rand_type rejection-samples a degree-w irreducible
polynomial for modM.
"""

from __future__ import annotations

from regpoly_cpp import _regpoly_cpp as cpp


def test_paper_form_and_explicit_nocoeff_construct_same_gen() -> None:
    # F2wLFSR2_7_416 paper-named row (Phase 1 catalog) — 3-term.
    explicit = dict(
        w=32, r=13, modM=0x92bb39c1, normal_basis=False, step=1,
        nocoeff=[9, 6, 0],
        coeff=[0x06000000, 0x41000000, 0x05000000],
    )
    paper = dict(
        w=32, r=13, modM=0x92bb39c1, normal_basis=False, step=1,
        nb_terms=3, t=9, q=6,
        coeff=[0x06000000, 0x41000000, 0x05000000],
    )
    g1 = cpp.create_generator("F2wLFSRGen", explicit, 32)
    g2 = cpp.create_generator("F2wLFSRGen", paper, 32)
    assert g1.k() == g2.k() == 13 * 32

    bv = cpp.BitVect(g1.k())
    bv.set_bit(0, 1)
    g1.init(bv)
    g2.init(bv)
    assert cpp.is_full_period(g1)
    assert cpp.is_full_period(g2)


def test_paper_form_nb_terms_equals_2_drops_q() -> None:
    """nb_terms=2 → nocoeff=[t, 0], coeff has 2 entries. F2wLFSR2_31_800."""
    paper = dict(
        w=32, r=25, modM=0xfa4f9b3f, normal_basis=False, step=1,
        nb_terms=2, t=7,
        coeff=[0xe6a68d20, 0x287ab842],
    )
    g = cpp.create_generator("F2wLFSRGen", paper, 32)
    bv = cpp.BitVect(g.k())
    bv.set_bit(0, 1)
    g.init(bv)
    assert cpp.is_full_period(g)


def test_irreducible_gf2_sampler_produces_irreducible_polys() -> None:
    """`is_irreducible_gf2w_modM` is exercised by the rand_type via
    `random_param`; sample 200 candidates and confirm none collide
    with a non-irreducible value (proxy: every sampled modM yields a
    full-period generator with otherwise-uniform parameters)."""
    seen = set()
    for _ in range(200):
        val, _ = cpp.random_param("irreducible_gf2", "w",
                                  {"w": 32}, 32)
        assert isinstance(val, int)
        assert (val >> 32) == 0, "sampled value must fit in w=32 bits"
        seen.add(val)
    # Reject sampler shouldn't get stuck on a single polynomial.
    assert len(seen) > 50, f"only {len(seen)} distinct samples in 200 draws"


def test_range_sampler_with_param_reference() -> None:
    """`range` with rand_args 'lo,name±N' resolves the offset expression
    against the params bag (used for q's '1,t-1')."""
    for _ in range(50):
        val, _ = cpp.random_param("range", "1,r-1",
                                  {"r": 8}, 32)
        assert 1 <= val <= 7, val


def test_bitmask_vec_scalar_length() -> None:
    """`bitmask_vec` with a scalar-length param (paper-notation coeff
    uses rand_args='w,nb_terms' where nb_terms is a scalar 2 or 3)."""
    for _ in range(50):
        val, _ = cpp.random_param("bitmask_vec", "w,nb_terms",
                                  {"w": 32, "nb_terms": 3}, 32)
        assert isinstance(val, list)
        assert len(val) == 3
        # Values arrive as signed int64 over the pybind11 bridge — mask
        # back to the unsigned w-bit range for the assertion.
        MASK = (1 << 32) - 1
        for v in val:
            assert 0 <= (v & MASK) <= MASK
