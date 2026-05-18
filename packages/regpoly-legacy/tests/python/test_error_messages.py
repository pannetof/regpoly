# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""Pin the error-message substrings the deleted C++ test_legacy_reader.cpp
relied on. Future drift in error wording will fail loudly here.
"""

from __future__ import annotations

import textwrap

import pytest

from regpoly_legacy.reader import parse_generator_specs, parse_transformation_specs


def test_carry_obsolete_well_type_6_rejected(tmp_path):
    """test_legacy_reader.cpp::RejectsObsoleteWellType6 substring-matched
    "obsolete". Hand-author a minimal fixture exercising that path."""
    # carry header: w r p; then 1 generator with 8 type-blocks
    # Use raw_type=6 (the multi-shift extension that has no paper Mi).
    # Layout per row: 3 m's + 8 × (1 type + 3 pi + 3 pu) = 59 tokens.
    row_blocks = []
    for k in range(8):
        # First block uses raw_type=6 (obsolete); rest use raw_type=7 (→ Mi=0)
        rt = "6" if k == 0 else "7"
        row_blocks.append(f"{rt} 0 0 0 00000000 00000000 00000000")
    row = " ".join(["1", "2", "3"] + row_blocks)
    p = tmp_path / "obsolete.dat"
    p.write_text(textwrap.dedent(f"""\
        carry
        32 16 0
        1
        {row}
    """))
    with pytest.raises(RuntimeError, match="obsolete"):
        parse_generator_specs(str(p), 32)


def test_unknown_tag_is_explicit(tmp_path):
    p = tmp_path / "nope.dat"
    p.write_text("not-a-real-tag\n0\n")
    with pytest.raises(RuntimeError, match="unknown generator tag"):
        parse_generator_specs(str(p), 32)


def test_unknown_transformation_type_is_explicit(tmp_path):
    p = tmp_path / "trans.dat"
    p.write_text("1\nweird_trans_type 32 0 0\n")
    with pytest.raises(RuntimeError, match="unknown transformation type"):
        parse_transformation_specs(str(p))


def test_marsaxorshift_unknown_type_is_explicit(tmp_path):
    p = tmp_path / "mx.dat"
    p.write_text("marsaxorshift\n1 1 32\n1\n999\n")
    with pytest.raises(RuntimeError, match="unknown marsaxorshift type"):
        parse_generator_specs(str(p), 32)
