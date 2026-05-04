"""
legacy_reader.py — thin Python shim over the C++ legacy .dat parsers.

Phase 3.3 moved the per-format dispatch tables, token tape readers,
hex-word packers, and MarsaXorshift variant expansion to C++. The
shim here only:

  * routes filename → C++ via legacy_read_generator_specs /
    legacy_read_transformation_specs,
  * wraps each returned (family, params_dict) tuple with
    Generator.create() / Transformation.create() so the Python
    wrapper's _params dict stays populated for downstream
    randomisation paths (e.g. the optimizer's safe_masks).

A C++-only consumer can call regpoly_legacy::read_generators /
regpoly_legacy::read_transformations directly and bypass this shim
entirely (see packages/regpoly-cpp/src/include/io/legacy_reader.h).
"""

from __future__ import annotations

from regpoly_cpp import _regpoly_cpp as _cpp


class LegacyReader:
    """Parsers for legacy text-based generator and transformation files."""

    @staticmethod
    def read_generators(filename: str, L: int) -> list:
        from regpoly.core.generator import Generator
        specs = _cpp.legacy_read_generator_specs(filename, L)
        return [Generator.create(family, L, **params) for family, params in specs]

    @staticmethod
    def read_transformations(filename: str) -> tuple:
        """Read transformations from the legacy text format.

        Returns (transformations, mk_opt).
        """
        from regpoly.core.transformation import Transformation
        specs, mk_opt = _cpp.legacy_read_transformation_specs(filename)
        return ([Transformation.create(t_type, **params) for t_type, params in specs],
                mk_opt)
