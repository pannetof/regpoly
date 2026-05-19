# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Single Python wrapper class for any C++ PRNG generator.

All computation is delegated to the C++ backend via the underlying
``_cpp_gen`` handle; no Python subclasses are needed — one `Generator`
class wraps any family.

The companion C++ class is `regpoly_cpp::Generator` (header:
`packages/regpoly-cpp/src/include/core/generator.h`).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import regpoly._regpoly_cpp as _cpp
from regpoly.core.bitvect import BitVect
from regpoly.core.parametric import Parametric

if TYPE_CHECKING:
    from regpoly.core.matrix import BitMatrix

# Family-name resolution: short tags and legacy class names → canonical
# C++ class names used by the factory. The factory itself (factory.cpp)
# accepts both old and new canonical names, so the post-rename names
# below are also valid keys for round-trip lookups.
_FAMILY_ALIASES: dict[str, str] = {
    # short tags / legacy lowercase
    "polylcg":       "PolyLCGGen",
    "taus":          "TauswortheGen",
    "taus2":         "TauswortheGen",
    "tgfsr":         "TGFSRGen",
    "MT":            "MTGen",
    "marsaxorshift": "MarsaXorshiftGen",
    "carry":         "WELLGen",
    "Carry2Gen":     "WELLGen",
    "tinymt":        "TinyMT32Gen",
    "TinyMT":        "TinyMT32Gen",
    "rmt":           "RMT64Gen",
    "rmt64":         "RMT64Gen",
    "dsfmt":         "DSFMTGen",
    "dSFMT":         "DSFMTGen",
    # pre-rename C++ class names → new canonical
    "PolyLCG":         "PolyLCGGen",
    "Tausworthe":      "TauswortheGen",
    "TGFSR":           "TGFSRGen",
    "MersenneTwister": "MTGen",
    "GenF2wBase":      "F2wBaseGen",
    "GenF2wLFSR":      "F2wLFSRGen",
    "GenF2wPolyLCG":   "F2wPolyLCGGen",
    "WELLRNG":         "WELLGen",
    "MELG":            "MELGGen",
    "SFMT":            "SFMTGen",
    "dSFMTGen":        "DSFMTGen",
    "MTGP":            "MTGPGen",
    "TinyMT32":        "TinyMT32Gen",
    "RMT64":           "RMT64Gen",
}


def resolve_family(family: str, params: dict | None = None) -> str:
    """Translate a legacy/short family name to the canonical C++ class name.

    The special case ``"genf2w"`` inspects ``params["type"]`` to choose
    between ``"F2wLFSRGen"`` and ``"F2wPolyLCGGen"``.

    Parameters
    ----------
    family
        A family name. Accepts canonical ``-Gen`` class names,
        pre-rename class names (e.g. ``"MersenneTwister"``), and
        legacy lowercase tags (e.g. ``"tgfsr"``).
    params
        Optional parameter dict. Only consulted for the
        ``"genf2w"`` disambiguation case.

    Returns
    -------
    str
        The canonical ``-Gen`` class name as registered in the C++ factory.
        Unknown names are returned unchanged (the factory will raise on
        them later, with a clearer error message).
    """
    if family == "genf2w":
        if params and params.get("type") == "lfsr":
            return "F2wLFSRGen"
        return "F2wPolyLCGGen"
    return _FAMILY_ALIASES.get(family, family)


class Generator(Parametric):
    """Python wrapper around any C++ generator.

    A `Generator` instance binds one C++ generator object and exposes
    the runtime API every search loop and analysis needs: state
    initialisation, single-step iteration, characteristic polynomial
    extraction, transition-matrix readout, primitivity test.

    Construct via the :meth:`regpoly.core.generator.Generator.create`
    class method (which routes through the C++ factory). The default
    constructor (taking a raw `cpp_gen` handle) is not part of the
    public API.

    `Generator` is itself a Python iterator: ``next(gen)`` advances one
    step and returns the next output word as a
    :class:`regpoly.core.bitvect.BitVect` of width ``L``.

    Attributes
    ----------
    k
        Total state size in bits (degree of the characteristic
        polynomial).
    L
        Output resolution in bits — every step emits an ``L``-bit
        word.

    See Also
    --------
    :cpp:class:`regpoly::core::Generator` : the C++ implementation this wraps
        (header: ``src/include/core/generator.h``).
    """

    @classmethod
    def _get_cpp_specs(cls, type_name: str) -> list[dict]:
        return _cpp.get_gen_param_specs(resolve_family(type_name))

    def __init__(self, cpp_gen) -> None:
        """Wrap an already-constructed C++ generator handle.

        Most callers should use :meth:`regpoly.core.generator.Generator.create`
        rather than this constructor directly — the factory is what
        validates parameters and fills randomized defaults.
        """
        self._cpp_gen = cpp_gen
        self.k: int = cpp_gen.k()
        self.L: int = cpp_gen.L()
        self._type_name: str = ""
        self._params: dict = {}

    @classmethod
    def create(cls, family: str, L: int, **params: Any) -> "Generator":
        """Build a generator, randomising any missing non-structural parameter.

        Structural parameters (those that define the period ``k``) must
        always be provided. Non-structural parameters are filled in
        randomly using the per-family hints declared in C++ when
        omitted from the call.

        Parameters
        ----------
        family
            Generator family — canonical ``-Gen`` class name
            or any accepted alias (see
            :func:`regpoly.core.generator.resolve_family`).
        L
            Output resolution in bits (1 ≤ L ≤ word width).
        **params
            Structural and non-structural parameters. Any
            non-structural parameter not provided will be
            randomised; structural omissions raise.

        Returns
        -------
        Generator
            A constructed :class:`regpoly.core.generator.Generator`
            instance bound to a fresh C++ object.

        Raises
        ------
        ValueError
            If a structural parameter is missing or any
            provided parameter is invalid for the family.
        RuntimeError
            If the C++ factory rejects the family name.

        Examples
        --------
        All params explicit — no randomisation:

        >>> gen = Generator.create(                       # doctest: +SKIP
        ...     "TGFSRGen", L=32, w=32, r=3, m=1, a=0x9908b0df,
        ... )

        `m` and `a` omitted — both randomised on each call:

        >>> gen = Generator.create("TGFSRGen", L=32, w=32, r=3)  # doctest: +SKIP
        """
        resolved = resolve_family(family, params)
        specs = _cpp.get_gen_param_specs(resolved)
        # Expose L to randomizers (e.g. tausworthe_poly) via the params
        # dict; strip it afterwards so we don't pass it to the C++ factory.
        work = dict(params)
        work.setdefault("L", L)
        full = cls.fill_params(specs, work)
        full.pop("L", None)

        cpp_gen = _cpp.create_generator(resolved, full, L)
        gen = cls(cpp_gen)
        gen._type_name = resolved
        gen._params = full
        return gen

    @classmethod
    def create_at_index(cls, family: str, L: int,
                        enumerator: Any, idx: int, **fixed: Any) -> "Generator":
        """Build a generator for one index of an exhaustive enumeration.

        The enumerator's `at(idx)` contribution wins over stale values
        in ``fixed`` — enumerated axes are authoritative. Skips the
        `fill_params` step because the enumerator emits every
        non-structural value the factory needs.

        Parameters
        ----------
        family
            Generator family (canonical or alias).
        L
            Output resolution in bits.
        enumerator
            A C++ enumerator with `at(idx) -> dict` semantics.
        idx
            Index into the enumeration (0-based).
        **fixed
            Parameter overrides. Enumerated parameters take
            precedence over equivalent keys here.

        Returns
        -------
        Generator
            A constructed :class:`regpoly.core.generator.Generator`
            at the given enumeration index.
        """
        resolved_family = resolve_family(family, fixed)
        enumerated = enumerator.at(idx)
        full = {**fixed, **enumerated}
        full.pop("L", None)
        cpp_gen = _cpp.create_generator(resolved_family, full, L)
        gen = cls(cpp_gen)
        gen._type_name = resolved_family
        gen._params = full
        return gen

    # -- Generator-specific methods ---------------------------------------

    def name(self) -> str:
        """Return the family's display name as reported by C++.

        Equivalent to `regpoly_cpp::Generator::name()`.

        Returns
        -------
        str
            The C++-side display string (typically the canonical
            ``-Gen`` class name).
        """
        return self._cpp_gen.name()

    def display(self) -> str:
        """Return a human-readable parameter dump.

        Used by `regpoly-cli show` and by the web UI's
        generator-detail page.

        Returns
        -------
        str
            A multi-line string describing the family and every
            parameter currently bound.
        """
        return self._cpp_gen.display_str()

    def set_param(self, name: str, value: int) -> None:
        """Always raises — generator parameters are immutable post-construction.

        Parameters
        ----------
        name
            Parameter name (ignored).
        value
            New value (ignored).

        Raises
        ------
        TypeError
            Always. Build a fresh generator with
            :meth:`regpoly.core.generator.Generator.create`
            instead.
        """
        raise TypeError("Generator parameters cannot be changed after creation")

    def initialize_state(self, init_bv: BitVect) -> BitVect:
        """Set the generator's internal state and return the first output.

        Parameters
        ----------
        init_bv
            A ``k``-bit `BitVect` providing the seed state.

        Returns
        -------
        BitVect
            The output word produced immediately after seeding (an
            ``L``-bit `BitVect`).
        """
        self._cpp_gen.init(_cpp.BitVect.from_int(init_bv.n, init_bv._val))
        return BitVect(self.L, self._cpp_gen.get_output().to_int())

    def __next__(self) -> BitVect:
        """Advance one step and return the next output word.

        Returns
        -------
        BitVect
            The ``L``-bit output word produced by this step.
        """
        self._cpp_gen.next()
        return BitVect(self.L, self._cpp_gen.get_output().to_int())

    def __iter__(self) -> "Generator":
        """Return self — `Generator` is its own iterator."""
        return self

    def copy(self) -> "Generator":
        """Return an independent deep copy of this generator.

        The new instance shares no mutable state with the original;
        advancing one does not affect the other.

        Returns
        -------
        Generator
            A new `Generator` with cloned C++ state plus a copy of the
            bound parameter dict.
        """
        g = Generator(self._cpp_gen.copy())
        g._type_name = self._type_name
        g._params = dict(self._params)
        return g

    def char_poly(self) -> BitVect:
        """Return the characteristic polynomial of the state-update matrix.

        Returns
        -------
        BitVect
            A ``k``-bit `BitVect` whose binary expansion encodes the
            polynomial coefficients (LSB = constant term).
        """
        bv = self._cpp_gen.char_poly()
        return BitVect(self.k, bv.to_int())

    def is_full_period(self) -> bool:
        """Return ``True`` iff the characteristic polynomial is primitive.

        A primitive characteristic polynomial implies the generator's
        period is the maximum 2^k − 1. The check delegates to the C++
        driver in `algebra/primitivity.cpp`, which has a Mersenne
        fast-path and a Cunningham-style factor-table lookup.

        Returns
        -------
        bool
            ``True`` if maximum period; ``False`` otherwise.
        """
        return _cpp.is_full_period(self._cpp_gen)

    def transition_matrix(self) -> "BitMatrix":
        """Return the state-update matrix as a `BitMatrix`.

        Returns
        -------
        BitMatrix
            A ``k × k`` `BitMatrix` ``A`` such that the next state is
            ``A * s`` for current state ``s`` (column vector).
        """
        from regpoly.core.matrix import BitMatrix
        K = self.k
        cpp_rows = self._cpp_gen.transition_matrix()
        A = BitMatrix(K, K)
        for i in range(K):
            A._rows[i] = cpp_rows[i].to_int()
        return A

    # -- YAML reader ------------------------------------------------------

    @classmethod
    def from_yaml(cls, filename: str, L: int) -> list["Generator"]:
        """Build a list of generators from a YAML file.

        The file must have the shape::

            family: <family_name>
            common: { ... }           # optional shared parameters
            generators: [ ... ]       # list of per-generator parameters

        Each entry in ``generators`` is merged with the file-level
        family params and the `common:` block, then passed to
        :meth:`regpoly.core.generator.Generator.create`.

        Parameters
        ----------
        filename
            Path to the YAML file.
        L
            Output resolution in bits, applied to every generator
            in the file.

        Returns
        -------
        list of Generator
            A list of constructed `Generator` instances, in the order
            they appear in the YAML.

        Raises
        ------
        FileNotFoundError
            If `filename` cannot be opened.
        KeyError
            If the YAML omits the required ``family`` or
            ``generators`` keys.
        """
        import yaml

        with open(filename) as f:
            data = yaml.safe_load(f)

        family = data["family"]
        family_params = {k: v for k, v in data.items()
                         if k not in ("family", "common", "generators")}
        common = {**family_params, **data.get("common", {})}

        generators = []
        for entry in data["generators"]:
            merged = {**common, **entry}
            generators.append(cls.create(family, L, **merged))
        return generators
