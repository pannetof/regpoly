"""
generateur.py — Single wrapper class for any C++ PRNG generator.

All computation is delegated to the C++ backend via _cpp_gen.
No Python subclasses needed — one class wraps any generator type.
"""

from __future__ import annotations

import regpoly._regpoly_cpp as _cpp
from regpoly.bitvect import BitVect
from regpoly.parametric import Parametric


# Legacy/short names → C++ class names used by the factory.
_FAMILY_ALIASES: dict[str, str] = {
    "polylcg":       "PolyLCG",
    "taus":          "Tausworthe",
    "taus2":         "Tausworthe",
    "tgfsr":         "TGFSR",
    "TGFSRGen":      "TGFSR",   # legacy C++ class name before rename
    "MT":            "MersenneTwister",
    "matsumoto":     "MatsumotoGen",
    "marsaxorshift": "MarsaXorshiftGen",
    "AC1D":          "AC1DGen",
    "carry":         "WELLRNG",
    "Carry2Gen":     "WELLRNG",
    "xorshift":      "XorShift128",
    "XORSHIFT128":   "XorShift128",
    "tinymt":        "TinyMT32",
    "TinyMT":        "TinyMT32",
    "rmt":           "RMT64",
    "rmt64":         "RMT64",
    "dsfmt":         "dSFMTGen",
    "dSFMT":         "dSFMTGen",
}


def resolve_family(family: str, params: dict | None = None) -> str:
    """
    Translate a legacy/short family name to the C++ class name.

    The special case ``"genf2w"`` inspects ``params["type"]`` to choose
    between ``"GenF2wLFSR"`` and ``"GenF2wPolyLCG"``.
    """
    if family == "genf2w":
        if params and params.get("type") == "lfsr":
            return "GenF2wLFSR"
        return "GenF2wPolyLCG"
    return _FAMILY_ALIASES.get(family, family)


class Generateur(Parametric):
    """
    Wrapper around a C++ generator object.

    Generateur is a Python iterator: call next(gen) to advance one step.

    Attributes
    ----------
    k : int — degree (total state bits)
    L : int — output resolution in bits
    """

    @classmethod
    def _get_cpp_specs(cls, type_name: str) -> list[dict]:
        return _cpp.get_gen_param_specs(resolve_family(type_name))

    def __init__(self, cpp_gen) -> None:
        self._cpp_gen = cpp_gen
        self.k: int = cpp_gen.k()
        self.L: int = cpp_gen.L()
        self._type_name: str = ""
        self._params: dict = {}

    @classmethod
    def create(cls, family: str, L: int, **params) -> "Generateur":
        """
        Create a generator, randomizing missing non-structural parameters.

        Structural parameters (those that define the period k) must always
        be provided.  Non-structural parameters are filled in randomly
        using the hints declared in C++ when omitted.

        Examples::

            # All params explicit — no randomization
            gen = Generateur.create("TGFSR", L=32, w=32, r=3, m=1,
                                    a=0x9908b0df)

            # m and a omitted — randomized
            gen = Generateur.create("TGFSR", L=32, w=32, r=3)

        Accepts both C++ class names and legacy aliases.
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
                        enumerator, idx: int, **fixed) -> "Generateur":
        """Build a generator for the given exhaustive-enumeration index.

        The enumerator's ``at(idx)`` contribution wins over stale
        values in ``fixed`` — enumerated axes are authoritative.  Skips
        ``fill_params`` because the enumerator emits every
        non-structural value the factory needs.
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
        return self._cpp_gen.name()

    def display(self) -> str:
        return self._cpp_gen.display_str()

    def set_param(self, name: str, value: int) -> None:
        """Not supported for generators (params are fixed at creation)."""
        raise TypeError("Generator parameters cannot be changed after creation")

    def initialize_state(self, init_bv: BitVect) -> BitVect:
        self._cpp_gen.init(_cpp.BitVect.from_int(init_bv.n, init_bv._val))
        return BitVect(self.L, self._cpp_gen.get_output().to_int())

    def __next__(self) -> BitVect:
        self._cpp_gen.next()
        return BitVect(self.L, self._cpp_gen.get_output().to_int())

    def __iter__(self) -> "Generateur":
        return self

    def copy(self) -> "Generateur":
        g = Generateur(self._cpp_gen.copy())
        g._type_name = self._type_name
        g._params = dict(self._params)
        return g

    def char_poly(self) -> BitVect:
        bv = self._cpp_gen.char_poly()
        return BitVect(self.k, bv.to_int())

    def is_full_period(self) -> bool:
        """Returns True if the characteristic polynomial is primitive,
        meaning the generator has maximum period 2^k - 1.

        For Mersenne prime exponents (2^k - 1 is prime), only an
        irreducibility test is needed.  For other k, uses precomputed
        prime factors of 2^k - 1 from the Cunningham tables.
        """
        from regpoly.primitivity import is_full_period as _is_full_period
        return _is_full_period(self._cpp_gen, self.k)

    def transition_matrix(self) -> "BitMatrix":
        from regpoly.matrix import BitMatrix
        K = self.k
        cpp_rows = self._cpp_gen.transition_matrix()
        A = BitMatrix(K, K)
        for i in range(K):
            A._rows[i] = cpp_rows[i].to_int()
        return A

    # -- YAML reader ------------------------------------------------------

    @classmethod
    def from_yaml(cls, filename: str, L: int) -> list["Generateur"]:
        """
        Read generators from a YAML file.

        The file must have:
            family: <family_name>
            common: { ... }           — optional shared parameters
            generators: [ ... ]       — list of per-generator parameters

        Each generator entry is merged with family-level and common params,
        then passed to create().
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
