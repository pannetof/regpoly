"""
generateur.py — Single wrapper class for any C++ PRNG generator.

All computation is delegated to the C++ backend via _cpp_gen.
No Python subclasses needed — one class wraps any generator type.
"""

import regpoly._regpoly_cpp as _cpp
from regpoly.bitvect import BitVect


# Legacy/short names → C++ class names used by the factory.
# The C++ factory expects the exact class name (e.g. "PolyLCG").
# "genf2w" is handled specially by resolve_family(): it inspects the
# params["type"] field to select GenF2wLFSR or GenF2wPolyLCG.
_FAMILY_ALIASES: dict[str, str] = {
    "polylcg":       "PolyLCG",
    "taus":          "Tausworthe",
    "taus2":         "Tausworthe",
    "tgfsr":         "TGFSRGen",
    "MT":            "MersenneTwister",
    "matsumoto":     "MatsumotoGen",
    "marsaxorshift": "MarsaXorshiftGen",
    "AC1D":          "AC1DGen",
    "carry":         "Carry2Gen",
}


def resolve_family(family: str, params: dict | None = None) -> str:
    """
    Translate a legacy/short family name to the C++ class name.

    If *family* is already a valid C++ class name it is returned as-is.
    The special case ``"genf2w"`` inspects ``params["type"]`` to choose
    between ``"GenF2wLFSR"`` and ``"GenF2wPolyLCG"``.
    """
    if family == "genf2w":
        if params and params.get("type") == "lfsr":
            return "GenF2wLFSR"
        return "GenF2wPolyLCG"
    return _FAMILY_ALIASES.get(family, family)


class Generateur:
    """
    Wrapper around a C++ generator object.

    Generateur is a Python iterator: call next(gen) to advance one step.

    Attributes
    ----------
    k : int — degree (total state bits)
    L : int — output resolution in bits
    """

    def __init__(self, cpp_gen) -> None:
        self._cpp_gen = cpp_gen
        self.k: int = cpp_gen.k()
        self.L: int = cpp_gen.L()

    def name(self) -> str:
        return self._cpp_gen.name()

    def display(self) -> str:
        return self._cpp_gen.display_str()

    def initialize_state(self, init_bv: BitVect) -> BitVect:
        self._cpp_gen.init(_cpp.BitVect.from_int(init_bv.n, init_bv._val))
        return BitVect(self.L, self._cpp_gen.get_output().to_int())

    def __next__(self) -> BitVect:
        self._cpp_gen.next()
        return BitVect(self.L, self._cpp_gen.get_output().to_int())

    def __iter__(self) -> "Generateur":
        return self

    def copy(self) -> "Generateur":
        return Generateur(self._cpp_gen.copy())

    def char_poly(self) -> BitVect:
        bv = self._cpp_gen.char_poly()
        return BitVect(self.k, bv.to_int())

    def is_full_period(self) -> bool:
        """Returns True if the characteristic polynomial is primitive,
        meaning the generator has maximum period 2^k - 1."""
        return self._cpp_gen.is_full_period()

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
        then passed to the C++ factory.
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
            params = {**common, **entry}
            cpp_gen = _cpp.create_generator(resolve_family(family, params), params, L)
            generators.append(cls(cpp_gen))
        return generators
