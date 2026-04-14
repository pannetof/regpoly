"""
generateur.py — Single wrapper class for any C++ PRNG generator.

All computation is delegated to the C++ backend via _cpp_gen.
No Python subclasses needed — one class wraps any generator type.
"""

from __future__ import annotations

import json
import os
import random as _random

import regpoly._regpoly_cpp as _cpp
from regpoly.bitvect import BitVect


# ═══════════════════════════════════════════════════════════════════════════
# Precomputed primitive factors of Φ_k(2) for primitivity testing
# ═══════════════════════════════════════════════════════════════════════════

_FACTORS_DATA: dict | None = None
_FACTORS_FILE = os.path.join(os.path.dirname(__file__), "data",
                              "primitive_factors.json")


def _load_factors() -> dict:
    global _FACTORS_DATA
    if _FACTORS_DATA is None:
        if os.path.exists(_FACTORS_FILE):
            with open(_FACTORS_FILE) as f:
                _FACTORS_DATA = json.load(f)
        else:
            _FACTORS_DATA = {}
    return _FACTORS_DATA


def _get_factors_for_k(k: int) -> list[int] | None:
    """
    Return all prime factors of 2^k - 1, or None if the factorization
    is incomplete.

    Uses precomputed factors of Φ_d(2) for every divisor d of k.
    If any divisor has an incomplete factorization, returns None
    (fallback to C++ trial division).
    """
    data = _load_factors()
    if not data:
        return None

    # For Mersenne prime exponents, 2^k - 1 is itself prime.
    # No factorization from the table needed.
    if _is_mersenne_prime_exponent(k):
        return [2**k - 1]

    # Collect factors of Phi_d(2) for every divisor d of k.
    # Phi_1(2) = 1 (no prime factors), so d=1 is always safe to skip.
    factors = set()
    for d in _divisors(k):
        if d == 1:
            continue
        entry = data.get(str(d))
        if entry is None:
            return None
        if not entry["complete"]:
            return None
        for p in entry["factors"]:
            factors.add(p)
    return sorted(factors)


_MERSENNE_PRIMES = {
    2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607,
    1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213,
    19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091,
    756839, 859433, 1257787,
}


def _is_mersenne_prime_exponent(k: int) -> bool:
    return k in _MERSENNE_PRIMES


def _divisors(n: int) -> list[int]:
    """Return all divisors of n in ascending order."""
    divs = []
    for i in range(1, int(n**0.5) + 1):
        if n % i == 0:
            divs.append(i)
            if i != n // i:
                divs.append(n // i)
    return sorted(divs)


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
    "carry":         "WELLRNG",
    "Carry2Gen":     "WELLRNG",
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


# ═══════════════════════════════════════════════════════════════════════════
# Random-generation helpers (interpret rand_type / rand_args from C++ specs)
# ═══════════════════════════════════════════════════════════════════════════

def _eval_expr(expr: str, params: dict) -> int:
    """Evaluate a simple expression: literal int, param name, or 'param±N'."""
    expr = expr.strip()
    # Pure integer literal
    if expr.lstrip("-").isdigit():
        return int(expr)
    # "param-N" or "param+N"
    for op in ("-", "+"):
        if op in expr and not expr.startswith(op):
            name, offset = expr.split(op, 1)
            val = params[name.strip()]
            n = int(offset.strip())
            return (val - n) if op == "-" else (val + n)
    # Bare param name
    return params[expr]


def _generate_random(spec: dict, params: dict) -> object:
    """Generate a random value for a non-structural parameter."""
    rt = spec["rand_type"]
    ra = spec["rand_args"]

    if rt == "bitmask":
        bits = _eval_expr(ra, params)
        return _random.getrandbits(bits)

    if rt == "range":
        lo_expr, hi_expr = ra.split(",", 1)
        lo = _eval_expr(lo_expr, params)
        hi = _eval_expr(hi_expr, params)
        return _random.randint(lo, hi)

    if rt == "poly_exponents":
        k = _eval_expr(ra, params)
        num_terms = _random.randint(1, min(k - 1, 10))
        return sorted(_random.sample(range(1, k), num_terms)) + [0]

    if rt == "bitmask_vec":
        bits_param, length_param = ra.split(",", 1)
        bits = _eval_expr(bits_param, params)
        length = len(params[length_param.strip()])
        return [_random.getrandbits(bits) for _ in range(length)]

    raise ValueError(
        f"Parameter '{spec['name']}' has rand_type='{rt}' which is not "
        f"randomizable — it must be provided explicitly"
    )


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
        self._family: str | None = None
        self._create_params: dict | None = None

    @classmethod
    def parameters(cls, family: str) -> list[dict]:
        """
        Return the parameter specifications for a generator family.

        Each entry is a dict with keys: name, type, structural,
        has_default, default, rand_type, rand_args.
        """
        return list(_cpp.get_param_specs(resolve_family(family)))

    @classmethod
    def create(cls, family: str, L: int, **params) -> "Generateur":
        """
        Create a generator, randomizing missing non-structural parameters.

        Structural parameters (those that define the period k) must always
        be provided.  Non-structural parameters are filled in randomly
        using the hints declared in C++ when omitted.

        Examples::

            # All params explicit — no randomization
            gen = Generateur.create("TGFSRGen", L=32, w=32, r=3, m=1,
                                    a=0x9908b0df)

            # m and a omitted — randomized
            gen = Generateur.create("TGFSRGen", L=32, w=32, r=3)

        Accepts both C++ class names and legacy aliases.
        """
        resolved = resolve_family(family, params)
        specs = _cpp.get_param_specs(resolved)
        full = dict(params)

        for spec in specs:
            name = spec["name"]
            if name in full:
                continue
            if spec["has_default"]:
                continue  # C++ factory will use its default
            if spec["structural"]:
                raise ValueError(
                    f"Structural parameter '{name}' is required for "
                    f"{resolved} (it defines the period)"
                )
            rt = spec["rand_type"]
            if not rt or rt == "none":
                raise ValueError(
                    f"Parameter '{name}' for {resolved} must be provided "
                    f"(no random generation available)"
                )
            full[name] = _generate_random(spec, full)

        cpp_gen = _cpp.create_generator(resolved, full, L)
        gen = cls(cpp_gen)
        gen._family = resolved
        gen._create_params = full
        return gen

    @property
    def family(self) -> str | None:
        return self._family

    @property
    def create_params(self) -> dict | None:
        return self._create_params

    def structural_params(self) -> dict:
        """Return only the structural parameters (those that define k)."""
        if not self._family or not self._create_params:
            return {}
        specs = _cpp.get_param_specs(self._family)
        return {s["name"]: self._create_params[s["name"]]
                for s in specs
                if s["structural"] and s["name"] in self._create_params}

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
        g = Generateur(self._cpp_gen.copy())
        g._family = self._family
        g._create_params = dict(self._create_params) if self._create_params else None
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
        cp = self._cpp_gen.char_poly()

        # Fast path: for Mersenne prime exponents, irreducible = primitive
        if _is_mersenne_prime_exponent(self.k):
            return _cpp.is_irreducible(cp, self.k)

        factors = _get_factors_for_k(self.k)
        if factors is None:
            raise ValueError(
                f"Cannot test primitivity for k={self.k}: the complete "
                f"factorization of 2^{self.k} - 1 is not available. "
                f"Regenerate primitive_factors.json with a larger range "
                f"or provide the missing factors."
            )
        return _cpp.is_primitive_with_factors(
            cp, self.k, [str(p) for p in factors]
        )

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
