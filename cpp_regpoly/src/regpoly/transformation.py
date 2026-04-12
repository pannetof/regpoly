"""
transformation.py — Single wrapper class for any C++ tempering transformation.

Python handles parameter originals, randomization, and w_original logic.
C++ handles only the apply computation and display.

Randomization specs are defined in the YAML file per parameter:
    b: { random: bitvect, bits: 32 }
    p: { random: coprime, mod: 32 }
    q: { random: int, min: 0, max: 31 }

Specs can reference other param names for their arguments:
    b: { random: bitvect, bits: w }
"""

from __future__ import annotations

import math
import random as _random

import regpoly._regpoly_cpp as _cpp


class Transformation:
    """
    Wrapper around a C++ transformation object.

    Python stores the original parameter values and random specs.
    On update_params(L), Python computes new values from the specs,
    then calls C++ update() to set them in place.
    """

    def __init__(self, trans_type: str, params: dict, w_original: int = 0) -> None:
        self._type = trans_type
        self._orig_params = dict(params)
        self._w_original = w_original
        self._random_specs: dict = {}  # key → {random: type, ...}
        self._cpp_trans = _cpp.create_transformation(trans_type, params)

    @property
    def name(self) -> str:
        return self._cpp_trans.name()

    @property
    def w(self) -> int:
        return self._cpp_trans.w()

    @property
    def w_original(self) -> int:
        return self._w_original

    def display(self) -> str:
        return self._cpp_trans.display_str()

    def update_params(self, L: int) -> None:
        """Compute w from L, generate random values from specs, update C++."""
        if self._w_original == -1:
            w = L
        else:
            w = min(self._w_original, L)

        params = dict(self._orig_params)
        params["w"] = w

        for key, spec in self._random_specs.items():
            params[key] = self._generate(spec, params)

        self._cpp_trans.update(params)

    def copy(self) -> Transformation:
        new = Transformation.__new__(Transformation)
        new._type = self._type
        new._orig_params = dict(self._orig_params)
        new._w_original = self._w_original
        new._random_specs = dict(self._random_specs)
        new._cpp_trans = self._cpp_trans.copy()
        return new

    @staticmethod
    def _generate(spec: dict, params: dict) -> int:
        """Generate a random value from a spec, resolving param references."""
        kind = spec["random"]

        if kind == "bitvect":
            bits = _resolve(spec["bits"], params)
            return _random.getrandbits(bits) & ((1 << bits) - 1)

        elif kind == "int":
            lo = _resolve(spec.get("min", 0), params)
            hi = _resolve(spec.get("max", 1), params)
            return _random.randint(lo, hi)

        elif kind == "coprime":
            mod = _resolve(spec["mod"], params)
            val = 0
            while math.gcd(val, mod) != 1:
                val = _random.randrange(1, mod)
            return val

        raise ValueError(f"Unknown random type: {kind}")

    # -- YAML reader ------------------------------------------------------

    @classmethod
    def from_yaml(cls, filename: str) -> tuple:
        """Read transformations from a YAML file. Returns (transformations, mk_opt)."""
        import yaml

        with open(filename) as f:
            data = yaml.safe_load(f)

        mk_opt = False
        transformations = []
        for entry in data["transformations"]:
            trans_type = entry["type"]

            params = {}
            random_specs = {}
            for k, v in entry.items():
                if k == "type":
                    continue
                if isinstance(v, dict) and "random" in v:
                    random_specs[k] = v
                    params[k] = 0  # placeholder
                else:
                    params[k] = v

            w_original = params.get("w", 0)
            t = cls(trans_type, params, w_original)
            t._random_specs = random_specs
            transformations.append(t)

            if entry.get("optimize", False):
                mk_opt = True
        return transformations, mk_opt


def _resolve(val, params: dict):
    """If val is a string, look it up in params. Otherwise return as-is."""
    if isinstance(val, str):
        return params[val]
    return val
