"""
parametric.py — Shared parameter handling for Generateur and Transformation.

Provides the Parametric base class with:
  - parameters(type_name): return ParamSpec list from C++
  - fill_params(specs, user_params): auto-randomize missing params
  - get_param / set_param: read/write parameter values
  - structural_params / optimizable_params: filter by spec flags
"""

from __future__ import annotations

import math as _math
import random as _random


class Parametric:
    """
    Base class for objects whose construction is driven by ParamSpec.

    Subclasses must:
      - Override ``_get_cpp_specs(type_name)`` to call the appropriate
        C++ function.
      - Override ``set_param(name, value)`` to push the new value to C++.
      - Set ``self._type_name`` (str) and ``self._params`` (dict) during
        construction.
    """

    _params: dict
    _type_name: str

    @classmethod
    def _get_cpp_specs(cls, type_name: str) -> list[dict]:
        """Return ParamSpec list from C++.  Override in subclass."""
        raise NotImplementedError

    @classmethod
    def parameters(cls, type_name: str) -> list[dict]:
        """
        Return the parameter specifications for a type.

        Each entry is a dict with keys: name, type, structural,
        has_default, default, rand_type, rand_args, optimizable.
        """
        return list(cls._get_cpp_specs(type_name))

    @staticmethod
    def fill_params(specs: list[dict], user_params: dict) -> dict:
        """
        Fill missing non-structural params with random values.

        Structural params that are missing raise ValueError.
        Non-structural params with a rand_type are auto-randomized.
        Params with has_default=True are left for C++ to fill.

        Returns a complete parameter dict.
        """
        full = dict(user_params)
        for spec in specs:
            name = spec["name"]
            if name in full:
                continue
            if spec["has_default"]:
                continue
            if spec["structural"]:
                raise ValueError(
                    f"Structural parameter '{name}' is required "
                    f"(it cannot be randomized)"
                )
            rt = spec["rand_type"]
            if not rt or rt == "none":
                raise ValueError(
                    f"Parameter '{name}' must be provided "
                    f"(no random generation available)"
                )
            full[name] = generate_random(spec, full)
        return full

    # -- Instance properties ------------------------------------------------

    @property
    def type_name(self) -> str:
        """The type name used to create this object."""
        return self._type_name

    @property
    def params(self) -> dict:
        """A copy of the current parameter dict."""
        return dict(self._params)

    # -- Parameter access ---------------------------------------------------

    def get_param(self, name: str) -> int:
        """Read the current value of a parameter."""
        return self._params[name]

    def set_param(self, name: str, value: int) -> None:
        """Set a parameter and push to C++.  Override in subclass."""
        raise NotImplementedError

    def randomize_params(self) -> None:
        """Re-randomize all non-structural params that have a rand_type.

        Structural params and params without a rand_type are left untouched.
        Each randomized param is pushed to C++ via set_param().
        """
        specs = self._get_cpp_specs(self._type_name)
        for spec in specs:
            if spec["structural"] or spec["has_default"]:
                continue
            rt = spec["rand_type"]
            if not rt or rt == "none":
                continue
            self.set_param(spec["name"], generate_random(spec, self._params))

    # -- Spec-based queries -------------------------------------------------

    def structural_params(self) -> dict:
        """Return only the structural parameters (those that define k)."""
        if not self._type_name:
            return {}
        specs = self._get_cpp_specs(self._type_name)
        return {s["name"]: self._params[s["name"]]
                for s in specs
                if s["structural"] and s["name"] in self._params}

    def optimizable_params(self) -> list[tuple[str, int]]:
        """Return parameters marked as optimizable: [(name, width), ...]."""
        if not self._type_name:
            return []
        specs = self._get_cpp_specs(self._type_name)
        w = self._params.get("w", 0)
        return [(s["name"], w) for s in specs
                if s["optimizable"] and s["name"] in self._params]


# ═══════════════════════════════════════════════════════════════════════════
# Random-generation helpers (interpret rand_type / rand_args from C++ specs)
# ═══════════════════════════════════════════════════════════════════════════

def eval_expr(expr: str, params: dict) -> int:
    """Evaluate a simple expression: literal int, param name, or 'param±N'."""
    expr = expr.strip()
    if expr.lstrip("-").isdigit():
        return int(expr)
    for op in ("-", "+"):
        if op in expr and not expr.startswith(op):
            name, offset = expr.split(op, 1)
            val = params[name.strip()]
            n = int(offset.strip())
            return (val - n) if op == "-" else (val + n)
    return params[expr]


def generate_random(spec: dict, params: dict) -> object:
    """Generate a random value for a non-structural parameter."""
    rt = spec["rand_type"]
    ra = spec["rand_args"]

    if rt == "bitmask":
        bits = eval_expr(ra, params)
        return _random.getrandbits(bits)

    if rt == "range":
        lo_expr, hi_expr = ra.split(",", 1)
        lo = eval_expr(lo_expr, params)
        hi = eval_expr(hi_expr, params)
        return _random.randint(lo, hi)

    if rt == "poly_exponents":
        k = eval_expr(ra, params)
        num_terms = _random.randint(1, min(k - 1, 10))
        return sorted(_random.sample(range(1, k), num_terms)) + [0]

    if rt == "bitmask_vec":
        bits_param, length_param = ra.split(",", 1)
        bits = eval_expr(bits_param, params)
        length = len(params[length_param.strip()])
        return [_random.getrandbits(bits) for _ in range(length)]

    if rt == "tausworthe_poly":
        # rand_args: "k,nb_terms".  When nb_terms is missing / 0 we
        # default to 3 (trinomial); when s is 0 / missing we pick a
        # random admissible s so both the poly and the decimation step
        # vary between tries.  We stash the chosen s into `params` so
        # the downstream C++ constructor uses exactly that value
        # instead of re-deriving it from the polynomial.
        k_expr, nbt_expr = ra.split(",", 1)
        k = int(eval_expr(k_expr.strip(), params))
        nbt_key = nbt_expr.strip()
        if nbt_key in params:
            nb_terms = int(params[nbt_key])
        elif nbt_key.lstrip("-").isdigit():
            nb_terms = int(nbt_key)
        else:
            nb_terms = 0
        if nb_terms <= 0:
            nb_terms = 3
        quicktaus = bool(params.get("quicktaus", True))
        L = int(params.get("L", 0))      # injected by Generateur.create
        s = int(params.get("s", 0))
        if s <= 0:
            s_max = (k - (nb_terms - 2)) if quicktaus else (k - 1)
            if s_max < 1:
                raise ValueError(
                    f"Tausworthe: nb_terms={nb_terms} too large for k={k}"
                )
            # Decimation step s must be coprime with 2^k - 1 for the
            # decimated LFSR sequence to have full period.  Re-roll
            # until we hit a coprime value; for most k the success
            # probability is > 50%, so a small retry budget suffices.
            period = (1 << k) - 1
            for _ in range(256):
                s = _random.randint(1, s_max)
                if _math.gcd(s, period) == 1:
                    break
            else:
                raise ValueError(
                    f"Tausworthe: could not find s in [1, {s_max}] with "
                    f"gcd(s, 2^{k}-1)=1 after 256 tries"
                )
            params["s"] = s    # propagate to the generator constructor
        import regpoly._regpoly_cpp as _cpp
        return list(_cpp.tausworthe_random_poly(
            k, nb_terms, quicktaus, L, s))

    raise ValueError(
        f"Parameter '{spec['name']}' has rand_type='{rt}' which is not "
        f"randomizable — it must be provided explicitly"
    )
