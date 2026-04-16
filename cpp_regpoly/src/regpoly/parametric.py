"""
parametric.py — Shared parameter handling for Generateur and Transformation.

Provides the Parametric base class with:
  - parameters(type_name): return ParamSpec list from C++
  - fill_params(specs, user_params): auto-randomize missing params
  - get_param / set_param: read/write parameter values
  - structural_params / optimizable_params: filter by spec flags
"""

from __future__ import annotations

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

    raise ValueError(
        f"Parameter '{spec['name']}' has rand_type='{rt}' which is not "
        f"randomizable — it must be provided explicitly"
    )
