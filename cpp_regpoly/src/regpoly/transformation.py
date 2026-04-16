"""
transformation.py — Single wrapper class for any C++ tempering transformation.

C++ handles the apply computation and display.
Python handles parameter storage and randomization via the Parametric base.

Missing non-structural parameters are auto-randomized from ParamSpec,
just like Generateur.create().
"""

from __future__ import annotations

import regpoly._regpoly_cpp as _cpp
from regpoly.parametric import Parametric


class Transformation(Parametric):
    """
    Wrapper around a C++ transformation object.

    Use ``Transformation.create(type, **params)`` to construct.
    Missing non-structural parameters with a rand_type are
    auto-randomized, just like Generateur.create().
    """

    @classmethod
    def _get_cpp_specs(cls, type_name: str) -> list[dict]:
        return _cpp.get_trans_param_specs(type_name)

    def __init__(self, cpp_trans, trans_type: str, params: dict) -> None:
        self._cpp_trans = cpp_trans
        self._type_name = trans_type
        self._params = dict(params)

    @classmethod
    def create(cls, trans_type: str, **params) -> "Transformation":
        """
        Create a transformation, randomizing missing non-structural parameters.

        Structural parameters must be provided.  Non-structural parameters
        with a rand_type are auto-randomized when omitted.

        Examples::

            # All params explicit
            t = Transformation.create("tempMK", w=32, eta=7, mu=15,
                                      b=0x9D2C5680, c=0xEFC60000)

            # b and c omitted — randomized
            t = Transformation.create("tempMK", w=32, eta=7, mu=15)
        """
        specs = _cpp.get_trans_param_specs(trans_type)
        full = cls.fill_params(specs, params)
        cpp_trans = _cpp.create_transformation(trans_type, full)
        return cls(cpp_trans, trans_type, full)

    # -- Transformation-specific methods ----------------------------------

    @property
    def name(self) -> str:
        return self._cpp_trans.name()

    @property
    def w(self) -> int:
        return self._cpp_trans.w()

    def display(self) -> str:
        return self._cpp_trans.display_str()

    def set_param(self, name: str, value: int) -> None:
        """Set a parameter and push to C++."""
        self._params[name] = value
        self._cpp_trans.update(self._params)

    def copy(self) -> Transformation:
        new = Transformation.__new__(Transformation)
        new._type_name = self._type_name
        new._params = dict(self._params)
        new._cpp_trans = self._cpp_trans.copy()
        return new

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
            for k, v in entry.items():
                if k == "type":
                    continue
                if isinstance(v, dict) and "random" in v:
                    continue  # omit — will be randomized by fill_params
                else:
                    params[k] = v

            t = cls.create(trans_type, **params)
            transformations.append(t)

            if entry.get("optimize", False):
                mk_opt = True
        return transformations, mk_opt
