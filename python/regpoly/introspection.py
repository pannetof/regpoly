# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Francois Panneton, Ph.D.

"""Read-only metadata helpers about generator families and tempering types.

Thin Python shim over the C++ factory specs. Lets `regpoly_web` and
other consumers ask *"what parameters does family X take?"* or *"can
family Y be enumerated?"* without importing `_regpoly_cpp` directly
(the web layer cannot import the C++ extension — that contract is
enforced by `import-linter`).

Every function returns plain Python dicts so callers don't have to
know about the underlying `ParamSpec` C++ struct.
"""

from __future__ import annotations

from typing import Any

from regpoly_cpp import _regpoly_cpp as _cpp


def get_gen_param_specs(family: str) -> list[dict[str, Any]]:
    """Return the parameter specs for a generator family.

    Parameters
    ----------
    family
        Canonical ``-Gen`` class name (or any alias accepted
        by the C++ factory).

    Returns
    -------
    list of dict
        List of parameter spec dicts. Each spec has the keys:

        - ``name`` (`str`): parameter name.
        - ``type`` (`str`): C++ type (`"int"`, `"uint32_t"`, etc.).
        - ``structural`` (`bool`): True if this parameter defines the
          period `k`.
        - ``has_default`` (`bool`), ``default`` (`Any`): default value
          when not provided.
        - ``rand_type`` (`str`), ``rand_args`` (`list`): randomisation
          hint for non-structural params.
        - ``optimizable`` (`bool`): True if the tempering optimiser
          may mutate this parameter (bitmask `b`, `c`, etc.).

    Raises
    ------
    RuntimeError
        If `family` is not registered in the C++ factory.

    See Also
    --------
    :cpp:func:`regpoly::core::get_gen_param_specs` : the C++ free function this wraps.
    """
    return _cpp.get_gen_param_specs(family)


def get_trans_param_specs(type_name: str) -> list[dict[str, Any]]:
    """Return the parameter specs for a transformation (tempering) type.

    Parameters
    ----------
    type_name
        Transformation type identifier (e.g. ``"tempMK2"``,
        ``"permut"``).

    Returns
    -------
    list of dict
        Same dict shape as
        :func:`regpoly.introspection.get_gen_param_specs`.

    Raises
    ------
    RuntimeError
        If `type_name` is not a registered transformation.
    """
    return _cpp.get_trans_param_specs(type_name)


def family_is_enumerable(family: str) -> bool:
    """Return True iff the family has a finite parameter enumerator.

    Used by the search CLI's enumerate-mode and by the tempering search
    worker to choose between random sampling and exhaustive enumeration.

    Parameters
    ----------
    family
        Canonical ``-Gen`` class name (or an alias).

    Returns
    -------
    bool
        True if `family` has a registered enumerator, False otherwise.
    """
    return bool(_cpp.family_is_enumerable(family))


__all__ = [
    "get_gen_param_specs",
    "get_trans_param_specs",
    "family_is_enumerable",
]
