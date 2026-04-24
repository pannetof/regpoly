"""Shared helpers for assembling a :class:`Combinaison` from a list
of component configs.

The tempering-search worker and the published-generators library both
translate ``{family, params, tempering: [...]}``-shaped descriptions
into lists of :class:`Generateur` and :class:`Transformation`
instances ready to feed into ``Combinaison.CreateFromFiles``.  This
module factors the straightforward single-generator-per-component
pattern out so both callers share one implementation.

The tempering-search worker has an additional pool-sharing layer (one
pool can carry many candidate generators and several components can
share the same pool); that logic stays in the worker.
"""

from __future__ import annotations

from regpoly.generateur import Generateur
from regpoly.transformation import Transformation


def build_tempering_chain(
    tempering_conf: list[dict],
) -> list[Transformation]:
    """Turn a list of ``{type, …params}`` dicts into a list of
    constructed :class:`Transformation` objects, in order.

    An empty input yields an empty chain.  Any missing ``type`` key
    raises :class:`KeyError` — the caller is expected to validate
    inputs before calling this.
    """
    chain: list[Transformation] = []
    for step in tempering_conf:
        ttype = step["type"]
        params = {k: v for k, v in step.items() if k != "type"}
        chain.append(Transformation.create(ttype, **params))
    return chain


def build_combinaison_inputs(
    components: list[dict], Lmax: int,
) -> tuple[list[list[Generateur]], list[list[Transformation]]]:
    """Convert library-shape component list into the pair of lists
    expected by :meth:`Combinaison.CreateFromFiles`.

    ``components``: a list of dicts, each with keys ``family``,
    ``params`` (dict), and ``tempering`` (list of step dicts); the
    ``L`` key is ignored because every component is rebuilt at
    ``Lmax`` (the combined generator must emit uniform output width).

    Returns ``(gen_pools, temperings)`` where each pool is a
    single-element list containing the component's generator.  The
    output matches the shape the tempering-search worker produces when
    pool-sharing is disabled.
    """
    gen_pools: list[list[Generateur]] = []
    temperings: list[list[Transformation]] = []
    for comp in components:
        gen = Generateur.create(comp["family"], Lmax, **comp["params"])
        gen_pools.append([gen])
        temperings.append(build_tempering_chain(comp.get("tempering") or []))
    return gen_pools, temperings


__all__ = ["build_tempering_chain", "build_combinaison_inputs"]
