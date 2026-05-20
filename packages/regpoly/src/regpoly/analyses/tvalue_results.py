# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Francois Panneton, Ph.D.

"""
tvalue_results.py — t-value test results container.

Holds the per-dimension t-value profile `tvals[s]`, aggregate `se`,
per-s caps `delta[s]`, the `max_t_sum` threshold, and the `verified`
flag. Mirrors EquidistributionResults' shape and predicates.
"""

from __future__ import annotations

from regpoly.analyses.abstract_results import AbstractTestResults


class TValueResults(AbstractTestResults):
    """Result of a (t,m,s)-net t-value profile.

    Attributes
    ----------
    s_max
        Largest dimension queried (profile is over `s = 2..s_max`).
    tvals
        Per-dimension t-value, indexed `0..s_max`. `tvals[0]` and
        `tvals[1]` are unused (t(1) = 0 by definition).
    se
        Cumulative t-value sum `Σ_{s=2..s_max} tvals[s]`. Lower is
        better; `se = 0` means every queried dimension achieved a
        `(0, m, s)`-net (perfect equidistribution).
    delta
        Per-s cap on `tvals[s]` indexed `0..s_max`; `delta[s] = INT_MAX`
        means "no cap".
    max_t_sum
        Cap on cumulative `se`; INT_MAX for no cap.
    method
        String name of the computation method (`"schmid"` by default;
        `"niederreiter_pirsic"` is registered but not yet implemented).
    _verified
        True iff the run completed without tripping a per-s `delta`
        or `max_t_sum` short-circuit. Exposed via the `verified`
        property declared on the base class.
    """

    def __init__(
        self,
        s_max: int,
        tvals: list[int],
        se: int,
        verified: bool,
        delta: list[int],
        max_t_sum: int,
        method: str = "schmid",
    ) -> None:
        self.s_max = s_max
        self.tvals = list(tvals)
        self.se = se
        self._verified = verified
        self.delta = list(delta)
        self.max_t_sum = max_t_sum
        self.method = method

    @property
    def verified(self) -> bool:
        return self._verified

    def is_optimal(self) -> bool:
        """True iff `tvals[s] == 0` for every `s = 2..s_max`."""
        return all(self.tvals[s] == 0 for s in range(2, self.s_max + 1))

    def display(self) -> str:
        """Return a human-readable summary."""
        lines = []
        lines.append(f"t-value profile  (method={self.method}, "
                     f"verified={self._verified})")
        lines.append(f"  s   | t(s)  | delta[s] | gap")
        lines.append(f"  ----+-------+----------+-----")
        for s in range(2, self.s_max + 1):
            cap = self.delta[s] if s < len(self.delta) else None
            gap = max(0, self.tvals[s] - (cap if cap is not None
                                           and cap < 2**30 else 0))
            cap_str = "inf" if (cap is None or cap > 2**30) else str(cap)
            lines.append(
                f"  {s:<3} | {self.tvals[s]:<5} | {cap_str:<8} | {gap}"
            )
        lines.append(f"  se = {self.se} "
                     f"(max_t_sum = {'inf' if self.max_t_sum > 2**30 else self.max_t_sum})")
        return "\n".join(lines)

    def to_dict(self) -> dict:
        """Round-trippable dict for YAML serialisation."""
        return {
            "method": self.method,
            "s_max": self.s_max,
            "tvals": list(self.tvals),
            "se": self.se,
            "verified": self._verified,
            "delta": list(self.delta),
            "max_t_sum": self.max_t_sum,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "TValueResults":
        return cls(
            s_max=d["s_max"],
            tvals=d["tvals"],
            se=d["se"],
            verified=d["verified"],
            delta=d["delta"],
            max_t_sum=d["max_t_sum"],
            method=d.get("method", "schmid"),
        )
