"""Rolling-rate helper for SSE progress payloads.

The v1 `rate` field is preserved as **instantaneous** (tries / elapsed)
for backwards compat with CLI consumers. The v2 `rate_rolling_5s` field
is computed from a per-run rolling window so the live UI sees a stable
tries-per-second value rather than a noisy instantaneous rate.

Worker globals: `_rates: dict[run_id, RollingRate]` populated lazily by
`for_run(run_id)`. Drop the entry on `event: end` via `drop_run(run_id)`.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from typing import Deque


@dataclass
class _Sample:
    t: float
    tries: int


@dataclass
class RollingRate:
    window_sec: float = 5.0
    samples: Deque[_Sample] = field(default_factory=deque)

    def observe(self, t: float, tries: int) -> None:
        self.samples.append(_Sample(t=t, tries=tries))
        cutoff = t - self.window_sec
        while self.samples and self.samples[0].t < cutoff:
            self.samples.popleft()

    def rate(self, now: float | None = None) -> float:
        if not self.samples:
            return 0.0
        if now is None:
            now = self.samples[-1].t
        cutoff = now - self.window_sec
        # Drop stale samples on read too in case observe() hasn't been
        # called recently.
        while self.samples and self.samples[0].t < cutoff:
            self.samples.popleft()
        if len(self.samples) < 2:
            return 0.0
        first, last = self.samples[0], self.samples[-1]
        dt = last.t - first.t
        if dt <= 0:
            return 0.0
        return float(last.tries - first.tries) / dt


_rates: dict[int, RollingRate] = {}


def for_run(run_id: int, window_sec: float = 5.0) -> RollingRate:
    """Get-or-create the per-run rolling-rate accumulator."""
    rr = _rates.get(run_id)
    if rr is None:
        rr = RollingRate(window_sec=window_sec)
        _rates[run_id] = rr
    return rr


def drop_run(run_id: int) -> None:
    """Discard the rolling-rate state for a finished run."""
    _rates.pop(run_id, None)
