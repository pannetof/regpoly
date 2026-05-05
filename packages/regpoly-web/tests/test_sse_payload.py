"""Phase 3 red — SSE dual-emit contract.

v1 (`data: …`) events stay unchanged; v2 adds a named `event: progress`
channel with `rate_rolling_5s`, `eta_seconds`, `cum_finds`, `t` for
primitive runs and `best_score`, `best_params`, `best_id` for tempering
runs. The `event: end` terminal is preserved.

This file uses the SSE parser in tests/_sse.py and seeds search_progress
rows directly, then opens the SSE endpoint with a short timeout to
collect a snapshot of events.
"""

from __future__ import annotations

import json

from tests._sse import parse_sse_events


def _stream_text(client, url: str, timeout: float = 0.6) -> str:
    """Pull a few seconds of SSE body off the endpoint, then close."""
    with client.stream("GET", url, timeout=timeout) as resp:
        assert resp.status_code == 200
        chunks: list[bytes] = []
        try:
            for chunk in resp.iter_bytes():
                chunks.append(chunk)
                # Read at most ~16 KB to avoid hanging.
                if sum(len(c) for c in chunks) > 16 * 1024:
                    break
        except Exception:
            pass
    return b"".join(chunks).decode("utf-8", errors="replace")


def test_v1_progress_unnamed_event_unchanged(seeded_client) -> None:
    """v1 SSE endpoint emits only `data:` blocks (no `event:` field on
    progress rows). The `end` event is the only named event."""
    body = _stream_text(seeded_client, "/api/primitive-searches/1/progress")
    events = parse_sse_events(body)
    progress_events = [e for e in events if e.name == "message"]
    end_events = [e for e in events if e.name == "end"]
    # If the run has no progress rows seeded, the SSE may emit only the
    # terminal `end` event. Either case is fine — we just assert the
    # contract.
    if progress_events:
        for e in progress_events:
            payload = json.loads(e.data)
            # v1 contract: tries_done present.
            assert "tries_done" in payload
    # Terminal event present (the seeded run has status='completed').
    assert end_events, "v1 SSE must still emit the named `end` event"
