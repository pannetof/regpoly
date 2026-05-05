"""SSE event-stream parser for contract tests.

Splits a text event-stream body into a list of (name, data) tuples
where unnamed `data: …` events appear as ("message", data) per the
EventSource default-channel convention.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Event:
    name: str
    data: str

    def is_default(self) -> bool:
        return self.name == "message"


def parse_sse_events(stream: str | bytes) -> list[Event]:
    """Parse a `text/event-stream` body into Event objects.

    Args:
        stream: The full body text (or bytes) — typically what
            httpx.Client.stream returns concatenated.

    Returns:
        List of Event(name, data) preserving order. Unnamed events
        carry name="message".
    """
    if isinstance(stream, bytes):
        text = stream.decode("utf-8", errors="replace")
    else:
        text = stream

    out: list[Event] = []
    cur_name: str | None = None
    cur_data_lines: list[str] = []
    for raw_line in text.splitlines():
        line = raw_line
        if line == "":
            if cur_data_lines:
                name = cur_name or "message"
                out.append(Event(name=name, data="\n".join(cur_data_lines)))
            cur_name = None
            cur_data_lines = []
            continue
        if line.startswith(":"):
            continue
        if ":" in line:
            field, _, value = line.partition(":")
            if value.startswith(" "):
                value = value[1:]
        else:
            field, value = line, ""
        if field == "event":
            cur_name = value
        elif field == "data":
            cur_data_lines.append(value)
    if cur_data_lines:
        out.append(Event(name=cur_name or "message",
                         data="\n".join(cur_data_lines)))
    return out
