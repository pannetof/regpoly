/* Small JS helpers used across the regpoly web UI.
 * Kept minimal — HTMX + Alpine.js handle most behavior declaratively. */

document.addEventListener('htmx:sseMessage', (evt) => {
    // Expose SSE payloads to Alpine.js components listening on this channel.
    try {
        const data = JSON.parse(evt.detail.data);
        const el = evt.detail.elt;
        if (el && el.__x) { el.__x.$data.progress = data; }
    } catch (_err) { /* non-JSON message — ignore */ }
});

function fmtHex(n, width) {
    if (n === null || n === undefined) return '';
    const v = BigInt(n);
    const nibbles = Math.max(1, Math.ceil((width || 32) / 4));
    return '0x' + v.toString(16).padStart(nibbles, '0');
}
