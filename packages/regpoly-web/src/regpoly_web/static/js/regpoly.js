/* regpoly-web v2 client-side glue.
 *
 * Loaded after app.js (legacy helpers) and before Alpine. Provides:
 *   - Theme toggle + persistence (already pre-applied in <head> to
 *     avoid a flash; this module wires the toggle button + dispatches
 *     a `regpoly:theme` event so other components can react).
 *   - KaTeX renderer on DOMContentLoaded.
 *   - Hex copy-on-click delegate (event delegation on body).
 *   - Keyboard shortcuts (?, g d/g g/.../n p/n t, /, j/k, Enter, Esc).
 *   - char_poly_card client-side exponent extraction.
 */
(function () {
    "use strict";

    const THEME_KEY = "regpoly.theme";

    /* ---------------- Theme ---------------- */
    function setTheme(theme) {
        document.documentElement.setAttribute("data-bs-theme", theme);
        try { localStorage.setItem(THEME_KEY, theme); } catch (e) { /* nope */ }
        document.dispatchEvent(new CustomEvent("regpoly:theme", { detail: { theme } }));
    }
    function currentTheme() {
        return document.documentElement.getAttribute("data-bs-theme") || "light";
    }
    function bindThemeToggle() {
        document.querySelectorAll("[data-theme-toggle]").forEach(btn => {
            btn.addEventListener("click", () => {
                setTheme(currentTheme() === "dark" ? "light" : "dark");
            });
        });
    }

    /* ---------------- KaTeX ---------------- */
    function renderMath(root) {
        if (typeof window.renderMathInElement !== "function") return;
        try {
            window.renderMathInElement(root || document.body, {
                delimiters: [
                    { left: "$$", right: "$$", display: true },
                    { left: "\\(", right: "\\)", display: false },
                ],
                throwOnError: false,
            });
        } catch (e) { console.warn("KaTeX render failed", e); }
    }

    /* ---------------- char_poly_card exponents -----------------
     * Reads data-poly-hex / data-poly-k / data-poly-mode and writes
     * the rendered KaTeX content into the placeholder span. Keeps
     * the macro static while dynamic computation lives here. */
    function bitsFromHex(hex) {
        // Returns a sorted descending list of bit positions where the
        // polynomial has set bits. The hex is interpreted as a
        // big-endian bitstring; bit 0 is the LSB.
        const s = String(hex || "").replace(/^0x/i, "");
        if (!s) return [];
        const bits = [];
        const total = s.length * 4;
        for (let i = 0; i < s.length; i++) {
            const nibble = parseInt(s[i], 16) || 0;
            for (let b = 0; b < 4; b++) {
                if (nibble & (1 << (3 - b))) {
                    // Position from MSB → bit index from LSB
                    bits.push(total - 1 - (i * 4 + b));
                }
            }
        }
        return bits.sort((a, b) => b - a);
    }
    function renderCharPolyCards() {
        document.querySelectorAll(".poly-katex").forEach(el => {
            const hex = el.getAttribute("data-poly-hex");
            const k = parseInt(el.getAttribute("data-poly-k") || "0", 10);
            const mode = el.getAttribute("data-poly-mode");
            if (!hex || !k) return;
            const bits = bitsFromHex(hex);
            // Filter to terms ≤ k (defensive) and drop 0 (handled below).
            const positive = bits.filter(b => b > 0 && b <= k);
            const has_one = bits.includes(0);
            let tex;
            if (mode === "exponents") {
                const terms = positive.map(b => "x^{" + b + "}");
                if (has_one) terms.push("1");
                tex = "\\(" + (terms.length ? terms.join(" + ") : "0") + "\\)";
            } else if (mode === "xnotation") {
                const terms = positive.map(b => "x^{" + b + "}");
                tex = "\\(" + (terms.length ? terms.join(" \\cdot ") : "1") + "\\)";
            } else {
                return;
            }
            el.textContent = tex;
            renderMath(el);
        });
    }

    /* ---------------- Hex copy delegate -----------------
     * Handles clicks anywhere on a [data-copy] element. Reads the full
     * hex value from data-hex (set by the hex_value macro) so truncated
     * displays still copy the full value. */
    async function copyToClipboard(text) {
        try {
            await navigator.clipboard.writeText(text);
            return true;
        } catch (e) {
            // Fallback for non-HTTPS / older browsers.
            const ta = document.createElement("textarea");
            ta.value = text;
            ta.style.position = "fixed";
            ta.style.opacity = "0";
            document.body.appendChild(ta);
            ta.focus(); ta.select();
            try { document.execCommand("copy"); document.body.removeChild(ta); return true; }
            catch (e2) { document.body.removeChild(ta); return false; }
        }
    }
    function bindCopyDelegate() {
        document.body.addEventListener("click", async (ev) => {
            const target = ev.target.closest("[data-copy], [data-copy-text]");
            if (!target) return;
            const text = target.getAttribute("data-copy-text") || target.getAttribute("data-hex") || target.textContent.trim();
            const ok = await copyToClipboard(text);
            if (ok) {
                target.classList.add("copied");
                setTimeout(() => target.classList.remove("copied"), 800);
                // Polite live region for screen readers.
                announce("hex value copied");
            }
        });
    }
    function announce(msg) {
        let live = document.getElementById("regpoly-live");
        if (!live) {
            live = document.createElement("div");
            live.id = "regpoly-live";
            live.setAttribute("aria-live", "polite");
            live.setAttribute("role", "status");
            live.style.position = "absolute";
            live.style.left = "-10000px";
            document.body.appendChild(live);
        }
        live.textContent = msg;
    }

    /* ---------------- Keyboard shortcuts -----------------
     * Chord state machine: first key arms a chord, second key within
     * 1.2 s consumes it. Pressing Esc clears any armed chord. Ignores
     * key presses inside form fields. */
    const CHORDS = {
        "g d": "/",
        "g g": "/generators",
        "g t": "/tested-generators",
        "g s": "/searches",
        "g l": "/library",
        "g o": "/tools",
        "n p": "/primitive-search?family=MTGen",
        "n t": "/tempering-search",
    };
    let armed = null;
    let armedTimer = null;
    function inFormField(target) {
        if (!target) return false;
        const tag = target.tagName;
        return tag === "INPUT" || tag === "TEXTAREA" || tag === "SELECT" ||
               (target.isContentEditable === true);
    }
    function jkRowNav(direction) {
        const rows = Array.from(document.querySelectorAll("table tbody tr"));
        if (!rows.length) return;
        const current = rows.findIndex(r => r.classList.contains("kbd-active"));
        let next;
        if (current < 0) {
            next = direction > 0 ? 0 : rows.length - 1;
        } else {
            rows[current].classList.remove("kbd-active");
            next = current + direction;
            if (next < 0) next = 0;
            if (next >= rows.length) next = rows.length - 1;
        }
        rows[next].classList.add("kbd-active");
        rows[next].scrollIntoView({ block: "nearest" });
    }
    function openActiveRow() {
        const row = document.querySelector("table tbody tr.kbd-active");
        if (!row) return;
        const link = row.querySelector("a[href]");
        if (link) location.href = link.getAttribute("href");
    }
    function bindKeyboard() {
        document.addEventListener("keydown", (e) => {
            if (e.metaKey || e.ctrlKey || e.altKey) return;
            if (inFormField(e.target)) {
                if (e.key === "Escape") e.target.blur();
                return;
            }
            // ? cheat sheet
            if (e.key === "?" || (e.shiftKey && e.key === "/")) {
                const modal = document.getElementById("kbd-help-modal");
                if (modal) {
                    modal.classList.add("show");
                    modal.style.display = "block";
                    modal.setAttribute("aria-hidden", "false");
                }
                e.preventDefault();
                return;
            }
            if (e.key === "Escape") {
                const modal = document.querySelector(".modal.show, [data-kbd-help-modal].show");
                if (modal) {
                    modal.classList.remove("show");
                    modal.style.display = "none";
                    modal.setAttribute("aria-hidden", "true");
                }
                if (armedTimer) { clearTimeout(armedTimer); armedTimer = null; }
                armed = null;
                return;
            }
            // / focuses the chip-toolbar / search input
            if (e.key === "/") {
                const target = document.querySelector("[data-filter-focus], input[type=search]");
                if (target) {
                    target.focus();
                    e.preventDefault();
                }
                return;
            }
            // j/k row nav
            if (e.key === "j") { jkRowNav(1); e.preventDefault(); return; }
            if (e.key === "k") { jkRowNav(-1); e.preventDefault(); return; }
            if (e.key === "Enter") {
                if (document.querySelector("table tbody tr.kbd-active")) {
                    openActiveRow();
                    e.preventDefault();
                }
                return;
            }
            // Chords
            if (armed) {
                const key = armed + " " + e.key;
                clearTimeout(armedTimer);
                armed = null;
                if (CHORDS[key]) {
                    location.href = CHORDS[key];
                    e.preventDefault();
                }
                return;
            }
            if (e.key === "g" || e.key === "n") {
                armed = e.key;
                armedTimer = setTimeout(() => { armed = null; armedTimer = null; }, 1200);
                e.preventDefault();
            }
        });
    }

    /* ---------------- DOMContentLoaded ---------------- */
    document.addEventListener("DOMContentLoaded", () => {
        bindThemeToggle();
        bindCopyDelegate();
        bindKeyboard();
        renderCharPolyCards();
        // KaTeX may not be loaded yet (it's deferred); poll briefly.
        let tries = 0;
        function tryRender() {
            if (typeof window.renderMathInElement === "function") {
                renderMath(document.body);
            } else if (tries++ < 30) {
                setTimeout(tryRender, 100);
            }
        }
        tryRender();
    });
})();
