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
    function paintThemeIcons() {
        const dark = currentTheme() === "dark";
        document.querySelectorAll("[data-theme-icon-light]").forEach(el => {
            el.hidden = dark;        // sun shows in dark mode (=> click for light)
        });
        document.querySelectorAll("[data-theme-icon-dark]").forEach(el => {
            el.hidden = !dark;       // moon shows in light mode (=> click for dark)
        });
    }
    function bindThemeToggle() {
        document.querySelectorAll("[data-theme-toggle]").forEach(btn => {
            btn.addEventListener("click", () => {
                setTheme(currentTheme() === "dark" ? "light" : "dark");
                paintThemeIcons();
            });
        });
        paintThemeIcons();
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

    /* ---------------- Alpine factories (P2) ---------------- */
    /* chipToolbar() — reads URL query string, exposes chip array, lets
     * x-for render chips with data-chip-* attributes the e2e tests
     * assert against. Removing a chip pushes a new URL via
     * history.pushState. Density toggle (?density=compact) is bundled.
     */
    window.chipToolbar = function chipToolbar(opts) {
        opts = opts || {};
        const supports = opts.supports || [];
        return {
            chips: [],
            density: "comfortable",
            quickFilter: "",
            init() {
                const params = new URLSearchParams(location.search);
                this.chips = supports
                    .filter(k => params.get(k))
                    .map(k => ({ key: k, value: params.get(k) }));
                // Default is comfortable; switch to compact only when
                // the URL or saved preference says so explicitly.
                let stored = null;
                try { stored = localStorage.getItem("regpoly.density"); }
                catch (e) {}
                if (params.has("density")) {
                    this.density = params.get("density") === "compact"
                        ? "compact" : "comfortable";
                } else if (stored === "compact") {
                    this.density = "compact";
                } else {
                    this.density = "comfortable";
                }
                window.addEventListener("popstate", () => this.init());
                window.addEventListener(
                    "regpoly:filters-changed", () => this.init());
            },
            removeChip(chip) {
                const params = new URLSearchParams(location.search);
                params.delete(chip.key);
                const url = location.pathname +
                    (params.toString() ? "?" + params.toString() : "");
                history.pushState({}, "", url);
                this.init();
            },
            toggleDensity() {
                const params = new URLSearchParams(location.search);
                if (this.density === "compact") {
                    params.delete("density");
                    this.density = "comfortable";
                    try { localStorage.setItem("regpoly.density", "comfortable"); } catch (e) {}
                } else {
                    params.set("density", "compact");
                    this.density = "compact";
                    try { localStorage.setItem("regpoly.density", "compact"); } catch (e) {}
                }
                const url = location.pathname +
                    (params.toString() ? "?" + params.toString() : "");
                history.pushState({}, "", url);
                // Immediate visual switch without reload.
                document.querySelectorAll("table.table").forEach(t => {
                    if (this.density === "compact") t.classList.add("table-sm");
                    else t.classList.remove("table-sm");
                });
            },
        };
    };

    /* ---------------- Toast region ----------------
     * Mounted by base.html. Public API:
     *   regpoly.toast({title, body, type})  // type: info|success|warn|error
     *   regpoly.toast("simple message")
     */
    function toast(opts) {
        const region = document.querySelector("[data-toast-region]");
        if (!region) return;
        const o = (typeof opts === "string") ? {body: opts} : (opts || {});
        const type = o.type || "info";
        const cls = {
            info: "alert-info", success: "alert-success",
            warn: "alert-warning", error: "alert-danger",
        }[type] || "alert-info";
        const el = document.createElement("div");
        el.className = "rp-toast alert " + cls;
        el.setAttribute("role", "status");
        el.innerHTML =
            (o.title ? '<strong>' + escapeHtml(o.title) + '</strong> ' : '') +
            (o.body ? escapeHtml(o.body) : '');
        region.appendChild(el);
        setTimeout(() => {
            el.style.transition = "opacity 220ms";
            el.style.opacity = "0";
            setTimeout(() => el.remove(), 240);
        }, o.duration || 4000);
    }
    function escapeHtml(s) {
        const div = document.createElement("div");
        div.textContent = s;
        return div.innerHTML;
    }

    /* ---------------- Command palette (⌘K / Ctrl+K) ----------------
     * Lightweight built-in (no cmdk dep) — fuzzy-find across nav
     * items + the visible page's headers. Opens the <dialog
     * data-command-palette> mounted by base.html.
     */
    function buildCommandPalette() {
        const dlg = document.querySelector("dialog[data-command-palette]");
        if (!dlg) return;
        if (dlg.dataset.cmdkBuilt) return;
        dlg.dataset.cmdkBuilt = "1";
        dlg.innerHTML = `
            <div class="rp-cmdk-root" data-cmdk-root role="dialog"
                 aria-label="Command palette">
                <div class="p-3 border-bottom">
                    <input type="text" class="form-control"
                           data-cmdk-input
                           placeholder="Jump to page or generator…"
                           autocomplete="off" autofocus>
                </div>
                <div class="rp-cmdk-list" data-cmdk-list role="listbox"
                     style="max-height: 50vh; overflow: auto;"></div>
            </div>
        `;
        const items = [
            ["Dashboard", "/", "home"],
            ["Generators", "/generators", "cpu"],
            ["Combined generators", "/tested-generators", "math-function"],
            ["Searches", "/searches", "player-play"],
            ["Library", "/library", "search"],
            ["Tools", "/tools", "settings"],
            ["New primitive search", "/primitive-search", "player-play"],
            ["New tempering search", "/tempering-search", "player-play"],
        ];
        const input = dlg.querySelector("[data-cmdk-input]");
        const list = dlg.querySelector("[data-cmdk-list]");
        function render(query) {
            const q = (query || "").toLowerCase();
            list.innerHTML = "";
            const matched = items.filter(([label]) =>
                !q || label.toLowerCase().includes(q));
            matched.forEach(([label, href, icon], idx) => {
                const a = document.createElement("a");
                a.href = href;
                a.className = "list-group-item list-group-item-action " +
                    "d-flex align-items-center gap-2";
                a.setAttribute("role", "option");
                a.dataset.cmdkItem = String(idx);
                a.innerHTML =
                    `<svg width="18" height="18">` +
                    `<use href="#tabler-${icon}"/></svg>` +
                    `<span>${escapeHtml(label)}</span>`;
                list.appendChild(a);
            });
        }
        render("");
        input.addEventListener("input", () => render(input.value));
        input.addEventListener("keydown", (ev) => {
            if (ev.key === "Enter") {
                const first = list.querySelector("[data-cmdk-item]");
                if (first) { ev.preventDefault(); first.click(); }
            }
            if (ev.key === "Escape") { dlg.close(); }
        });
        dlg.addEventListener("click", (ev) => {
            if (ev.target === dlg) dlg.close();
        });
    }
    function openCommandPalette() {
        const dlg = document.querySelector("dialog[data-command-palette]");
        if (!dlg) return;
        buildCommandPalette();
        if (typeof dlg.showModal === "function") {
            try { dlg.showModal(); }
            catch (e) { dlg.setAttribute("open", ""); }
        } else {
            dlg.setAttribute("open", "");
        }
        const input = dlg.querySelector("[data-cmdk-input]");
        if (input) { input.value = ""; input.focus(); }
    }
    function bindCommandPalette() {
        document.addEventListener("keydown", (ev) => {
            if ((ev.metaKey || ev.ctrlKey) && ev.key.toLowerCase() === "k") {
                ev.preventDefault();
                openCommandPalette();
            }
        });
    }

    /* ---------------- uPlot helpers ----------------
     * Tiny wrappers over uPlot for sparklines + activity charts.
     * No-ops gracefully when uPlot hasn't loaded yet.
     */
    function sparklineUplot(el, points, opts) {
        if (typeof window.uPlot !== "function") return null;
        if (!el || !points || points.length === 0) return null;
        const xs = points.map((_, i) => i);
        const data = [xs, points.slice()];
        const o = opts || {};
        const u = new window.uPlot({
            width: o.width || 120,
            height: o.height || 32,
            scales: {x: {time: false}},
            axes: [{show: false}, {show: false}],
            legend: {show: false},
            cursor: {show: false},
            series: [
                {},
                {
                    stroke: o.color || "var(--tblr-primary, #5B5BD6)",
                    width: 1.5,
                    fill: o.fill || "rgba(91, 91, 214, 0.10)",
                },
            ],
        }, data, el);
        return u;
    }

    /* ---------------- DOMContentLoaded ---------------- */
    document.addEventListener("DOMContentLoaded", () => {
        bindThemeToggle();
        bindCopyDelegate();
        bindKeyboard();
        bindCommandPalette();
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

    // Public API for inline templates / Alpine factories.
    window.regpoly = window.regpoly || {};
    window.regpoly.toast = toast;
    window.regpoly.openCommandPalette = openCommandPalette;
    window.regpoly.sparklineUplot = sparklineUplot;
})();
