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

/* ───────────── Shared client-side table sort / filter helpers ───────────
 * Apply to any list of row objects.  Per-column filters do a
 * case-insensitive substring match on the value returned by `getValue`.
 * Sorting detects numeric values and compares numerically when both
 * sides parse as numbers, falling back to locale-aware string compare.
 */

/* Canonical string key used for both filter matching and distinct-list
 * population.  Arrays are serialised with JSON.stringify so two equal
 * lists hash to the same key. */
function _rngKey(v) {
    if (v === null || v === undefined) return '';
    if (Array.isArray(v)) return JSON.stringify(v);
    return String(v);
}

window.rngTableApply = function (rows, sort, filters, getValue) {
    let out = rows || [];
    if (filters) {
        for (const [key, allowed] of Object.entries(filters)) {
            // Only Set-typed filters are active; anything else is
            // treated as "no filter on this column".
            if (!(allowed instanceof Set)) continue;
            out = out.filter(r => allowed.has(_rngKey(getValue(r, key))));
        }
    }
    if (sort && sort.key && sort.dir) {
        const key = sort.key;
        const cmp = (a, b) => {
            const va = getValue(a, key);
            const vb = getValue(b, key);
            const na = Number(va);
            const nb = Number(vb);
            const aNum = va !== '' && va !== null && va !== undefined && !isNaN(na);
            const bNum = vb !== '' && vb !== null && vb !== undefined && !isNaN(nb);
            if (aNum && bNum) return na - nb;
            return String(va == null ? '' : va)
                .localeCompare(String(vb == null ? '' : vb));
        };
        out = [...out].sort(cmp);
        if (sort.dir === 'desc') out.reverse();
    }
    return out;
};

/* Return every distinct value a column takes, in sort order.  Each
 * entry is {value: canonical-string-key, display: original-value}.
 * Used to populate the per-column filter dropdown. */
window.rngTableDistinct = function (rows, key, getValue) {
    const seen = new Map();
    for (const r of rows || []) {
        const v = getValue(r, key);
        const k = _rngKey(v);
        if (!seen.has(k)) seen.set(k, v);
    }
    const out = Array.from(seen.entries()).map(
        ([value, display]) => ({ value, display })
    );
    out.sort((a, b) => {
        const na = Number(a.value), nb = Number(b.value);
        const isA = a.value !== '' && !isNaN(na);
        const isB = b.value !== '' && !isNaN(nb);
        if (isA && isB) return na - nb;
        return a.value.localeCompare(b.value);
    });
    return out;
};

/* Shared filter-state behaviour for Alpine components that render a
 * tabular view.  Spread the return value into the component's state
 * object and implement two helpers on the component itself:
 *   rowsForFilter()        — source rows used to enumerate distinct
 *                            column values (typically `this.items`).
 *   rowValue(row, key)     — getter shared with rngTableApply. */
window.rngTableFilterMethods = function () {
    return {
        openFilter: null,          // column key whose popover is open
        _filterAnchor: null,       // DOM element the popover is anchored to

        toggleFilter(key, ev) {
            if (ev) ev.stopPropagation();
            if (this.openFilter === key) {
                this.openFilter = null;
                this._filterAnchor = null;
            } else {
                this.openFilter = key;
                this._filterAnchor =
                    (ev && ev.currentTarget) ? ev.currentTarget : null;
            }
        },

        filterPopoverStyle() {
            const el = this._filterAnchor;
            if (!el || typeof el.getBoundingClientRect !== 'function') return '';
            const r = el.getBoundingClientRect();
            const top  = r.bottom + window.scrollY;
            const left = r.left   + window.scrollX;
            return `position:absolute;top:${top}px;left:${left}px;`;
        },

        _distinctValues(key) {
            return rngTableDistinct(
                this.rowsForFilter(), key,
                (r, k) => this.rowValue(r, k)).map(d => d.value);
        },

        filterDistinctDisplay(key) {
            return rngTableDistinct(
                this.rowsForFilter(), key,
                (r, k) => this.rowValue(r, k));
        },

        filterActive(key) {
            const s = this.colFilters[key];
            if (!(s instanceof Set)) return false;
            const all = this._distinctValues(key);
            return all.some(v => !s.has(v)) || s.size === 0;
        },

        allChecked(key) {
            return !(this.colFilters[key] instanceof Set);
        },

        isChecked(key, value) {
            const s = this.colFilters[key];
            if (!(s instanceof Set)) return true;
            return s.has(value);
        },

        toggleAll(key) {
            if (this.allChecked(key)) {
                // Currently "all" → switch to "none" (filter hides everything).
                this.colFilters = { ...this.colFilters, [key]: new Set() };
            } else {
                // Any other state → back to "all".
                const cf = { ...this.colFilters };
                delete cf[key];
                this.colFilters = cf;
            }
        },

        toggleOne(key, value) {
            let s = this.colFilters[key];
            const all = this._distinctValues(key);
            if (!(s instanceof Set)) {
                // First tick → everything but the toggled value is allowed.
                s = new Set(all);
                s.delete(value);
            } else {
                if (s.has(value)) s.delete(value);
                else              s.add(value);
            }
            // If every distinct value is allowed again, drop the filter.
            if (s.size === all.length && all.every(v => s.has(v))) {
                const cf = { ...this.colFilters };
                delete cf[key];
                this.colFilters = cf;
            } else {
                this.colFilters = { ...this.colFilters, [key]: s };
            }
        },
    };
};

window.rngTableToggleSort = function (sort, key) {
    const s = sort || {};
    if (s.key !== key) return { key, dir: 'asc' };
    if (s.dir === 'asc')  return { key, dir: 'desc' };
    return { key: null, dir: null };
};

window.rngTableSortArrow = function (sort, key) {
    const s = sort || {};
    if (s.key !== key) return '⇅';
    if (s.dir === 'asc')  return '▲';
    if (s.dir === 'desc') return '▼';
    return '⇅';
};
