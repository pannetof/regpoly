# regpoly-web — UI redesign spec

## Context

This spec is for Claude Code. The goal is to redesign the regpoly-web frontend from
plain HTML/CSS to a polished, professional dashboard interface using the Tabler
open-source template. The app is a single-user scientific research tool for studying
GF(2)-linear pseudorandom number generators.

## Current stack (do NOT change)

- **Backend**: FastAPI + Uvicorn + Jinja2 templates
- **Frontend**: server-rendered HTML, Alpine.js for interactivity, static CSS
- **Data**: SQLite (aiosqlite async + sqlite3 for workers)
- **Live updates**: SSE (Server-Sent Events) for search progress
- **Heavy compute**: ProcessPoolExecutor → C++ binaries
- **No build step**: no bundler, no npm build, no Tailwind — static files served as-is

## Design framework: Tabler via CDN

Add these two lines to the base Jinja2 template `<head>`:

```html
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/@tabler/core@latest/dist/css/tabler.min.css">
<script src="https://cdn.jsdelivr.net/npm/@tabler/core@latest/dist/js/tabler.min.js" defer></script>
```

Alpine.js must still be loaded AFTER Tabler JS. No jQuery.

---

## Design references

The target aesthetic draws from these tools (screenshots provided to the designer):

| Reference        | What to take from it                                             |
|------------------|------------------------------------------------------------------|
| Grafana          | Dark sidebar, card-based content area, dark mode overall         |
| Cloudflare       | Clean light mode, metric summary cards, breadcrumb nav, spacing  |
| Smartlook        | Data-dense tables with row hover, left sidebar filters           |
| Wolfram Alpha    | Scientific input/output layout, section-based results, math font |
| Modernize        | Sidebar nav with icons, metric cards grid, table styling         |

---

## Global layout

Use Tabler's **vertical layout** (sidebar left + main content right):

```
┌──────────┬─────────────────────────────────────┐
│          │  Breadcrumb: Home > Section > Page   │
│  SIDEBAR │─────────────────────────────────────│
│          │                                     │
│  • Home  │  Page title                         │
│  • Gen.  │  Subtitle / description             │
│  • Search│                                     │
│  • ...   │  ┌─ Card ─────┐  ┌─ Card ─────┐   │
│          │  │ Metric      │  │ Metric      │   │
│          │  └─────────────┘  └─────────────┘   │
│          │                                     │
│          │  ┌─ Card ──────────────────────────┐│
│          │  │ Data table / results             ││
│          │  └─────────────────────────────────┘│
└──────────┴─────────────────────────────────────┘
```

### Sidebar navigation

The sidebar must:
- Be collapsible (Tabler built-in)
- Show the app name "regpoly" and a small logo/icon at the top
- Use Tabler icons (loaded with the CDN) for each nav item
- Highlight the current page

Navigation items (adapt to actual app routes):

| Label                | Icon suggestion          | Route                  |
|----------------------|--------------------------|------------------------|
| Dashboard            | `icon-home`              | `/`                    |
| Tested generators    | `icon-cpu`               | `/tested-generators/`  |
| Tempering search     | `icon-search`            | `/tempering-search/`   |
| Polynomial tools     | `icon-math-function`     | `/polynomial/`         |
| Running jobs         | `icon-player-play`       | `/jobs/`               |
| Settings             | `icon-settings`          | `/settings/`           |

*Adjust items and routes to match actual app structure.*

### Breadcrumbs

Every page must show a breadcrumb trail in the page header area:
`Home > Tested generators > MELG607-64`

Use Tabler's `.page-header` + `.breadcrumb` component.

---

## Theme and colors

### Default: dark mode

The primary mode is **dark** (like Grafana). Use Tabler's built-in dark theme:

```html
<body class="theme-dark">
```

Optionally add a light/dark toggle in the navbar (Tabler has a built-in switcher component).

### Accent color

Use Tabler's **azure** (default blue) as the primary accent. This matches the
Cloudflare/Grafana aesthetic without customization effort.

For status indicators:
- **Primitive**: green badge (`badge bg-success`)
- **Not primitive**: red badge (`badge bg-danger`)
- **Running**: blue badge with spinner (`badge bg-primary`)
- **Pending**: yellow badge (`badge bg-warning`)

---

## Component patterns

### 1. Metric cards (generator summary)

At the top of detail pages, show key properties as a row of metric cards.
Use Tabler's `card` with `card-sm`:

```html
<div class="row row-deck row-cards mb-3">
  <div class="col-sm-6 col-lg-3">
    <div class="card card-sm">
      <div class="card-body">
        <div class="row align-items-center">
          <div class="col-auto">
            <span class="bg-primary text-white avatar">k</span>
          </div>
          <div class="col">
            <div class="font-weight-medium">Period exponent</div>
            <div class="text-secondary">607</div>
          </div>
        </div>
      </div>
    </div>
  </div>
  <!-- more cards: word size, primitivity status, D_max -->
</div>
```

Generator detail pages should have metric cards for:
- Period exponent (k)
- Word size (w)
- Primitivity status (badge)
- Best known D_max
- Number of search results

### 2. Data tables (results, generator lists)

Use Tabler's `.table .table-vcenter .card-table`:

```html
<div class="card">
  <div class="card-header">
    <h3 class="card-title">Search results</h3>
    <div class="card-actions">
      <!-- filter controls with Alpine.js -->
    </div>
  </div>
  <div class="table-responsive">
    <table class="table table-vcenter card-table table-hover">
      <thead>
        <tr>
          <th>Rank</th>
          <th>Parameters</th>
          <th>D_max</th>
          <th>Merit</th>
          <th>Time</th>
        </tr>
      </thead>
      <tbody>
        <!-- rows -->
      </tbody>
    </table>
  </div>
</div>
```

**Important for this app**:
- Hex parameter values must use `font-family: var(--tblr-font-monospace)` (monospace)
- Tables may have many rows — add `table-hover` for row highlighting
- Consider Tabler's `.table-sort` for client-side sorting via Alpine.js
- For very long tables, keep the header sticky if Tabler supports it, or add via CSS:
  `thead th { position: sticky; top: 0; z-index: 1; }`

### 3. Forms (search configuration)

Use Tabler's form components inside a card:

```html
<div class="card">
  <div class="card-header">
    <h3 class="card-title">Tempering search configuration</h3>
  </div>
  <div class="card-body">
    <div class="row mb-3">
      <div class="col-md-4">
        <label class="form-label">Generator</label>
        <select class="form-select" x-model="generator">
          <!-- options -->
        </select>
      </div>
      <div class="col-md-4">
        <label class="form-label">Tempering parameter b₁</label>
        <input type="text" class="form-control font-monospace"
               x-model="b1" placeholder="0x...">
      </div>
      <div class="col-md-4">
        <label class="form-label">Search depth D_max</label>
        <input type="number" class="form-control"
               x-model="dmax" min="1" max="64">
      </div>
    </div>
  </div>
  <div class="card-footer text-end">
    <button class="btn btn-primary" @click="runSearch">
      <svg><!-- tabler icon: play --></svg>
      Run search
    </button>
    <button class="btn btn-outline-secondary" @click="cancel"
            x-show="isRunning">
      Cancel
    </button>
  </div>
</div>
```

**Hex input fields** must always use `font-monospace` class.

### 4. Progress / SSE streaming

During a running search, show progress inside a card:

```html
<div class="card" x-show="isRunning">
  <div class="card-body">
    <div class="d-flex align-items-center mb-2">
      <div class="spinner-border spinner-border-sm text-primary me-2"></div>
      <span>Searching... <span x-text="progress"></span></span>
    </div>
    <div class="progress">
      <div class="progress-bar" :style="`width: ${pct}%`"
           x-text="`${pct}%`"></div>
    </div>
    <div class="mt-2 text-secondary small" x-text="statusMessage"></div>
  </div>
</div>
```

The SSE events should update Alpine.js state, which in turn updates these components.
No changes needed to the SSE backend — only the frontend presentation changes.

### 5. Empty states

When a page has no data yet, use Tabler's empty state pattern:

```html
<div class="empty">
  <p class="empty-title">No search results yet</p>
  <p class="empty-subtitle text-secondary">
    Configure and run a tempering search to see results here.
  </p>
  <div class="empty-action">
    <a href="/tempering-search/" class="btn btn-primary">
      Start a search
    </a>
  </div>
</div>
```

---

## Typography rules

| Element                          | Style                                      |
|----------------------------------|-------------------------------------------|
| Page titles                      | Tabler default `h1` / `.page-title`       |
| Card titles                      | `h3.card-title` (Tabler default)          |
| Body text                        | Tabler default sans-serif                  |
| Hex values (params, addresses)   | `.font-monospace`, same size as body       |
| Mathematical notation            | If needed, load KaTeX via CDN for formulas |
| Table data                       | 13-14px, monospace for hex, regular for labels |

### KaTeX (optional, for polynomial display)

If any page needs to render polynomials like `x^607 + x^273 + 1`:

```html
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16/dist/katex.min.css">
<script src="https://cdn.jsdelivr.net/npm/katex@0.16/dist/katex.min.js" defer></script>
<script src="https://cdn.jsdelivr.net/npm/katex@0.16/dist/contrib/auto-render.min.js" defer
        onload="renderMathInElement(document.body)"></script>
```

Wrap math in `\( ... \)` for inline or `$$ ... $$` for display.

---

## Page-by-page spec

### Dashboard (`/`)

- Row of 4 metric cards: total generators tested, total searches run, best global merit found, active jobs
- Card with recent activity table (last 10 searches)
- Card with quick links to common actions

### Tested generators list (`/tested-generators/`)

- Card with data table: name, k, w, primitivity badge, best D_max, number of searches
- Click a row → navigate to detail page
- Optional: filter/search bar in card header (Alpine.js)

### Generator detail (`/tested-generators/<id>/`)

- Row of metric cards (k, w, primitivity, D_max)
- Card: generator parameters (displayed as a definition list or horizontal key-value pairs)
- Card: characteristic polynomial (rendered with KaTeX if available)
- Card: search results table for this generator
- Card: tempering parameters (if applicable)

### Tempering search form (`/tempering-search/`)

- Card: search configuration form (see §3 above)
- Card: progress indicator (see §4 above, shown when running)
- Card: results table (populated as results stream in via SSE)

### Running jobs (`/jobs/`)

- Card with table: job ID, generator, status badge, progress bar, elapsed time, cancel button
- Auto-refresh via SSE or polling with Alpine.js

---

## Custom CSS overrides

Create a file `static/css/regpoly.css` loaded AFTER Tabler, for app-specific tweaks:

```css
/* Monospace for hex values everywhere */
.hex-value {
  font-family: var(--tblr-font-monospace);
}

/* Sticky table headers for long result tables */
.table-sticky thead th {
  position: sticky;
  top: 0;
  z-index: 1;
  background: var(--tblr-card-bg);
}

/* Compact metric cards for generator detail pages */
.card-metric .card-body {
  padding: 0.75rem 1rem;
}

/* SSE log output area */
.sse-log {
  font-family: var(--tblr-font-monospace);
  font-size: 0.8125rem;
  max-height: 300px;
  overflow-y: auto;
  background: var(--tblr-bg-surface);
  padding: 0.75rem;
  border-radius: var(--tblr-border-radius);
}
```

---

## Jinja2 template structure

Restructure templates to use Tabler's layout system:

```
templates/
├── base.html              ← Tabler shell: <html>, sidebar, scripts
├── partials/
│   ├── sidebar.html       ← Sidebar nav (included in base.html)
│   ├── breadcrumb.html    ← Breadcrumb macro
│   ├── metric_card.html   ← Reusable metric card macro
│   └── flash_messages.html← Tabler alert-based flash messages
├── dashboard/
│   └── index.html
├── tested_generators/
│   ├── list.html
│   └── detail.html
├── tempering_search/
│   ├── form.html
│   └── results.html
└── jobs/
    └── list.html
```

### base.html skeleton

```html
<!DOCTYPE html>
<html lang="en" data-bs-theme="dark">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{% block title %}regpoly{% endblock %}</title>
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/@tabler/core@latest/dist/css/tabler.min.css">
  <link rel="stylesheet" href="{{ url_for('static', path='css/regpoly.css') }}">
  {% block extra_head %}{% endblock %}
</head>
<body class="theme-dark">
  <div class="page">
    <!-- Sidebar -->
    <aside class="navbar navbar-vertical navbar-expand-lg">
      <div class="container-fluid">
        <h1 class="navbar-brand">
          <span class="navbar-brand-text">regpoly</span>
        </h1>
        {% include "partials/sidebar.html" %}
      </div>
    </aside>

    <div class="page-wrapper">
      <!-- Page header with breadcrumbs -->
      <div class="page-header d-print-none">
        <div class="container-xl">
          <div class="page-pretitle">
            {% block pretitle %}{% endblock %}
          </div>
          <h2 class="page-title">
            {% block page_title %}{% endblock %}
          </h2>
        </div>
      </div>

      <!-- Page body -->
      <div class="page-body">
        <div class="container-xl">
          {% block content %}{% endblock %}
        </div>
      </div>
    </div>
  </div>

  <script src="https://cdn.jsdelivr.net/npm/@tabler/core@latest/dist/js/tabler.min.js"></script>
  <script src="https://cdn.jsdelivr.net/npm/alpinejs@3/dist/cdn.min.js" defer></script>
  {% block extra_scripts %}{% endblock %}
</body>
</html>
```

---

## Migration strategy

1. **Start with base.html**: set up the Tabler shell, sidebar, CDN links
2. **Convert one page at a time**: start with the simplest page (dashboard or generator list)
3. **Keep Alpine.js bindings intact**: only change the HTML structure and CSS classes around them
4. **Test SSE integration**: make sure progress/streaming still works after the template change
5. **Add polish last**: sticky headers, KaTeX, empty states, hover effects

Each page conversion should be a self-contained commit. Don't try to convert everything at once.

---

## What NOT to change

- Backend routes, API endpoints, data models — unchanged
- Alpine.js logic (x-data, x-model, @click handlers) — keep all, only adjust surrounding HTML
- SSE event format — unchanged
- SQLite schema — unchanged
- C++ binaries / Python wrappers — unchanged
- Static file serving approach — still plain files, no bundler
