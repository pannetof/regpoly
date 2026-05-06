# REGPOLY web — user stories

This document enumerates the concrete user-facing scenarios the
`regpoly-web` application supports as of v2.0.0. It is the
operator-side companion to [Web usage](web.md): if `web.md` answers
"how do I run the server?", this file answers "what can a user
*do* once it is running?". Each story is grounded in an actual
route or template, not aspirational.

The stories are grouped by **persona**:

- **R** — Researcher: looking for new high-quality F₂-linear PRNG
  parameter sets, wants to launch and monitor searches.
- **C** — Cataloger: maintains the published-generator library;
  imports parameter sets from external sources; publishes
  search results worth keeping.
- **A** — Analyst: browses already-tested generators, compares
  equidistribution profiles, exports for downstream use.
- **D** — Developer: integrates the catalog or test runner into
  another tool; needs the JSON API surface.

---

## 1. Onboarding & navigation

### S1.1 — First-visit landing (any persona)

**As any user**, when I navigate to `/`, **I see** a dashboard
listing every supported generator family grouped by category
(MT-derived, WELL, MELG, Tausworthe, …) and the catalog of
published configurations grouped by paper, **so that** I can
orient myself without prior knowledge of the codebase.

*Route:* `GET /` → `templates/index.html`. *Backed by:*
`/api/families`, `/api/library/papers`.

### S1.2 — Family detail page (R, A)

**As a researcher or analyst**, when I click a family name on the
landing page, **I see** the family's mathematical role
(structural vs search parameters, period exponent formula, links
to the matching theory page), **so that** I can decide whether
this family is a fit for my application.

*Route:* `GET /family/{family}` → `templates/family/landing.html`.
*Backed by:* `/api/families/{family}`, `/api/families/{family}/docs`.

### S1.3 — Persistent global nav (any persona)

**As any user**, every page exposes the same global navigation
(Dashboard, Generators, Combined-generator search, Searches,
Combined generators, Library), **so that** I can move between
the four main workflows in one click.

*Implementation:* `templates/base.html`.

---

## 2. Browsing the live generator pool (A)

### S2.1 — List every full-period generator on file

**As an analyst**, when I open `/generators`, **I see** a
filterable, paginated table of every primitive (full-period)
generator stored in the SQLite DB, with columns for family, $L$,
$k$, structural params, found-at-try, and the analysis snapshot
(characteristic polynomial hex, Hamming weight, $\sum \delta(\ell)$).

*Route:* `GET /generators` → `templates/generators/list.html`.
*Backed by:* `/api/generators` with query filters.

### S2.2 — Filter by family / parameter values

**As an analyst**, the list page lets me filter by family,
$k$ range, and any structural parameter value, **so that** I can
narrow to the candidates that match an external reference
implementation.

*Backed by:* `/api/generators/families`, `/api/generators/param-values`.

### S2.3 — Drill into a single generator

**As an analyst**, clicking a row opens a detail page showing
the full parameter set, the equidistribution analysis result if
computed, the source search-run id, and a button to fetch the
transition-matrix coordinates for visual inspection.

*Route:* `GET /generators/{gen_id}` →
`templates/generators/detail.html`. *Backed by:*
`/api/generators/{gen_id}`,
`/api/generators/{gen_id}/transition-matrix-coords`.

### S2.4 — Delete a generator

**As an analyst**, I can delete a generator I no longer need
from its detail page, **so that** the analysis worker doesn't
keep retrying a known-bad row in a tight loop.

*Backed by:* `DELETE /api/generators/{gen_id}`.

---

## 3. Launching a primitive (full-period) search (R)

### S3.1 — Start a new full-period search

**As a researcher**, I navigate to `/primitive-search`, choose a
family, set the structural params (e.g. MT $w=32, r=624, p=31$),
and pick a search budget (`max_tries`, `max_seconds`), **so that**
I can hunt for new full-period parameter sets without writing a
YAML config by hand.

*Route:* `GET /primitive-search` →
`templates/primitive_search/form.html`. *Backed by:* `POST
/api/primitive-searches` (returns the new run id).

### S3.2 — Estimate before launching

**As a researcher**, before I commit to a long search budget I
can ask the server to estimate how many tries the parameter
space contains in exhaustive mode, **so that** I pick `random`
vs `exhaustive` mode wisely.

*Backed by:* `POST /api/primitive-searches/estimate`.

### S3.3 — Watch progress in real time

**As a researcher**, after kicking off a search the detail page
shows a live progress bar (tries done / max tries / found count)
that updates over Server-Sent Events, **so that** I don't have
to refresh manually.

*Route:* `GET /primitive-search/{run_id}` →
`templates/primitive_search/detail.html`. *Backed by:*
`/api/primitive-searches/{run_id}/progress` (SSE).

### S3.4 — Pause / resume / cancel / restart

**As a researcher**, I can pause an ongoing search, resume it
later, cancel it permanently, or restart it from scratch with
the same config, **so that** I can yield the workers to a
higher-priority job and pick this one back up later.

*Backed by:* `POST /api/primitive-searches/{run_id}/{pause,resume,cancel,restart}`.

### S3.5 — Inspect the generators a run produced

**As a researcher**, the search detail page lists every
generator that the run accepted (each row a primitive generator),
**so that** I can pick the most promising ones for the tempering
search.

*Backed by:* `GET /api/primitive-searches/{run_id}/generators`.

### S3.6 — Bulk-delete uncommitted hits

**As a researcher**, I can bulk-delete every generator a run
produced in one click, **so that** I can clean up an exploratory
run that turned out badly.

*Backed by:* `DELETE /api/primitive-searches/{run_id}/generators`.

---

## 4. Launching a tempering search (R)

### S4.1 — Start a combined-generator (tempering) search

**As a researcher**, I navigate to `/tempering-search`, pick
one or more components (each component is a pool of primitive
generators + a tempering parameter shape), choose a test type
(equidistribution / collision_free / tuplets) and budget
(`nb_tries`, optional `optimizer_config`), **so that** I can
search for a combined generator that passes the chosen test.

*Route:* `GET /tempering-search` →
`templates/tempering_search/form.html`. *Backed by:* `POST
/api/tempering-searches`.

### S4.2 — Live progress per component

**As a researcher**, the detail page shows the live state of
each combo (combos_done / combos_total) and the best-so-far
sum of equidistribution gaps, streamed over SSE, **so that** I
can decide whether to let it run to completion or kill it
early.

*Route:* `GET /tempering-search/{run_id}` →
`templates/tempering_search/detail.html`. *Backed by:*
`/api/tempering-searches/{run_id}/progress` (SSE).

### S4.3 — Pause / resume / cancel / restart (tempering)

**As a researcher**, the same lifecycle controls available for
primitive searches apply: pause, resume, cancel, restart with
the same config.

*Backed by:* `POST /api/tempering-searches/{run_id}/{pause,resume,cancel,restart}`.

### S4.4 — Inspect every tested generator the run produced

**As a researcher**, I can browse every combined generator the
search produced, sorted by best score, **so that** I can pick
the winner and publish it.

*Backed by:* `GET /api/tempering-searches/{run_id}/results`.

### S4.5 — Bulk-delete a run's results

**As a researcher**, I can purge every tested generator a
tempering run produced in one click, **so that** I can clear
the noise from an exploratory run before the next attempt.

*Backed by:* `DELETE /api/tempering-searches/{run_id}/results`.

---

## 5. Tracking searches across the lifetime of the server (R, A)

### S5.1 — Unified searches dashboard

**As a researcher or analyst**, `/searches` lists every
full-period and combined-generator search the server has ever
seen, with status (running / pending / paused / completed /
cancelled / failed), elapsed time, and per-row inline
pause/resume/cancel/restart buttons, **so that** I have a
single page from which to triage active and historical work.

*Route:* `GET /searches` → `templates/searches.html`.
*Backed by:* `/api/primitive-searches`, `/api/tempering-searches`.

### S5.2 — Filter the searches list

**As any user**, the dashboard accepts query-string filters by
type (full-period / combined), family, $k$, status, and offers a
configurable refresh interval (2s / 5s / 10s / off), **so that**
I can keep an eye on long runs without the page flapping.

### S5.3 — Searches are durable across server restarts

**As an operator**, when I restart the server, any search that
was `running` or `pending` at shutdown is reaped to `cancelled`
with an explanatory `error_message`, so the dashboard never
shows phantom "running" rows that no worker is actually
processing.

*Implementation:* `regpoly_web.database.init_sync` runs an
orphan-reap on startup.

---

## 6. Browsing tested generators (A)

### S6.1 — List every tested generator

**As an analyst**, `/tested-generators` lists every combined
generator that has been tested (single-component or
multi-component), filterable by family, $k_g$ range, test type,
maximum allowed $\sum \delta$, and "is ME?", **so that** I can
find prior work matching the shape I need.

*Route:* `GET /tested-generators` →
`templates/tested_generators/list.html`. *Backed by:*
`/api/tested-generators` with query filters.

### S6.2 — Detail page for a tested generator

**As an analyst**, clicking a row opens a detail page showing
the components (family, $L$, $k$, structural params, tempering
chain), the per-test results (equidistribution / collision_free
/ tuplets), and the link back to the source search run, **so
that** I can audit a published reference implementation without
leaving the browser.

*Route:* `GET /tested-generators/{tg_id}` →
`templates/tested_generators/detail.html`. *Backed by:*
`/api/tested-generators/{tg_id}` — which since v2.0 reads from
the typed result tables introduced in [schema v2](web.md#database-schema).

### S6.3 — Delete a tested generator

**As an analyst**, I can delete a tested generator from its
detail page (cascade-deletes its components and per-test
results), **so that** I can prune the DB.

*Backed by:* `DELETE /api/tested-generators/{tg_id}`.

### S6.4 — Detail page surfaces the typed equidistribution result

**As an analyst**, the detail page renders the sparse
$\delta(\ell)$ map (only the resolutions where the gap is
non-zero), the verified flag, and the elapsed analysis time,
**so that** I can compare against published reference numbers
without staring at JSON.

*Backed by:* the v2 typed table read path
(`regpoly_web.results.read_typed_results_async`).

---

## 7. Library / catalog browsing (C, A)

### S7.1 — Browse all published parameter sets

**As a cataloger or analyst**, `/library` lists every published
generator across every paper in the catalog, filterable by
family, **so that** I can find a known-good MT / WELL / MELG
parametrisation by name.

*Route:* `GET /library` → `templates/library/list.html`.
*Backed by:* `/api/library/generators`.

### S7.2 — Browse by paper

**As a cataloger**, `/library/papers` lists every paper in the
catalog (Matsumoto & Nishimura 1998, Saito & Matsumoto 2008,
Harase 2014, …), and `/library/{paper_id}` shows the paper's
metadata and the parameter sets it introduced, **so that** I
can navigate the catalog the way the literature is organised.

*Routes:* `GET /library/papers` → `templates/library/papers.html`,
`GET /library/{paper_id}` → `templates/library/paper.html`.
*Backed by:* `/api/library/papers`, `/api/library/papers/{paper_id}`.

### S7.3 — Per-generator library page

**As a cataloger**, `/library/{paper_id}/{gen_id}` shows one
published parameter set in full — components, tempering chain,
notes from the YAML — **so that** I can reproduce or cite it.

*Route:* `GET /library/{paper_id}/{gen_id}` →
`templates/library/detail.html`. *Backed by:*
`/api/library/generators/{gen_id}`.

### S7.4 — Instantiate a library entry

**As a cataloger or analyst**, from a library detail page I can
click "instantiate" to create a `tested_generator` row in the
DB pointing at the library entry's parameters, **so that** the
analysis machinery can test it without retyping the params.

*Backed by:* `POST /api/library/generators/{gen_id}/instantiate`.

### S7.5 — "Run test" on a library entry without leaving the page

**As an analyst**, the library detail page exposes a "Run test"
action that kicks off a one-shot equidistribution test against
the library entry's parameters in a background worker, with a
job id I poll for completion, **so that** I get a quick
verification number without going through the full search UI.

*Backed by:* `POST /api/library/generators/{gen_id}/run-test`,
polled via `GET /api/library/run-test-jobs/{job_id}`.

---

## 8. Importing & exporting parameter sets (C)

### S8.1 — Import a single YAML file of generator parameters

**As a cataloger**, I can upload a YAML file describing one or
more generator parameter sets through `POST /api/import/generators`,
**so that** I can grow the DB from external reference
implementations without writing per-row SQL.

*Backed by:* `POST /api/import/generators` (multipart upload).

### S8.2 — Import every YAML in a directory

**As a cataloger** doing a bulk migration, I can ask the server
to import every YAML under a server-local directory, **so that**
I don't upload one file at a time.

*Backed by:* `POST /api/import/generators-dir` (server-side path).

### S8.3 — Export the live generator pool to YAML

**As a cataloger**, I can export every generator currently in
the DB as a YAML stream via `GET /api/export/generators`, **so
that** I can move parameter sets between two REGPOLY instances
or check them into a sibling repo.

*Backed by:* `GET /api/export/generators` (text/plain YAML).

### S8.4 — Imports are tracked in the DB

**As an operator**, every import is recorded with the source
filename, type, row count, and timestamp in the `yaml_import`
table, **so that** I can audit when a particular parameter set
landed.

---

## 9. Publishing (C)

### S9.1 — Mark a tested generator as published

**As a cataloger**, when a tested generator passes my
acceptance bar, I can attach a `library_id` to it, **so that**
the web UI surfaces the "published" marker and the generator
shows up in catalog browses.

*Note:* publishing **to** a paper YAML (writing the
`generators:` block) is intentionally a C++-only operation per
[Phase 4.3](../theory/search_format.md). The web app reads the
flag via `tested_generator.library_id`; the canonical
write-side is `regpoly-cli publish FILE.yaml --paper PAPER_ID
--gen-id GEN_ID`.

### S9.2 — Unpublish a tested generator

**As a cataloger**, I can clear the `library_id` to retract a
publication if the parameter set turns out to be incorrect.

---

## 10. JSON API for downstream tooling (D)

### S10.1 — Family introspection endpoint

**As a developer** integrating REGPOLY into another tool,
`GET /api/families` returns every supported family with its
canonical name, parameter schema (structural vs search,
random-sampler type), and "enumerable" flag, **so that** I can
build a UI on top without hard-coding family names.

*Backed by:* `regpoly.introspection.get_gen_param_specs`.

### S10.2 — Transformation introspection

**As a developer**, `GET /api/transformations` lists every
tempering transformation type (`tempMK`, `tempMK2`, `permut`,
`laggedTempering`, …) with its parameter schema, **so that** I
can build a tempering-config form that stays in sync with the
C++ catalog.

### S10.3 — Test catalog

**As a developer**, `GET /api/tests` lists the three test
kinds (`equidistribution`, `collision_free`, `tuplets`) the
search drivers support, **so that** I can render a "pick a
test" UI that doesn't drift from the backend.

### S10.4 — Stable list+detail+delete pattern

**As a developer**, every entity (generators, primitive
searches, tempering searches, tested generators) follows the
same `GET /api/{entity}` (list) / `GET /api/{entity}/{id}`
(detail) / `DELETE /api/{entity}/{id}` shape, **so that** I
can write a single CRUD client.

### S10.5 — Server-Sent Events for live progress

**As a developer**, the per-run `progress` endpoints return
`text/event-stream`, with one event per write to
`search_progress`, **so that** I can integrate live progress
into another dashboard without polling.

*Endpoints:* `GET /api/primitive-searches/{run_id}/progress`,
`GET /api/tempering-searches/{run_id}/progress`.

---

## 11. Operational characteristics (any persona)

### S11.1 — Single-user, single-machine assumption

The web app does not authenticate; it expects to run on
`127.0.0.1` or behind a reverse proxy that handles authentication.
This is intentional given REGPOLY's research/academic use case.

### S11.2 — DB upgrades happen automatically

**As an operator**, the first time I open an older v1 SQLite
DB with a v2.0 server, the schema is upgraded in place at
startup — typed result tables are created and the legacy
`test_result` rows are backfilled. The migration is idempotent.

*See:* [Web usage → Database schema](web.md#database-schema).

### S11.3 — Workers are out-of-process

**As an operator**, every search runs in a `ProcessPoolExecutor`
worker so the FastAPI request-handling thread never blocks on a
multi-second analysis. The pool size is configured via the
`Settings.pool_size` field.

### S11.4 — Search runs survive a server crash gracefully

**As an operator**, if the server is killed mid-search, the
next startup reaps the orphaned `running`/`pending` rows to
`cancelled` (with an explanatory `error_message`) — there are
no zombie progress rows that pretend a worker is still alive.

---

## See also

- [Web usage](web.md) — operator-side: install, flags, schema migration.
- [Architecture](../dev/architecture.md) — how the web app talks to
  the C++ core through `regpoly`.
- [Python usage](python.md) — the same algorithms without the web layer.
