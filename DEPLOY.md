# Deploying regpoly

Three-container Docker stack. One command to bring it up. Runs on any host that has Docker — Intel mini-PC, ARM NAS (Synology/QNAP/Asustor), Raspberry Pi 5, Linux server, macOS dev box.

## TL;DR

```sh
git clone --recursive https://github.com/pannetof/regpoly.git
cd regpoly
cp deploy/.env.example deploy/.env
$EDITOR deploy/.env          # set DB_PASSWORD
chmod 600 deploy/.env
docker compose -f deploy/compose.yml up -d
```

Open `http://localhost:8000` on the host. Done.

To update: `docker compose -f deploy/compose.yml pull && docker compose -f deploy/compose.yml up -d`.

---

## What runs

Three containers on one Docker network.

| Service  | Image                                        | Role                                                         |
| -------- | -------------------------------------------- | ------------------------------------------------------------ |
| `db`     | `postgres:16.6-alpine`                       | PostgreSQL 16 server. Persists state in the `pgdata` volume. |
| `web`    | `ghcr.io/pannetof/regpoly:latest`            | FastAPI app on `:8000`. Applies DB migrations on startup.    |
| `worker` | `ghcr.io/pannetof/regpoly:latest` (same tag) | Search + analysis schedulers. Reaps orphan rows at startup.  |

`web` and `worker` share one image, run with different `command:` per service. The image is published multi-arch (`linux/amd64`, `linux/arm64`); `docker compose pull` picks the right one for your host.

Startup ordering is `db (healthy) → web (healthy) → worker (started)` via `depends_on`. `web` is "healthy" only after migrations apply, so the worker never sees a pre-migration schema.

---

## Configuration

All config lives in `deploy/.env`. Only one variable is required:

```
DB_PASSWORD=replace-me
```

Optional variables:

| Variable        | Default        | Purpose                                                                            |
| --------------- | -------------- | ---------------------------------------------------------------------------------- |
| `WEB_BIND_ADDR` | `127.0.0.1`    | Host interface the web service binds to (see [Access](#access) below).             |
| `REGPOLY_TAG`   | `latest`       | Image tag. Pin to `sha-<short>` to roll back to a specific commit.                 |
| `POSTGRES_TAG`  | `16.6-alpine`  | Postgres image tag. Override only when validated.                                  |

> ⚠️ `chmod 600 deploy/.env` after writing — `docker compose config` will otherwise print the file (with `DB_PASSWORD`) to anyone with shell access on the host.

---

## Access

The web service binds to a single host port; **the stack adds no TLS or auth**. Whatever fronts it must.

### Host-local only — default

`WEB_BIND_ADDR=127.0.0.1` (or unset). Reach the UI at `http://localhost:8000` from a browser on the host. Anything not on the host is locked out.

### Remote, single user — SSH tunnel

Same `WEB_BIND_ADDR=127.0.0.1`. From your workstation:

```sh
ssh -L 8000:127.0.0.1:8000 user@<host>
# Then open http://localhost:8000 in your local browser.
```

### Remote, multiple users — front the stack

Set `WEB_BIND_ADDR=0.0.0.0` and put a fronting layer between users and the host port. Common choices:

- **Cloudflare Tunnel (`cloudflared`)** on a separate machine, with a Cloudflare Access policy gating who can sign in. Public hostname, TLS at the edge, no inbound port on the regpoly host.
- **Tailscale** on the regpoly host. Collaborators on the tailnet hit `http://<tailnet-name>:8000`.
- A LAN-only reverse proxy (Caddy, nginx, Traefik) with its own auth (Authentik, Authelia, …).

Make sure the host firewall doesn't also expose `:8000` to anywhere you don't want.

---

## Operations

### Status

```sh
docker compose -f deploy/compose.yml ps
docker compose -f deploy/compose.yml logs -f web         # or worker, db
```

### Restart

```sh
docker compose -f deploy/compose.yml restart web         # safe; worker keeps running
docker compose -f deploy/compose.yml restart worker      # waits up to 75s grace
```

A web restart **never reaps in-flight rows** — the reap lives in the worker startup. A worker restart goes through its SIGTERM handler (cancels in-flight rows cleanly), then the new worker reaps anything orphaned (e.g. by a previous SIGKILL).

### Update to a newer image

```sh
docker compose -f deploy/compose.yml pull
docker compose -f deploy/compose.yml up -d
```

Compose's graceful-stop + the worker's `stop_grace_period: 75s` give the SIGTERM handler enough time to drain. No external orchestration needed.

### Roll back

In `deploy/.env`, set `REGPOLY_TAG=sha-<previous-sha>`, then `docker compose up -d`. The `:sha-<short>` tags on GHCR are immutable.

**Caveat:** rolling back to an image whose schema differs from the current PG state can break. If a deploy includes a non-idempotent migration, document a rollback procedure in the release notes.

### Backups

Out of scope of this stack — the named `pgdata` volume lives at the host path printed by:

```sh
docker volume inspect regpoly_pgdata --format '{{ .Mountpoint }}'
```

Snapshot or replicate that directory externally (NAS snapshot, ZFS / btrfs snapshot, off-site sync). For a live-PG filesystem snapshot to be consistent, the underlying filesystem needs atomic-snapshot semantics. If your filesystem doesn't support that, take logical dumps instead:

```sh
docker compose -f deploy/compose.yml exec -T db \
    pg_dump -U regpoly --format=custom regpoly \
    > regpoly-$(date -u +%F).dump
```

### Migrate from a legacy SQLite database (optional)

If you have rows in a pre-PostgreSQL `regpoly.db`:

```sh
bash deploy/scripts/initial-data-import.sh /path/to/regpoly.db --dry-run
# Review counts; rerun without --dry-run for real.
```

---

## Architecture rationale

- **Migrations in the FastAPI lifespan.** No separate init container. Web's `/healthz` only returns 200 after migrations apply, which is what gates the worker via `depends_on: service_healthy`. If migrations fail, web stays unhealthy and worker stays in `created` — the correct failure mode (refuse to write to an unknown schema).
- **Reap in the worker startup, never the web lifespan.** Single-worker invariant: at worker boot, no other worker exists, so any row in `running` is from a dead predecessor and is safe to cancel (`stale_seconds=0`). The opposite is unsafe — reaping from the web lifespan would clobber rows the live worker is currently processing during a casual `docker compose restart web`.
- **No init / Caddy / cloudflared in the stack.** Init was redundant once the lifespan applies migrations. TLS / auth / public exposure belong out of stack — whichever front layer you pick handles them.
- **Named `pgdata` volume.** Zero host config: no directory to create, no UID-mismatch surprises. Inspect the volume's host path for external backup.
- **Portable image.** Built with `-O3` only (no `-march=native`, no CPU-family tuning). One image runs on x86 and ARM.

---

## Source

- `deploy/Dockerfile` — multi-stage, multi-arch image build.
- `deploy/compose.yml` — three services.
- `deploy/.env.example` — env template.
- `deploy/scripts/initial-data-import.sh` — one-shot SQLite → PG migration tool.
- `.github/workflows/docker.yml` — CI build + multi-arch push to GHCR.

Plan-of-record: `~/.claude/plans/the-architecture-of-the-streamed-lollipop.md`.
