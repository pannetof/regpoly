# regpoly deployment — operator guide

Three-container stack: **db** (PostgreSQL 16), **web** (FastAPI), **worker** (search + analysis schedulers). One `docker compose up -d` brings it up. No init container, no Caddy, no host scripts.

The full design context lives in `~/.claude/plans/the-architecture-of-the-streamed-lollipop.md`. This README covers ops only.

---

## Setup

```sh
cd /path/to/regpoly_monorepo
cp deploy/.env.example deploy/.env
# Edit DB_PASSWORD (required) and WEB_BIND_ADDR (optional) in deploy/.env.
chmod 600 deploy/.env
```

`docker compose config` will otherwise print `deploy/.env` to anyone with shell access — `chmod 600` is what stops that.

## Run

```sh
docker compose -f deploy/compose.yml up -d
```

That's the whole thing. The stack:
- creates the named `pgdata` volume on first run,
- applies any pending DB migrations during the web container's startup (FastAPI lifespan),
- starts the worker once `/healthz` returns 200,
- self-reaps orphan `running` rows when the worker starts (covers crashes / OOM / SIGKILL of a previous worker).

## Update

```sh
docker compose -f deploy/compose.yml pull
docker compose -f deploy/compose.yml up -d
```

Compose's graceful-stop + the worker's `stop_grace_period: 75s` give the worker's SIGTERM handler time to cancel in-flight rows cleanly before the new container takes over. No external orchestration needed.

## Roll back

In `deploy/.env`, set `REGPOLY_TAG=sha-<previous-sha>`, then `docker compose up -d`. The `:sha-<short>` tags on GHCR are immutable.

Caveat: rolling back to an image whose schema differs from the current PG state can break. If a deploy includes a non-idempotent migration, document a rollback procedure in the release notes.

---

## Access

The web service binds to a single host port, controlled by `WEB_BIND_ADDR` in `deploy/.env`.

### Host-local only (default)

`WEB_BIND_ADDR=127.0.0.1` (or unset). Reach the UI at `http://localhost:8000` from a browser on the host.

### Remote, single user

Same `WEB_BIND_ADDR=127.0.0.1`. From your workstation:

```sh
ssh -L 8000:127.0.0.1:8000 user@<host>
# Now open http://localhost:8000 in your local browser.
```

### Remote, multiple users / production

Set `WEB_BIND_ADDR=0.0.0.0` and put something in front. The stack has no TLS or auth — whatever fronts it must:

- **cloudflared** on a separate machine (e.g. a NAS) with a Cloudflare Tunnel pointing at `http://<host>:8000`. Add a Cloudflare Access policy to require Google/email login. No public port on the regpoly host.
- **Tailscale** on the regpoly host. Collaborators on the tailnet reach `http://<tailnet-name>:8000`.
- A LAN-only reverse proxy with its own auth (Authentik, Authelia, Pomerium, …).

Whichever you pick, ensure the host firewall does not also expose :8000 to anywhere you don't want.

---

## Architectures

The published image is multi-arch (`linux/amd64` + `linux/arm64`). `docker compose pull` picks the right one for your host automatically. Runs on x86 mini-PCs, ARM NASes (Synology / QNAP), Raspberry Pi 5s, …

---

## Backups

Out of scope of this stack. The named `pgdata` volume is at the path printed by:

```sh
docker volume inspect regpoly_pgdata --format '{{ .Mountpoint }}'
```

Arrange snapshots externally (NAS filesystem snapshots, ZFS/btrfs, off-site sync, …). For a live-PG snapshot to be consistent, the underlying filesystem needs atomic-snapshot semantics. If your filesystem doesn't support that, take logical dumps instead:

```sh
docker compose -f deploy/compose.yml exec -T db \
    pg_dump -U regpoly --format=custom regpoly > regpoly-$(date -u +%F).dump
```

---

## Routine ops

### Status

```sh
docker compose -f deploy/compose.yml ps
docker compose -f deploy/compose.yml logs -f web      # or worker, db
```

### Restart a single service

```sh
docker compose -f deploy/compose.yml restart web      # safe; worker keeps running
docker compose -f deploy/compose.yml restart worker   # waits up to 75s grace
```

A web restart **never reaps** in-flight rows — reap lives in the worker startup, not the web lifespan, specifically so this is safe.

### Manual reap (after a crashed worker that didn't drain)

Worker startup reaps automatically; you usually don't need to do this manually. If you do:

```sh
docker compose -f deploy/compose.yml stop --timeout 75 worker
docker compose -f deploy/compose.yml start worker
```

The new worker's startup reap finds anything orphaned.

### (Optional) Migrate existing SQLite data

If you have rows in a legacy SQLite `regpoly.db`:

```sh
bash deploy/scripts/initial-data-import.sh /path/to/regpoly.db --dry-run
# Review counts; rerun without --dry-run for real.
```

---

## File index

- `Dockerfile` — multi-stage, multi-arch image build.
- `compose.yml` — three services: db, web, worker.
- `.env.example` — env template (real `.env` is gitignored, `chmod 600`).
- `scripts/initial-data-import.sh` — one-shot SQLite → PG migration tool.
