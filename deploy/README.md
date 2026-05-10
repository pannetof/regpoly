# regpoly deployment — operator guide

Production deployment lives at **https://regpoly.frpanneton.ca** on a
single Hetzner CAX21 (ARM64). Architecture: 5-container Docker
Compose stack (`db`, `init`, `web`, `worker`, `caddy`) backed by
PostgreSQL 16, fronted by Caddy with HTTP basic auth, with daily
NAS-pulled backups via SSH.

The full design context lives in
`~/.claude/plans/i-want-to-dockerize-twinkly-rocket.md`. This README
covers ops only.

---

## First-time setup

### 1. Provision the CAX21

In the [Hetzner Cloud Console](https://console.hetzner.cloud):

1. **Security → SSH Keys**: add your dev workstation's pubkey
   (`~/.ssh/id_ed25519.pub`). Name it `frpan-laptop`.
2. **Servers → Add Server**:
   - Location: **Falkenstein (FSN1)**.
   - Image: **Ubuntu 24.04 LTS (ARM64)**.
   - Type: **CAX21**.
   - Networking: IPv4 (skip IPv6; not in scope today).
   - **No Cloud Firewall** (ufw on the host is the single firewall).
   - SSH keys: select the one from step 1.
3. Note the public IPv4. SSH in once to accept the host key:
   `ssh root@<ip>` → `exit`.

### 2. DNS (Cloudflare)

In the Cloudflare dashboard for `frpanneton.ca`:

1. **DNS → Records**:
   - A record: `regpoly` → `<CAX21 IPv4>`, **proxy: DNS only** (gray
     cloud), TTL Auto.
   - CAA record: `frpanneton.ca` → `0 issue "letsencrypt.org"`.
2. Verify: `dig +short regpoly.frpanneton.ca` returns the CAX21 IP
   (NOT a Cloudflare 104.x.x.x address).

### 3. Bootstrap the host

```sh
ssh root@<ip>
curl -fsSL https://raw.githubusercontent.com/pannetof/regpoly/master/deploy/scripts/bootstrap-cax21.sh -o /tmp/bootstrap.sh
bash /tmp/bootstrap.sh
```

Idempotent. Installs Docker + compose, hardens sshd, configures ufw,
enables NTP, creates the `regpoly` user (in the `docker` group), and
creates the named volumes.

### 4. Clone the repo as `regpoly`

```sh
sudo -iu regpoly bash -c '
    git clone --recursive https://github.com/pannetof/regpoly.git ~/regpoly
'
```

The `--recursive` flag pulls the yaml-cpp submodule.

### 5. Configure secrets

```sh
sudo -iu regpoly bash
cd ~/regpoly
cp deploy/.env.example deploy/.env
vim deploy/.env  # fill in DB_PASSWORD + BASIC_AUTH_HASH
```

Generate the basic-auth hash:

```sh
docker run --rm caddy:2.8-alpine caddy hash-password
# Type your password at the prompt; copy the bcrypt hash into
# BASIC_AUTH_HASH in deploy/.env.
```

Use a strong password (≥20 random chars) — there is no rate limiting.

### 6. Backup-pull script

Copy the wrapper script to a stable path:

```sh
mkdir -p ~/bin
cp ~/regpoly/deploy/scripts/backup-pg-dump.sh ~/bin/
chmod +x ~/bin/backup-pg-dump.sh
```

Add the NAS pubkey to `~/.ssh/authorized_keys` (one line) restricted
to the wrapper:

```
command="/home/regpoly/bin/backup-pg-dump.sh",no-port-forwarding,no-X11-forwarding,no-agent-forwarding,no-pty ssh-ed25519 AAAA... regpoly-backup@nas
```

(Generate the NAS keypair on the NAS with
`ssh-keygen -t ed25519 -f ~/.ssh/regpoly_backup -C regpoly-backup@nas`.)

### 7. First deploy

```sh
bash ~/regpoly/deploy/scripts/deploy.sh
```

Wait ~3-5 min for first-deploy: `initdb` (~25s) + image pull (~30s) +
LE cert acquisition (~30s) + healthchecks (~30s).

**One-time GHCR visibility toggle** after the first GHA push: visit
https://github.com/users/pannetof/packages/container/regpoly/settings
→ "Change package visibility" → Public. Until you do this,
`docker pull` from the CAX21 fails with `denied`.

### 8. UptimeRobot

After the site is up, register a free monitor at
[uptimerobot.com](https://uptimerobot.com):

- Type: **HTTPS**
- URL: `https://regpoly.frpanneton.ca/healthz`
- Keyword: `"db":"up"`
- Interval: 5 min
- Alert: email `frpanneton@gmail.com`, "down for 10 min".

### 9. (Optional) Migrate existing SQLite data

If you have rows in `packages/regpoly-web/var/regpoly.db` you want to
preserve, from your dev machine:

```sh
bash deploy/scripts/initial-data-import.sh \
    /path/to/local/regpoly.db \
    regpoly@regpoly.frpanneton.ca \
    --dry-run
```

Review counts; rerun without `--dry-run` for real.

---

## Routine ops

### Deploy a new version

```sh
# On your dev machine, master branch all committed:
git tag v0.X.Y -m "release notes"
git push --tags
# Wait ~10 min for GHA build (+ smoke test + push to GHCR).
ssh regpoly@regpoly.frpanneton.ca 'bash ~/regpoly/deploy/scripts/deploy.sh'
```

### Roll back

For an already-pushed sha:

```sh
ssh regpoly@regpoly.frpanneton.ca
vim ~/regpoly/deploy/.env   # set REGPOLY_TAG=sha-<previous-sha>
bash ~/regpoly/deploy/scripts/deploy.sh
```

`:sha-<short>` tags are immutable on GHCR.

**Caveat:** rolling back to an image whose schema differs from the
current PG state can break. If a deploy includes a non-idempotent
migration, document a rollback procedure in the release notes.

### Check status

```sh
ssh regpoly@regpoly.frpanneton.ca
cd ~/regpoly
docker compose -f deploy/compose.yml ps
docker compose -f deploy/compose.yml logs -f web      # or worker, db, caddy
```

### Restart a single service

```sh
docker compose -f deploy/compose.yml restart web
# or:
docker compose -f deploy/compose.yml restart worker
```

The `db` service should rarely need a manual restart — restarting it
forces web/worker to reconnect (asyncpg pool handles this within ~10s).

---

## Backups

Backups are pulled by the NAS daily at 03:00 local time. Verify:

```sh
ls -lt /volume1/backups/regpoly/ | head -5
cat /volume1/backups/regpoly/.last-success     # heartbeat
```

If `.last-success` is more than 36h old, the cron stopped working;
investigate the NAS task scheduler.

### Restore SOP

1. From the NAS, identify the dump:
   ```sh
   ls -lt /volume1/backups/regpoly/regpoly-*.dump.gz | head -5
   ```
2. SCP to the CAX21:
   ```sh
   scp /volume1/backups/regpoly/regpoly-2026-MM-DD.dump.gz \
       regpoly@regpoly.frpanneton.ca:/tmp/
   ssh regpoly@regpoly.frpanneton.ca 'gunzip /tmp/regpoly-*.dump.gz'
   ```
3. On the CAX21, drop and recreate:
   ```sh
   ssh regpoly@regpoly.frpanneton.ca
   cd ~/regpoly
   docker compose -f deploy/compose.yml stop web worker
   docker compose -f deploy/compose.yml exec db psql -U regpoly -d postgres \
       -c "DROP DATABASE regpoly; CREATE DATABASE regpoly OWNER regpoly;"
   docker compose -f deploy/compose.yml exec -T db pg_restore -U regpoly -d regpoly --no-owner \
       < /tmp/regpoly-*.dump
   docker compose -f deploy/compose.yml start web worker
   ```
4. Verify with `/healthz` and a UI smoke check.

**Schema-version reconciliation:** the restored DB carries the
schema_version of the dump. If the current image expects a higher
version, the next `init` run applies the missing migrations
automatically. If LOWER (rolling back the image AND restoring an old
dump), restore the old image first.

### Quarterly restore drill

Pick a recent NAS dump, restore it into a throwaway local Docker
stack on your dev machine, and verify row counts match. An untested
backup is a wish.

---

## DR runbook (CAX21 dies entirely)

1. Provision a new CAX21 in Hetzner Console (steps 1.1-1.3 above).
2. Re-add the dev SSH key to the new server (Hetzner UI does this at
   provision time).
3. **Update Cloudflare DNS**: change the A record for `regpoly` to
   the new IP. TTL is Auto (~60s default), so propagation is fast.
4. Bootstrap + clone + configure (steps 1.3-1.6 above).
5. Bring up the stack:
   ```sh
   bash ~/regpoly/deploy/scripts/deploy.sh
   ```
6. SCP the latest dump from NAS to the new CAX21 and run the restore
   SOP above.
7. (Optional) restore `caddy-data` from the most recent
   `caddy-*.tgz.gz` to skip Let's Encrypt re-issuance:
   ```sh
   scp /volume1/backups/regpoly/caddy-2026-MM-DD.tgz \
       regpoly@<new-ip>:/tmp/
   ssh regpoly@<new-ip> "
       docker compose -f ~/regpoly/deploy/compose.yml stop caddy
       docker run --rm -v regpoly_caddy-data:/data alpine sh -c \
           'cd /data && tar xzf -' < /tmp/caddy-*.tgz
       docker compose -f ~/regpoly/deploy/compose.yml start caddy
   "
   ```
8. Verify `/healthz` over HTTPS.

RPO: ~24h (last NAS backup). RTO: ~30-60 min.

---

## Future work (not in scope for first deploy)

- Switch to Cloudflare orange-cloud + DNS-01 ACME via the
  `caddy-dns/cloudflare` plugin (custom xcaddy build) if the origin
  IP is ever scraped/abused.
- Add `worker_heartbeat` table when scaling to multiple worker
  containers (the current `cancel-all-running-on-shutdown` signal
  handler is single-worker-correct only).
- Tighter SSE keepalive cadence if browsers report disconnects.
- Multi-region resilience: managed PG (Hetzner Cloud DB / Supabase)
  + a second CAX21 in HEL1.
- Hourly backups instead of daily if research output volume warrants.
- AAAA / IPv6 plumbing if an IPv6-only client materializes.

---

## File index

- `Dockerfile` — multi-stage image build.
- `compose.yml` — 5-service stack definition.
- `Caddyfile` — TLS + basic auth + reverse proxy.
- `.env.example` — env template (real `.env` is gitignored).
- `scripts/bootstrap-cax21.sh` — one-time host setup.
- `scripts/deploy.sh` — pull + recreate + smoke test.
- `scripts/initial-data-import.sh` — one-shot SQLite → PG migration.
- `scripts/backup-pg-dump.sh` — wrapper invoked by the NAS SSH key.
