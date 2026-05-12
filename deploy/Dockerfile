# syntax=docker/dockerfile:1.7
#
# regpoly_monorepo Docker image.
#
# Single image, two roles. The same `ghcr.io/pannetof/regpoly:<tag>`
# runs as the `web` or `worker` service per the compose command (the
# official `postgres:16-alpine` covers the third service).
#
# Built portably (no -march/-mcpu) so the same image can be published
# as linux/amd64 + linux/arm64 and run on whatever Docker host the
# operator points it at.

# ─── Builder stage ───────────────────────────────────────────────────
FROM debian:12-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake git ca-certificates curl \
        libntl-dev libgmp-dev libgf2x-dev libpq-dev \
        python3 python3-dev python3-venv \
    && rm -rf /var/lib/apt/lists/*

# Install uv via the official installer — faster + more current than apt.
RUN curl -LsSf https://astral.sh/uv/install.sh | sh \
    && mv /root/.local/bin/uv /usr/local/bin/uv

# Force uv to use system Python 3.11; never download a managed
# interpreter. Without this, uv may pick a python-build-standalone
# build whose pybind11 ABI doesn't match the runtime stage's libpython.
# UV_PREFER_BINARY=1 favours wheels over sdist (matters for psycopg
# binary + asyncpg, both of which ship aarch64 wheels).
ENV UV_PYTHON_DOWNLOADS=never \
    UV_PYTHON_PREFERENCE=only-system \
    UV_PREFER_BINARY=1

WORKDIR /app

# Layer cache strategy: COPY in increasing-change-frequency order so
# pure-Python edits don't invalidate the C++ extension's build cache.
#
# 1. Lockfile + workspace pyproject (stable; changes only on dep bumps).
COPY pyproject.toml uv.lock ./

# 2. Per-package pyproject.toml (stable).
COPY packages/regpoly-cpp/pyproject.toml      packages/regpoly-cpp/pyproject.toml
COPY packages/regpoly/pyproject.toml          packages/regpoly/pyproject.toml
COPY packages/regpoly-web/pyproject.toml      packages/regpoly-web/pyproject.toml

# 3. C++ source + vendored yaml-cpp (changes occasionally; rebuild
#    expensive). yaml-cpp is a git submodule under third_party/, so
#    GHA must check out with `submodules: recursive`.
COPY third_party/yaml-cpp third_party/yaml-cpp
COPY packages/regpoly-cpp packages/regpoly-cpp

# 4. Pure-Python source (changes often; rebuild cheap).
COPY packages/regpoly      packages/regpoly
COPY packages/regpoly-web  packages/regpoly-web

# 5. Catalog/data (changes rarely).
COPY shared shared
COPY docs/library docs/library
COPY docs/papers docs/papers

# 6. Build the workspace. Cache uv's download cache and scikit-build's
#    intermediate build dir between builds for fast iteration on
#    BuildKit-enabled hosts.
#    `--no-editable` installs workspace members fully into .venv so the
#    runtime image only needs /app/.venv (no source dirs).
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=cache,target=/app/packages/regpoly-cpp/_skbuild \
    CMAKE_ARGS="-DREGPOLY_NATIVE_ARCH=OFF" \
    uv sync --no-dev --frozen --no-editable --python /usr/bin/python3.11

# ─── Runtime stage ──────────────────────────────────────────────────
FROM debian:12-slim AS runtime

# Runtime libs only — no -dev packages, no compilers, no git/uv.
# - libntl44 / libgmp10 / libgf2x3: regpoly_cpp's transitive deps.
# - libgomp1: defensive in case any C++ source uses #pragma omp.
# - libpq5: psycopg's runtime client lib.
# - tini: PID 1 + signal forwarding so worker SIGTERMs reach the
#   ProcessPoolExecutor children (which inherit the parent's pgrp).
# - coreutils: stat / date for the worker container's healthcheck
#   (touchfile age check — see compose.yml).
RUN apt-get update && apt-get install -y --no-install-recommends \
        libntl44 libgmp10 libgf2x3 libgomp1 libpq5 \
        python3 tini ca-certificates coreutils \
    && rm -rf /var/lib/apt/lists/* \
    && groupadd -g 10001 app \
    && useradd -u 10001 -g 10001 -M -d /app -s /usr/sbin/nologin app

# PYTHONDONTWRITEBYTECODE=1 keeps the venv read-only at runtime (no
# __pycache__ writes from uid 10001 into root-owned dirs).
# PYTHONUNBUFFERED=1 makes stdout/stderr line-buffered so `docker logs`
# is timely.
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PATH=/app/.venv/bin:$PATH \
    REGPOLY_LIBRARY_DIR=/app/share/regpoly/library \
    REGPOLY_PAPERS_DIR=/app/share/regpoly/papers

WORKDIR /app

# Multi-stage COPY with explicit ownership: without --chown, files
# land as root:root and the runtime user (10001) can't write any
# adjacent state. The .venv contents are read by uid 10001; package
# data + catalog files likewise.
# With `--no-editable` the workspace members (regpoly_cpp, regpoly,
# regpoly_web) are fully installed under /app/.venv/lib/.../site-packages,
# so we no longer COPY the source tree at runtime.
COPY --from=builder --chown=10001:10001 /app/.venv                       /app/.venv
COPY --from=builder --chown=10001:10001 /app/packages/regpoly-web/scripts /app/packages/regpoly-web/scripts
COPY --from=builder --chown=10001:10001 /app/docs/library                /app/share/regpoly/library
COPY --from=builder --chown=10001:10001 /app/docs/papers                 /app/share/regpoly/papers
COPY --from=builder --chown=10001:10001 /app/shared                      /app/share/regpoly/shared

# Strip dev-only shell launchers (compose passes its own command).
RUN rm -f /app/packages/regpoly-web/scripts/wipe-db.sh \
       && rm -f /app/packages/regpoly-web/scripts/regpoly-web.sh

USER 10001:10001

# tini -s: subreaper (harmless when PID 1; future-proof if invoked
# under another init).
# tini -g: forward signals to the entire process group, so SIGTERM
# reaches the ProcessPoolExecutor children spawned by the worker.
ENTRYPOINT ["/usr/bin/tini", "-s", "-g", "--"]
EXPOSE 8000
# No CMD — compose's `command:` selects the role per service.
