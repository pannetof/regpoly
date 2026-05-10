#!/usr/bin/env bash
# bootstrap-cax21.sh — one-time host setup for the regpoly Hetzner CAX21
# (dockerize-plan Phase 7).
#
# Idempotent: safe to re-run. Run as root (sudo) on a fresh Ubuntu
# 24.04 LTS ARM64 server provisioned via the Hetzner Cloud Console.
#
# Side effects:
#   - apt: ca-certificates, curl, ufw + Docker Engine + compose plugin
#   - useradd regpoly (in `docker` group)
#   - sshd hardening drop-in at /etc/ssh/sshd_config.d/99-regpoly.conf
#   - ufw: deny incoming, limit 22/tcp, allow 80/443/tcp, enable
#   - timedatectl: NTP on
#   - docker volumes created: pgdata, caddy-data, caddy-config

set -euo pipefail

if [ "$(id -u)" -ne 0 ]; then
    echo "error: run as root (sudo bash $0)" >&2
    exit 1
fi

if [ "$(uname -m)" != "aarch64" ]; then
    echo "warning: this script targets ARM64 (uname -m=$(uname -m)); proceeding anyway"
fi

echo "== apt update + base packages"
apt-get update
DEBIAN_FRONTEND=noninteractive apt-get install -y \
    ca-certificates curl ufw

echo "== Docker Engine + compose plugin"
if ! command -v docker >/dev/null; then
    curl -fsSL https://get.docker.com | sh
fi
if ! dpkg -l | grep -q docker-compose-plugin; then
    apt-get install -y docker-compose-plugin
fi

echo "== regpoly user"
if ! id -u regpoly >/dev/null 2>&1; then
    useradd -m -s /bin/bash regpoly
fi
usermod -aG docker regpoly

echo "== SSH hardening (sshd_config.d/99-regpoly.conf)"
cat >/etc/ssh/sshd_config.d/99-regpoly.conf <<'EOF'
# Managed by deploy/scripts/bootstrap-cax21.sh
PermitRootLogin no
PasswordAuthentication no
KbdInteractiveAuthentication no
MaxAuthTries 3
LoginGraceTime 20
EOF
systemctl reload ssh || systemctl reload sshd || true

echo "== ufw firewall"
ufw --force default deny incoming
ufw --force default allow outgoing
ufw --force limit  22/tcp
ufw --force allow  80/tcp
ufw --force allow  443/tcp
ufw --force enable

echo "== NTP"
timedatectl set-ntp true
systemctl enable --now systemd-timesyncd

echo "== Docker volumes"
docker volume create pgdata >/dev/null
docker volume create caddy-data >/dev/null
docker volume create caddy-config >/dev/null

echo
echo "✓ bootstrap complete."
echo
echo "Next steps:"
echo "  1. As the regpoly user, clone the repo:"
echo "       sudo -iu regpoly bash -c 'git clone --recursive https://github.com/pannetof/regpoly.git ~/regpoly'"
echo "  2. Copy deploy/.env.example → deploy/.env and fill in values."
echo "  3. Run deploy/scripts/deploy.sh."
echo "  4. After first GHCR push, toggle the package visibility to public:"
echo "       https://github.com/users/pannetof/packages/container/regpoly/settings"
