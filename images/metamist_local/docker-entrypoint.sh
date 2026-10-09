#!/usr/bin/env bash
#
# Boot the whole metamist stack inside one container, in order (these are the
# steps of metamist's docs/installation.md):
#   1. start MariaDB (init data dir on first run)
#   2. create the sm_dev database, sm_api user and role
#   3. run the liquibase migrations
#   4. grant the local default user project-creators / members-admin
#   5. start fake-gcs-server
#   6. start the API (uvicorn --reload)
#   7. seed once, only if SEED=1 (running SEED_SCRIPT with the image's python)
# then log "[metamist_local] ready" and hand control to the given command
# (defaults to tailing the API log).
#
# The "[metamist_local N/7] ..." phase lines and the final "[metamist_local]
# ready" line are an interface: host tooling waits on them. See README.md.

set -euo pipefail

LOG_DIR=/var/log/metamist
DATADIR=/var/lib/mysql
VENV_PY=/opt/venv/bin/python
SEED_MARKER="${DATADIR}/.metamist_seeded"

log() { echo "[metamist_local] $*"; }
# Numbered boot phases, so a host CLI (or a human) can follow progress;
# "[metamist_local] ready" at the end is the line to wait for.
phase() { echo "[metamist_local $1] $2"; }

mkdir -p "${LOG_DIR}"

# Run from the mounted working tree if present, otherwise the baked-in copy.
METAMIST_DIR="${METAMIST_DIR:-/app/metamist}"
if [[ ! -f "${METAMIST_DIR}/api/server.py" ]]; then
    log "no source mounted at ${METAMIST_DIR}, using baked /build" >&2
    METAMIST_DIR=/build
fi
cd "${METAMIST_DIR}"
log "metamist source: ${METAMIST_DIR}"

# The default user is interpolated into bootstrap SQL; keep it to a safe charset.
if [[ ! "${SM_LOCALONLY_DEFAULTUSER}" =~ ^[A-Za-z0-9_.@-]+$ ]]; then
    log "SM_LOCALONLY_DEFAULTUSER contains unexpected characters; refusing to start"
    exit 1
fi

# There is no default seed: asking for one without naming it is a config error.
if [[ "${SEED:-0}" == "1" && -z "${SEED_SCRIPT:-}" ]]; then
    log "SEED=1 needs SEED_SCRIPT (absolute path, or relative to ${METAMIST_DIR}); refusing to start"
    exit 1
fi

# --- 1. MariaDB -------------------------------------------------------------
if [[ ! -d "${DATADIR}/mysql" ]]; then
    log "initialising MariaDB data directory"
    mariadb-install-db --user=mysql --datadir="${DATADIR}" \
        --auth-root-authentication-method=normal >/dev/null
fi
chown -R mysql:mysql "${DATADIR}"

# mariadbd needs its runtime socket dir; /run is tmpfs so create it each boot.
mkdir -p /run/mysqld && chown mysql:mysql /run/mysqld

phase 1/7 "starting MariaDB"
gosu mysql mariadbd --datadir="${DATADIR}" --bind-address=0.0.0.0 --port=3306 \
    >"${LOG_DIR}/mariadb.log" 2>&1 &

log "waiting for MariaDB to accept connections"
for _ in $(seq 1 60); do
    if mariadb-admin --protocol=socket ping >/dev/null 2>&1; then break; fi
    sleep 1
done
mariadb-admin --protocol=socket ping >/dev/null 2>&1 \
    || { log "MariaDB did not come up"; cat "${LOG_DIR}/mariadb.log"; exit 1; }

# --- 2. database / user / role (idempotent) ---------------------------------
phase 2/7 "ensuring database, user and role exist"
mariadb -u root <<SQL
CREATE DATABASE IF NOT EXISTS ${SM_DEV_DB_NAME};
CREATE USER IF NOT EXISTS '${SM_DEV_DB_USER}'@'%';
CREATE USER IF NOT EXISTS '${SM_DEV_DB_USER}'@'localhost';
CREATE ROLE IF NOT EXISTS sm_api_role;
GRANT sm_api_role TO '${SM_DEV_DB_USER}'@'%';
GRANT sm_api_role TO '${SM_DEV_DB_USER}'@'localhost';
SET DEFAULT ROLE sm_api_role FOR '${SM_DEV_DB_USER}'@'%';
SET DEFAULT ROLE sm_api_role FOR '${SM_DEV_DB_USER}'@'localhost';
GRANT ALL PRIVILEGES ON ${SM_DEV_DB_NAME}.* TO sm_api_role;
SQL

# --- 3. migrations ----------------------------------------------------------
phase 3/7 "running liquibase migrations"
( cd "${METAMIST_DIR}/db" && liquibase \
    --classpath=/opt/mariadb-java-client.jar \
    --changeLogFile project.xml \
    --url "jdbc:mariadb://127.0.0.1:3306/${SM_DEV_DB_NAME}" \
    --driver org.mariadb.jdbc.Driver \
    --username "${SM_DEV_DB_USER}" \
    --password "" \
    update )

# --- 4. grant the local default user project-creator rights (idempotent) ----
phase 4/7 "granting ${SM_LOCALONLY_DEFAULTUSER} project-creators / members-admin"
mariadb -h 127.0.0.1 -P 3306 -u "${SM_DEV_DB_USER}" "${SM_DEV_DB_NAME}" <<SQL
INSERT INTO group_member (group_id, member)
SELECT g.id, '${SM_LOCALONLY_DEFAULTUSER}'
FROM \`group\` g
WHERE g.name IN ('project-creators', 'members-admin')
  AND NOT EXISTS (
      SELECT 1 FROM group_member gm
      WHERE gm.group_id = g.id AND gm.member = '${SM_LOCALONLY_DEFAULTUSER}'
  );
SQL

# --- 5. fake GCS emulator ---------------------------------------------------
# Filesystem backend on a volume so uploaded objects survive container recreation
# and stay in sync with the persisted DB. (With a memory backend, a recreate
# wipes the objects while the DB + seed marker persist and the seed, which
# uploads them, never re-runs, so analysis outputs would look missing.)
FAKEGCS_ROOT=/data/fakegcs
mkdir -p "${FAKEGCS_ROOT}"
# Binds all container interfaces so the port can be published to the host.
# The URLs it hands out (e.g. resumable-upload sessions) use the *host-published*
# port, which HOST_GCS_PORT carries in; it differs from 4443 when the container
# port is published on another host port.
GCS_HOST_PORT="${HOST_GCS_PORT:-4443}"
phase 5/7 "starting fake-gcs-server on :4443 (filesystem backend at ${FAKEGCS_ROOT})"
fake-gcs-server -scheme http -host 0.0.0.0 -port 4443 \
    -backend filesystem -filesystem-root "${FAKEGCS_ROOT}" \
    -external-url "http://localhost:${GCS_HOST_PORT}" -public-host "localhost:${GCS_HOST_PORT}" \
    >"${LOG_DIR}/fakegcs.log" 2>&1 &
for _ in $(seq 1 30); do
    if curl -fs "http://localhost:4443/storage/v1/b?project=${GOOGLE_CLOUD_PROJECT}" >/dev/null 2>&1; then
        break
    fi
    sleep 1
done

# --- 6. API server ----------------------------------------------------------
phase 6/7 "starting API on :8000 (uvicorn --reload)"
( cd "${METAMIST_DIR}" && exec "${VENV_PY}" -m uvicorn \
    --host 0.0.0.0 --port 8000 --reload api.server:app ) \
    >"${LOG_DIR}/api.log" 2>&1 &

log "waiting for the API to respond"
for _ in $(seq 1 60); do
    if curl -fs http://localhost:8000/openapi.json >/dev/null 2>&1; then break; fi
    sleep 1
done
curl -fs http://localhost:8000/openapi.json >/dev/null 2>&1 \
    || { log "API did not come up"; tail -n 50 "${LOG_DIR}/api.log"; exit 1; }

# --- 7. seed once, if asked -------------------------------------------------
if [[ "${SEED:-0}" == "1" && ! -f "${SEED_MARKER}" ]]; then
    # An absolute SEED_SCRIPT runs as-is; a relative one is resolved against the
    # metamist source tree.
    seed_path="${SEED_SCRIPT}"
    [[ "${seed_path}" != /* ]] && seed_path="${METAMIST_DIR}/${seed_path}"

    # A missing or failing seed is fatal: the container never reports ready with
    # less data than was asked for. No marker is written, so the next boot retries.
    if [[ ! -f "${seed_path}" ]]; then
        log "seed script ${seed_path} not found; refusing to start"
        exit 1
    fi
    phase 7/7 "seeding: ${seed_path}"
    if ! ( cd "${METAMIST_DIR}" && "${VENV_PY}" "${seed_path}" ); then
        log "seed script failed; refusing to start"
        exit 1
    fi
    touch "${SEED_MARKER}"
    log "seed complete"
else
    phase 7/7 "seed skipped (already seeded, or SEED=0)"
fi

log "ready"

API_HOST_PORT="${HOST_API_PORT:-8000}"
cat <<BANNER

  metamist_local is up:
    Swagger    http://localhost:${API_HOST_PORT}/docs
    GraphiQL   http://localhost:${API_HOST_PORT}/graphql
    fake GCS   http://localhost:${GCS_HOST_PORT}  (STORAGE_EMULATOR_HOST)
    user       ${SM_LOCALONLY_DEFAULTUSER}  (no token needed in local mode)

  Query it:    mmquery '{ myProjects { id name } }'
  Logs:        ${LOG_DIR}/api.log  ${LOG_DIR}/mariadb.log

BANNER

exec "$@"
