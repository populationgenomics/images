#!/usr/bin/env bash
#
# Boot the whole metamist stack inside one container, in order (these are the
# steps of metamist's docs/installation.md):
#   1. start MariaDB (init data dir on first run)
#   2. create the sm_dev database, sm_api user and role
#   3. run the liquibase migrations
#   4. grant the local default user project-creators / members-admin
#   5. start fake-gcs-server
#   6. start the API (uvicorn; --reload only when a checkout is mounted)
#   7. seed once, only if SEED=1 (running SEED_SCRIPT with the image's python)
# then log "[metamist_local_server] ready" and run the given command (defaults to
# tailing the API log) as a child. This script stays PID 1 and supervises: when
# the command or a service exits, or on SIGTERM/SIGINT, it stops what it started
# (command, API, fake-gcs, then MariaDB with a clean shutdown) and exits.
#
# The "[metamist_local_server N/7] ..." phase lines and the final "[metamist_local_server]
# ready" line are an interface: host tooling waits on them. See README.md.

set -euo pipefail

LOG_DIR=/var/log/metamist
DATADIR=/var/lib/mysql
VENV_PY=/opt/venv/bin/python
SEED_MARKER="${DATADIR}/.metamist_seeded"
SEED_IN_PROGRESS="${DATADIR}/.metamist_seed_in_progress"

# PIDs of what this script started, for cleanup(). STEP_PID is the boot step in
# flight (liquibase or the seed script).
MARIADB_PID=
GCS_PID=
API_PID=
CMD_PID=
STEP_PID=
CURRENT_PHASE=

log() { echo "[metamist_local_server] $*"; }
# Numbered boot phases, so a host CLI (or a human) can follow progress;
# "[metamist_local_server] ready" at the end is the line to wait for.
phase() {
    CURRENT_PHASE="$1 $2"
    echo "[metamist_local_server $1] $2"
}
# A fatal error: say why and exit non-zero (cleanup() then stops what is running).
die() {
    log "$*"
    exit 1
}

# An unexpected failure of a top-level command (set -e) names its phase. There is
# no `set -E`, so the trap is not inherited by functions or subshells: a failure
# inside a function is reported with die, and a failing ( ... ) subshell fires
# this trap once, in the shell that ran it.
trap 'log "failed during phase ${CURRENT_PHASE:-startup}"' ERR

# Stop a process this script started: TERM, then KILL if it is still there after 5 s.
stop_process() {
    local name="$1" pid="$2" i
    if [[ -z "${pid}" ]] || ! kill -0 "${pid}" 2>/dev/null; then
        return 0
    fi
    log "stopping ${name}"
    kill -TERM "${pid}" 2>/dev/null || true
    for ((i = 0; i < 50; i++)); do
        kill -0 "${pid}" 2>/dev/null || return 0
        sleep 0.1
    done
    log "${name} did not stop within 5 s, killing it"
    kill -KILL "${pid}" 2>/dev/null || true
}

# Ask MariaDB to shut down cleanly (so the next boot needs no crash recovery) and
# wait for mariadbd to exit. docker stop's timeout is the only bound on the wait.
stop_mariadb() {
    if [[ -z "${MARIADB_PID}" ]] || ! kill -0 "${MARIADB_PID}" 2>/dev/null; then
        return 0
    fi
    log "stopping MariaDB"
    if ! mariadb-admin --protocol=socket shutdown >/dev/null 2>&1; then
        log "mariadb-admin shutdown failed, sending SIGTERM to mariadbd"
        kill -TERM "${MARIADB_PID}" 2>/dev/null || true
    fi
    wait "${MARIADB_PID}" 2>/dev/null || true
}

# Runs on every exit path (normal exit, fatal error, set -e failure, TERM, INT) and
# keeps the exit status it was entered with.
cleanup() {
    local status=$?
    # Finish the shutdown even if another signal arrives, and do not let a failing
    # shutdown command abort it or fire the ERR trap.
    trap '' TERM INT
    trap - ERR
    set +e
    if [[ -n "${STEP_PID}" ]]; then
        # liquibase is a shell script that runs java as a child, so stop children too.
        pkill -TERM -P "${STEP_PID}" 2>/dev/null
        stop_process "boot step" "${STEP_PID}"
    fi
    stop_process "command" "${CMD_PID}"
    stop_process "API" "${API_PID}"
    stop_process "fake-gcs-server" "${GCS_PID}"
    stop_mariadb
    exit "${status}"
}
trap cleanup EXIT
trap 'log "received SIGTERM"; exit 143' TERM
trap 'log "received SIGINT"; exit 130' INT

# wait_for <name> <timeout_s> <logfile> <check command...>
# Poll the check command once a second until it succeeds. On timeout, say so, show
# the end of the service's log and exit 1.
wait_for() {
    local name="$1" timeout="$2" logfile="$3" i
    shift 3
    log "waiting for ${name}"
    for ((i = 0; i < timeout; i++)); do
        if "$@" >/dev/null 2>&1; then
            return 0
        fi
        sleep 1
    done
    log "${name} did not come up"
    tail -n 50 "${logfile}" || true
    exit 1
}

# run_step <dir> <command...>
# Run the command in <dir> as a child and wait for it, returning its status. The
# exec makes the child's PID the command's own, for cleanup(). Waiting on a
# background job lets the TERM/INT traps fire at once; a foreground command
# (liquibase, the seed script) would delay them until it finished.
run_step() {
    local dir="$1" status=0
    shift
    ( cd "${dir}" && exec "$@" ) &
    STEP_PID=$!
    wait "${STEP_PID}" || status=$?
    STEP_PID=
    return "${status}"
}

mkdir -p "${LOG_DIR}"

# Run from the mounted working tree if present, otherwise the baked-in copy.
# Live reload is only useful for a mounted checkout; the baked copy never changes,
# so it is served without the file watcher.
METAMIST_DIR="${METAMIST_DIR:-/app/metamist}"
UVICORN_RELOAD=(--reload)
if [[ ! -f "${METAMIST_DIR}/api/server.py" ]]; then
    log "no source mounted at ${METAMIST_DIR}, using baked /build" >&2
    METAMIST_DIR=/build
    UVICORN_RELOAD=()
fi
cd "${METAMIST_DIR}"
log "metamist source: ${METAMIST_DIR}"

# The default user is interpolated into bootstrap SQL; keep it to a safe charset.
if [[ ! "${SM_LOCALONLY_DEFAULTUSER}" =~ ^[A-Za-z0-9_.@-]+$ ]]; then
    die "SM_LOCALONLY_DEFAULTUSER contains unexpected characters; refusing to start"
fi

# SEED is a switch, not a truthy string: "true" or "yes" must not quietly mean "no".
SEED="${SEED:-0}"
if [[ "${SEED}" != "0" && "${SEED}" != "1" ]]; then
    die "SEED must be 0 or 1, got '${SEED}'; refusing to start"
fi

# There is no default seed: asking for one without naming it is a config error.
# An absolute SEED_SCRIPT runs as-is; a relative one is resolved against the
# metamist source tree.
SEED_PATH=
if [[ "${SEED}" == "1" ]]; then
    if [[ -z "${SEED_SCRIPT:-}" ]]; then
        die "SEED=1 needs SEED_SCRIPT (absolute path, or relative to ${METAMIST_DIR}); refusing to start"
    fi
    SEED_PATH="${SEED_SCRIPT}"
    if [[ "${SEED_PATH}" != /* ]]; then
        SEED_PATH="${METAMIST_DIR}/${SEED_PATH}"
    fi
fi

# The volumes only make sense together with the seed that filled them, so an
# unfinished seed, or a different one than the volume holds, is fatal rather than
# something to run on top of. Checked before anything is started.
RESET_HINT="remove both volumes (${DATADIR} and /data/fakegcs) and start again"
if [[ -f "${SEED_IN_PROGRESS}" ]]; then
    die "seed $(<"${SEED_IN_PROGRESS}") started on this volume and did not finish, so the data may be partial; ${RESET_HINT}; refusing to start"
fi
if [[ "${SEED}" == "1" && -f "${SEED_MARKER}" && "$(<"${SEED_MARKER}")" != "${SEED_PATH}" ]]; then
    die "this volume was seeded with '$(<"${SEED_MARKER}")', not ${SEED_PATH}; ${RESET_HINT}; refusing to start"
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
MARIADB_PID=$!

wait_for MariaDB 60 "${LOG_DIR}/mariadb.log" mariadb-admin --protocol=socket ping

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
LIQUIBASE_ARGS=(
    --classpath=/opt/mariadb-java-client.jar
    --changeLogFile project.xml
    --url "jdbc:mariadb://127.0.0.1:3306/${SM_DEV_DB_NAME}"
    --driver org.mariadb.jdbc.Driver
    --username "${SM_DEV_DB_USER}"
    --password ""
)
# A stop during an earlier update leaves liquibase's changelog lock held, and the
# next update then waits on it. This container is the only migrator, so any held
# lock is stale.
run_step "${METAMIST_DIR}/db" liquibase "${LIQUIBASE_ARGS[@]}" releaseLocks
run_step "${METAMIST_DIR}/db" liquibase "${LIQUIBASE_ARGS[@]}" update

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
GCS_PID=$!

wait_for fake-gcs-server 30 "${LOG_DIR}/fakegcs.log" \
    curl -fs "http://localhost:4443/storage/v1/b?project=${GOOGLE_CLOUD_PROJECT}"

# --- 6. API server ----------------------------------------------------------
phase 6/7 "starting API on :8000 (uvicorn ${UVICORN_RELOAD[*]:-without reload})"
( cd "${METAMIST_DIR}" && exec "${VENV_PY}" -m uvicorn \
    --host 0.0.0.0 --port 8000 "${UVICORN_RELOAD[@]}" api.server:app ) \
    >"${LOG_DIR}/api.log" 2>&1 &
API_PID=$!

wait_for API 60 "${LOG_DIR}/api.log" curl -fs http://localhost:8000/openapi.json

# --- 7. seed once, if asked -------------------------------------------------
if [[ "${SEED}" == "0" ]]; then
    phase 7/7 "seed skipped (SEED=0)"
elif [[ -f "${SEED_MARKER}" ]]; then
    # Its content was checked against SEED_PATH before anything started.
    phase 7/7 "seed skipped (already seeded with ${SEED_PATH})"
else
    # A missing or failing seed is fatal: the container never reports ready with
    # less data than was asked for.
    if [[ ! -f "${SEED_PATH}" ]]; then
        die "seed script ${SEED_PATH} not found; refusing to start"
    fi
    phase 7/7 "seeding: ${SEED_PATH}"
    # The in-progress marker stays if the seed fails or is interrupted: the data is
    # then partial, and rerunning the script on top of it is not safe.
    echo "${SEED_PATH}" >"${SEED_IN_PROGRESS}"
    seed_status=0
    run_step "${METAMIST_DIR}" "${VENV_PY}" "${SEED_PATH}" || seed_status=$?
    if ((seed_status != 0)); then
        die "seed script failed (exit ${seed_status}), so the data may be partial; ${RESET_HINT}; refusing to start"
    fi
    # Success marker first: a stop between the two writes leaves both files, which
    # the next boot refuses, rather than neither, which it would seed over.
    echo "${SEED_PATH}" >"${SEED_MARKER}"
    rm -f "${SEED_IN_PROGRESS}"
    log "seed complete"
fi

log "ready"

API_HOST_PORT="${HOST_API_PORT:-8000}"
cat <<BANNER

  metamist_local_server is up:
    Swagger    http://localhost:${API_HOST_PORT}/docs
    GraphiQL   http://localhost:${API_HOST_PORT}/graphql
    fake GCS   http://localhost:${GCS_HOST_PORT}  (STORAGE_EMULATOR_HOST)
    user       ${SM_LOCALONLY_DEFAULTUSER}  (no token needed in local mode)

  Query it:    mmquery '{ myProjects { id name } }'
  Logs:        ${LOG_DIR}/api.log  ${LOG_DIR}/mariadb.log

BANNER

# Run the command as a child (not exec) so this script can stop the services when
# it ends. <&0 keeps an interactive command's stdin; a background job would
# otherwise get /dev/null.
"$@" <&0 &
CMD_PID=$!

# Block until the command or a service exits. wait -n returns that job's status;
# || keeps a non-zero status from tripping set -e and the ERR trap.
exited_pid=
status=0
wait -n -p exited_pid "${CMD_PID}" "${MARIADB_PID}" "${GCS_PID}" "${API_PID}" || status=$?

case "${exited_pid}" in
    "${CMD_PID}")
        log "command exited with status ${status}"
        exit "${status}"
        ;;
    "${MARIADB_PID}")
        service=MariaDB logfile="${LOG_DIR}/mariadb.log"
        ;;
    "${GCS_PID}")
        service=fake-gcs-server logfile="${LOG_DIR}/fakegcs.log"
        ;;
    "${API_PID}")
        service=API logfile="${LOG_DIR}/api.log"
        ;;
    *)
        die "wait returned without a known process (status ${status})"
        ;;
esac
log "${service} exited unexpectedly (status ${status})"
tail -n 50 "${logfile}" || true
exit 1
