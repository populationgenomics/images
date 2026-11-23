#!/usr/bin/env bash
#
# Boot the whole metamist stack inside one container, in order (these are the
# steps of metamist's docs/installation.md):
#   1. start PostgreSQL (initdb on first run)
#   2. run the dbmate migrations (this also creates the metamist database)
#   3. let the metamist_server role log in, and grant the local default user
#      project-creators / members-admin
#   4. start fake-gcs-server
#   5. start the API (uvicorn; --reload only when a checkout is mounted)
#   6. seed once, only if SEED=1 (running SEED_SCRIPT with the image's python)
# then log "[metamist_local_server] ready" and run the given command (defaults to
# tailing the API log) as a child. This script stays PID 1 and supervises: when
# the command or a service exits, or on SIGTERM/SIGINT, it stops what it started
# (command, API, fake-gcs, then PostgreSQL with a clean shutdown) and exits.
#
# The "[metamist_local_server N/6] ..." phase lines and the final "[metamist_local_server]
# ready" line are an interface: host tooling waits on them. See README.md.

set -euo pipefail

LOG_DIR=/var/log/metamist
VENV_PY=/opt/venv/bin/python

# PGDATA and PG_MAJOR come from the image. The seed markers live inside the data
# directory, so they exist exactly as long as the data they describe: a volume
# that has lost its data has lost its markers too.
: "${PGDATA:?PGDATA is not set (the image sets it)}"
: "${PG_MAJOR:?PG_MAJOR is not set (the image sets it)}"
PG_BIN="/usr/lib/postgresql/${PG_MAJOR}/bin"
PG_VOLUME=/var/lib/postgresql
PG_SOCKET_DIR=/var/run/postgresql
PG_HBA_FILE=/etc/metamist/pg_hba.conf
SEED_MARKER="${PGDATA}/.metamist_seeded"
SEED_IN_PROGRESS="${PGDATA}/.metamist_seed_in_progress"

# Names metamist's migrations hard-code: the bootstrap superuser that owns the
# schemas (ALTER DEFAULT PRIVILEGES FOR ROLE ...), the database they grant CONNECT
# on, and the application login role they create (as NOLOGIN).
PG_SUPERUSER=metamist_superuser
PG_DATABASE=metamist
PG_APP_ROLE=metamist_server

# PIDs of what this script started, for cleanup(). STEP_PID is the boot step in
# flight (dbmate or the seed script); SLEEP_PID is the sleep of the supervision
# loop at the end.
POSTGRES_PID=
GCS_PID=
API_PID=
CMD_PID=
STEP_PID=
SLEEP_PID=
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

# wait_gone <pid> <tenths of a second>: return 0 once the process is gone, 1 if it
# is still there after that long.
wait_gone() {
    local pid="$1" i
    for ((i = 0; i < $2; i++)); do
        kill -0 "${pid}" 2>/dev/null || return 0
        sleep 0.1
    done
    return 1
}

# Ask PostgreSQL to shut down cleanly (so the next boot needs no crash recovery) and
# wait for it to exit. SIGINT is its "fast" shutdown: it disconnects clients, rolls
# back open transactions and checkpoints. (SIGTERM, "smart", would wait for every
# client to disconnect first.) If that takes more than 20 s, SIGQUIT is its
# "immediate" shutdown: no checkpoint, so the next boot runs crash recovery, but it
# is quick. SIGKILL is the last resort, as for the other processes.
stop_postgres() {
    if [[ -z "${POSTGRES_PID}" ]] || ! kill -0 "${POSTGRES_PID}" 2>/dev/null; then
        return 0
    fi
    log "stopping PostgreSQL"
    kill -INT "${POSTGRES_PID}" 2>/dev/null || true
    wait_gone "${POSTGRES_PID}" 200 && return 0
    log "PostgreSQL did not shut down within 20 s, forcing an immediate shutdown"
    kill -QUIT "${POSTGRES_PID}" 2>/dev/null || true
    wait_gone "${POSTGRES_PID}" 50 && return 0
    log "PostgreSQL did not stop within 5 s of the immediate shutdown, killing it"
    kill -KILL "${POSTGRES_PID}" 2>/dev/null || true
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
        # A seed script may have started children of its own, so stop them too.
        pkill -TERM -P "${STEP_PID}" 2>/dev/null
        stop_process "boot step" "${STEP_PID}"
    fi
    # The supervision loop's sleep may be in flight when a signal ends the loop.
    if [[ -n "${SLEEP_PID}" ]]; then
        kill "${SLEEP_PID}" 2>/dev/null
    fi
    stop_process "command" "${CMD_PID}"
    stop_process "API" "${API_PID}"
    stop_process "fake-gcs-server" "${GCS_PID}"
    stop_postgres
    exit "${status}"
}
trap cleanup EXIT
trap 'log "received SIGTERM"; exit 143' TERM
trap 'log "received SIGINT"; exit 130' INT

# wait_for <name> <pid> <timeout_s> <logfile> <check command...>
# Poll the check command once a second until it succeeds. If the service's process
# <pid> exits first, or on timeout, say so, show the end of the service's log and
# exit 1.
wait_for() {
    local name="$1" pid="$2" timeout="$3" logfile="$4" i status
    shift 4
    log "waiting for ${name}"
    for ((i = 0; i < timeout; i++)); do
        if ! kill -0 "${pid}" 2>/dev/null; then
            status=0
            wait "${pid}" 2>/dev/null || status=$?
            log "${name} exited before it came up (status ${status})"
            tail -n 50 "${logfile}" || true
            exit 1
        fi
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
# (dbmate, the seed script) would delay them until it finished.
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

# metamist 8 moved from MariaDB/liquibase to PostgreSQL/dbmate. A mounted 7.x
# checkout has no db/migrations, and this image has no MariaDB to run it on.
if [[ ! -d "${METAMIST_DIR}/db/migrations" ]]; then
    die "${METAMIST_DIR}/db/migrations not found: this image runs metamist 8 or later (PostgreSQL); a 7.x checkout needs an older image; refusing to start"
fi

# The image never creates /var/lib/mysql, so it exists only if a 7.x setup's
# MariaDB volume is still mounted there. This image cannot read it, and the
# PostgreSQL data goes elsewhere, so a stale mount would be silently ignored.
if [[ -e /var/lib/mysql ]]; then
    die "/var/lib/mysql exists, so an old MariaDB (7.x) volume is still mounted; this image keeps its data in PostgreSQL at ${PG_VOLUME} and cannot use it (see 'Upgrading from 7.14.x' in the README); refusing to start"
fi

# The migrations create the metamist_server role and grant CONNECT on the metamist
# database by name, so other names cannot work. Fail here, rather than later with
# a connection error, if a 7.x setup still passes its old SM_DEV_DB_* overrides.
if [[ "${SM_DEV_DB_NAME}" != "${PG_DATABASE}" || "${SM_DEV_DB_USER}" != "${PG_APP_ROLE}" ]]; then
    die "SM_DEV_DB_NAME must be ${PG_DATABASE} and SM_DEV_DB_USER ${PG_APP_ROLE} (metamist's migrations fix those names), got '${SM_DEV_DB_NAME}' and '${SM_DEV_DB_USER}'; refusing to start"
fi

# The default user ends up in SQL and in log lines; keep it to a safe charset even
# though psql quotes it as a literal.
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
RESET_HINT="remove both volumes (${PG_VOLUME} and /data/fakegcs) and start again"
if [[ -f "${SEED_IN_PROGRESS}" ]]; then
    die "seed $(<"${SEED_IN_PROGRESS}") started on this volume and did not finish, so the data may be partial; ${RESET_HINT}; refusing to start"
fi
if [[ "${SEED}" == "1" && -f "${SEED_MARKER}" && "$(<"${SEED_MARKER}")" != "${SEED_PATH}" ]]; then
    die "this volume was seeded with '$(<"${SEED_MARKER}")', not ${SEED_PATH}; ${RESET_HINT}; refusing to start"
fi
# A missing script is only an error when the seed is going to run: a volume that is
# already seeded with this path skips it.
if [[ "${SEED}" == "1" && ! -f "${SEED_MARKER}" && ! -f "${SEED_PATH}" ]]; then
    die "seed script ${SEED_PATH} not found; refusing to start"
fi

# --- 1. PostgreSQL ----------------------------------------------------------
# PGDATA is only ever created whole: initdb runs in a scratch directory beside it,
# which is renamed to PGDATA (atomic, on the same volume) once initdb has finished.
# A stop in the middle of initdb therefore leaves PGDATA absent, never half-written,
# and the next boot starts over. The postgres user owns the volume root and the
# directory above PGDATA (so it can make the scratch directory there) and everything
# it creates; nothing is chowned recursively.
PGDATA_PARENT="$(dirname "${PGDATA}")"
PGDATA_SCRATCH="${PGDATA}.initdb"
mkdir -p "${PGDATA_PARENT}"
chown postgres:postgres "${PG_VOLUME}" "${PGDATA_PARENT}"

if [[ -e "${PGDATA}" && ! -d "${PGDATA}" ]]; then
    die "${PGDATA} exists but is not a directory; ${RESET_HINT}; refusing to start"
elif [[ -d "${PGDATA}" && -n "$(ls -A "${PGDATA}")" ]]; then
    # A cluster has both of these files. Without them this was not made by the
    # initdb above (which only renames a finished directory into place), so do not
    # start PostgreSQL on it, and do not delete what is there either.
    if [[ ! -f "${PGDATA}/PG_VERSION" || ! -f "${PGDATA}/global/pg_control" ]]; then
        die "${PGDATA} is not a complete PostgreSQL data directory (PG_VERSION or global/pg_control is missing); ${RESET_HINT}; refusing to start"
    fi
    chown postgres:postgres "${PGDATA}"
else
    log "initialising PostgreSQL data directory"
    if [[ -e "${PGDATA_SCRATCH}" ]]; then
        log "removing ${PGDATA_SCRATCH}, left by an interrupted first boot"
        rm -rf "${PGDATA_SCRATCH}"
    fi
    # The superuser is named for what metamist's migrations expect (see above). The
    # locale is the one metamist's own db image (postgres:18) uses.
    if ! gosu postgres "${PG_BIN}/initdb" -D "${PGDATA_SCRATCH}" -U "${PG_SUPERUSER}" \
        --encoding=UTF8 --locale=en_US.utf8 --auth-local=trust --auth-host=reject \
        >"${LOG_DIR}/initdb.log" 2>&1; then
        tail -n 50 "${LOG_DIR}/initdb.log" || true
        rm -rf "${PGDATA_SCRATCH}"
        die "initdb failed; refusing to start"
    fi
    # An empty PGDATA (a mount point cannot be replaced) is removed so the rename can
    # take its place; rmdir fails, and so refuses, if it is a mount point.
    if [[ -d "${PGDATA}" ]]; then
        rmdir "${PGDATA}" || die "${PGDATA} is an empty directory that cannot be replaced (a volume mounted on it?); mount the volume at ${PG_VOLUME} instead; refusing to start"
    fi
    mv "${PGDATA_SCRATCH}" "${PGDATA}"
fi

# Who may connect, written on every boot so it always matches this image rather than
# whatever the volume was initialised with. There is no password anywhere, so:
#   - the superuser (and any other role) can connect only over the local socket,
#     which means from inside the container;
#   - the application login is the only role allowed over TCP, which is what a
#     published port or the API itself uses. It is not a superuser and only sees
#     the metamist database.
# Publish the port on loopback only; see README.md.
cat >"${PG_HBA_FILE}" <<HBA
local  all             all                                   trust
host   ${PG_DATABASE}  ${PG_APP_ROLE}  0.0.0.0/0             trust
host   ${PG_DATABASE}  ${PG_APP_ROLE}  ::/0                  trust
HBA

# The server needs its runtime socket dir; /run does not survive a container
# recreate, so create it each boot.
mkdir -p "${PG_SOCKET_DIR}" && chown postgres:postgres "${PG_SOCKET_DIR}"

phase 1/6 "starting PostgreSQL"
gosu postgres "${PG_BIN}/postgres" -D "${PGDATA}" \
    -c listen_addresses='*' \
    -c unix_socket_directories="${PG_SOCKET_DIR}" \
    -c hba_file="${PG_HBA_FILE}" \
    >"${LOG_DIR}/postgres.log" 2>&1 &
POSTGRES_PID=$!

# The superuser and the postgres database always exist, so the probe is not logged
# by the server as a failed login (as pg_isready's default root/root would be).
wait_for PostgreSQL "${POSTGRES_PID}" 60 "${LOG_DIR}/postgres.log" \
    "${PG_BIN}/pg_isready" -h "${PG_SOCKET_DIR}" -U "${PG_SUPERUSER}" -d postgres

# --- 2. migrations ----------------------------------------------------------
# dbmate creates the metamist database if it is missing, then applies whatever
# migrations are not yet recorded (a rerun on an up-to-date database does nothing).
# It connects as the superuser over the socket. The search_path is the one metamist's
# own db/docker-compose.yaml uses, so dbmate's bookkeeping table lands in public.
# Unlike liquibase there is no changelog lock to release after an interrupted run:
# dbmate's lock is a session lock, which PostgreSQL drops with the connection.
phase 2/6 "running dbmate migrations"
DBMATE_URL="postgres://${PG_SUPERUSER}@/${PG_DATABASE}?socket=${PG_SOCKET_DIR}&sslmode=disable&search_path=public,main,history"
run_step "${METAMIST_DIR}" dbmate --url "${DBMATE_URL}" \
    --migrations-dir "${METAMIST_DIR}/db/migrations" --no-dump-schema up

# --- 3. app login + local default user's project-creator rights (idempotent) -
phase 3/6 "enabling ${PG_APP_ROLE} login and granting ${SM_LOCALONLY_DEFAULTUSER} project-creators / members-admin"
# The migrations create metamist_server as NOLOGIN, because they cannot hold a
# password; here it just gets LOGIN, with no password (the pg_hba rules above are the
# access control). Its search_path is set so a client that connects with psql from
# the host sees metamist's tables without schema-qualifying them (the API sets
# the same path itself on every connection).
"${PG_BIN}/psql" -X -q -w -v ON_ERROR_STOP=1 -h "${PG_SOCKET_DIR}" -U "${PG_SUPERUSER}" -d "${PG_DATABASE}" <<SQL
ALTER ROLE ${PG_APP_ROLE} LOGIN;
ALTER ROLE ${PG_APP_ROLE} SET search_path = public, main, history;
SQL

# The grant is metamist's own, from its docs/installation.md, run as the app login
# over TCP: the same path the API takes, so a broken login fails the boot here.
"${PG_BIN}/psql" -X -q -w -v ON_ERROR_STOP=1 -v local_user="${SM_LOCALONLY_DEFAULTUSER}" \
    -h "${SM_DEV_DB_HOST}" -p "${SM_DEV_DB_PORT}" -U "${PG_APP_ROLE}" -d "${PG_DATABASE}" <<'SQL'
INSERT INTO main.group_member (group_id, member)
SELECT g.id, :'local_user'
FROM main."group" g
WHERE g.name IN ('project-creators', 'members-admin')
  AND NOT EXISTS (
      SELECT 1 FROM main.group_member gm
      WHERE gm.group_id = g.id AND gm.member = :'local_user'
  );
SQL

# --- 4. fake GCS emulator ---------------------------------------------------
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
phase 4/6 "starting fake-gcs-server on :4443 (filesystem backend at ${FAKEGCS_ROOT})"
fake-gcs-server -scheme http -host 0.0.0.0 -port 4443 \
    -backend filesystem -filesystem-root "${FAKEGCS_ROOT}" \
    -external-url "http://localhost:${GCS_HOST_PORT}" -public-host "localhost:${GCS_HOST_PORT}" \
    >"${LOG_DIR}/fakegcs.log" 2>&1 &
GCS_PID=$!

wait_for fake-gcs-server "${GCS_PID}" 30 "${LOG_DIR}/fakegcs.log" \
    curl -fs "http://localhost:4443/storage/v1/b?project=${GOOGLE_CLOUD_PROJECT}"

# --- 5. API server ----------------------------------------------------------
phase 5/6 "starting API on :8000 (uvicorn ${UVICORN_RELOAD[*]:-without reload})"
( cd "${METAMIST_DIR}" && exec "${VENV_PY}" -m uvicorn \
    --host 0.0.0.0 --port 8000 "${UVICORN_RELOAD[@]}" api.server:app ) \
    >"${LOG_DIR}/api.log" 2>&1 &
API_PID=$!

wait_for API "${API_PID}" 60 "${LOG_DIR}/api.log" curl -fs http://localhost:8000/openapi.json

# --- 6. seed once, if asked -------------------------------------------------
if [[ "${SEED}" == "0" ]]; then
    phase 6/6 "seed skipped (SEED=0)"
elif [[ -f "${SEED_MARKER}" ]]; then
    # Its content was checked against SEED_PATH before anything started.
    phase 6/6 "seed skipped (already seeded with ${SEED_PATH})"
else
    # A missing (checked before boot) or failing seed is fatal: the container never
    # reports ready with less data than was asked for.
    phase 6/6 "seeding: ${SEED_PATH}"
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
    Postgres   ${PG_APP_ROLE}@localhost:5432/${PG_DATABASE}  (no password; publish 5432 to reach it from the host)

  Query it:    mmquery '{ myProjects { id name } }'
  Logs:        ${LOG_DIR}/api.log  ${LOG_DIR}/postgres.log

BANNER

# Run the command as a child (not exec) so this script can stop the services when
# it ends. <&0 keeps an interactive command's stdin; a background job would
# otherwise get /dev/null.
"$@" <&0 &
CMD_PID=$!

# Block until the command or a service exits, by polling with kill -0 rather than
# `wait -n`. bash's wait -n looks for finished jobs, then makes a blocking wait; a
# child that exits between the two is reaped without being reported and wait -n
# then blocks for good (seen about once in several thousand runs, in bash 5.3 too),
# which here would leave the container up with its command already gone. kill -0
# and `wait <pid>` have no such gap. The sleep runs in the background and is waited
# for, so the TERM and INT traps fire at once instead of after the sleep.
exited_pid=
status=0
while [[ -z "${exited_pid}" ]]; do
    for pid in "${CMD_PID}" "${POSTGRES_PID}" "${GCS_PID}" "${API_PID}"; do
        if ! kill -0 "${pid}" 2>/dev/null; then
            exited_pid="${pid}"
            break
        fi
    done
    if [[ -z "${exited_pid}" ]]; then
        sleep 0.5 &
        SLEEP_PID=$!
        wait "${SLEEP_PID}" || true
        SLEEP_PID=
    fi
done
# The process is gone, so this returns its status at once (|| keeps a non-zero
# status from tripping set -e and the ERR trap).
wait "${exited_pid}" || status=$?

case "${exited_pid}" in
    "${CMD_PID}")
        log "command exited with status ${status}"
        exit "${status}"
        ;;
    "${POSTGRES_PID}")
        service=PostgreSQL logfile="${LOG_DIR}/postgres.log"
        ;;
    "${GCS_PID}")
        service=fake-gcs-server logfile="${LOG_DIR}/fakegcs.log"
        ;;
    "${API_PID}")
        service=API logfile="${LOG_DIR}/api.log"
        ;;
    *)
        die "unknown process ${exited_pid} exited (status ${status})"
        ;;
esac
log "${service} exited unexpectedly (status ${status})"
tail -n 50 "${logfile}" || true
exit 1
