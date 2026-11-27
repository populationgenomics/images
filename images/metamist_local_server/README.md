# metamist_local_server

Metamist in local mode, in one container: an automation of metamist's own
development setup (`docs/installation.md`). It gives a working, empty metamist
with no GCP credentials:

- a PostgreSQL 18 server with the `temporal_tables` extension, holding the
  `metamist` database
- metamist's dbmate migrations, run at boot
- the metamist API (uvicorn), serving Swagger and GraphiQL
- the generated metamist python client, installed into the image's python
- fake-gcs-server, a local stand-in for Google Cloud Storage (the emulator
  metamist's own tests use)

It is credential-free: in local mode the API authenticates every request as
`SM_LOCALONLY_DEFAULTUSER` and the client sends no token. A `sitecustomize.py`
points the unmodified `google.cloud.storage.Client()` at the fake GCS with
anonymous credentials, and an import hook in it serves `hail` / `hailtop` as
import-only placeholders, so metamist scripts that import them still run. Code
that tries to run hail raises a `RuntimeError`.

The React web UI is not built. The data is empty unless you seed it (see below).

## Versions

| Component | Version |
| --- | --- |
| metamist | 8.0.0, commit `785e3b02337850d736ad9864cc49d8f1d1897704` (metamist has no tags) |
| PostgreSQL | 18, newest point release from the PGDG apt repository |
| temporal_tables | 1.2.2, built from the PGXN release |
| dbmate | 2.38.1 |
| fake-gcs-server | 1.54.0 |
| Python | 3.14, dependencies from metamist's `uv.lock` (installed with uv 0.8.22 at build time; the image has no uv) |

The image tag is `8.0.0-N`: metamist's version plus CI's build counter.

Pinned to the 8.0.0 release commit, on metamist's `pg-main` branch (PostgreSQL).
Metamist 7.x, on the `dev` (default) and `main` branches, used MariaDB and is not
what this image runs.

## Upgrading from 7.14.x

Version 8 replaces MariaDB with PostgreSQL, so:

- **Old database volumes are incompatible.** A volume that held `/var/lib/mysql`
  cannot be used, at either path: the container refuses to start if it finds
  MariaDB files in `/var/lib/postgresql`. Remove it, and the fake GCS volume with
  it (the two are kept in step), then start again. The data is rebuilt by the seed.

  ```bash
  docker rm -f metamist_local_server
  docker volume rm metamist_local_server_db metamist_local_server_gcs
  ```

  Removing only the database volume is refused too: a new database does not start
  beside a fake GCS volume that already holds objects.
- The database volume is now mounted at `/var/lib/postgresql`, not `/var/lib/mysql`.
  The container refuses to start if `/var/lib/mysql` exists (a compose file that
  still mounts the old volume there), rather than ignore the mount and start empty.
- The database port is 5432, not 3306. Connect with a PostgreSQL client to the
  `metamist` database (formerly `sm_dev`) as `metamist_server` (formerly `sm_api`).
  The tables are in the `main` schema, which is on that user's `search_path`, so
  queries need no schema prefix (see [PostgreSQL](#postgresql)).
- The boot log has six numbered phases, not seven (see [Boot log](#boot-log)).
- `SM_DEV_DB_NAME` and `SM_DEV_DB_USER` can no longer be overridden: metamist's
  migrations fix them as `metamist` and `metamist_server`, and the container
  refuses to start with other values. It also refuses an `SM_DEV_DB_PORT` other
  than 5432, such as a leftover 3306.
- The log file `mariadb.log` is now `postgres.log`.
- The seed markers `.metamist_seeded` and `.metamist_seed_in_progress` moved with
  the data, from `/var/lib/mysql` to `/var/lib/postgresql/18/docker`. Their names
  and meaning are unchanged.
- A mounted metamist checkout must be version 8 or later (it needs `db/migrations`).

Everything else is as before: the other environment variables, the seed's rules
for refusing to start, the `ready` line, the health check, `mmquery`, and the
`--build-context` override.

## Pull

Authenticate docker to the registry once, with an account that can read the
`cpg-common` images repository:

```bash
gcloud auth configure-docker australia-southeast1-docker.pkg.dev
```

List the available tags, then pull one (the commands below show `<tag>`):

```bash
gcloud artifacts docker tags list \
    australia-southeast1-docker.pkg.dev/cpg-common/images/metamist_local_server
docker pull australia-southeast1-docker.pkg.dev/cpg-common/images/metamist_local_server:<tag>
```

## Run

```bash
docker run -d --name metamist_local_server \
    -p 127.0.0.1:8000:8000 -p 127.0.0.1:4443:4443 \
    -v metamist_local_server_db:/var/lib/postgresql \
    -v metamist_local_server_gcs:/data/fakegcs \
    australia-southeast1-docker.pkg.dev/cpg-common/images/metamist_local_server:<tag>
```

Then `http://localhost:8000/docs` is Swagger and `http://localhost:8000/graphql`
is GraphiQL. If you publish a container port on a different host port, tell the
container with `HOST_API_PORT` / `HOST_GCS_PORT` (below). To reach PostgreSQL
from the host as well, add `-p 127.0.0.1:5432:5432` (see [PostgreSQL](#postgresql)).

The API is unauthenticated and the database login `metamist_server` has no
password, so publish the ports on loopback only, as above.

Day to day:

```bash
docker logs -f metamist_local_server    # follow the boot log, then the API log
docker stop metamist_local_server       # clean shutdown; the data stays on the volumes
docker start metamist_local_server      # boot again on the same data
docker rm -f metamist_local_server      # remove the container; the volumes stay
```

To reset to an empty metamist, remove the container and both volumes:

```bash
docker rm -f metamist_local_server
docker volume rm metamist_local_server_db metamist_local_server_gcs
```

`mmquery` is installed in the image: a small stdlib-only GraphQL client that
prints JSON and exits non-zero on GraphQL errors.

```bash
docker exec metamist_local_server mmquery '{ myProjects { name } }'
```

It also runs on the host (copy it out with `docker cp`); point it at the
published port with `SM_URL=http://localhost:8000`.

`docker exec` starts in `/app/metamist`, which is empty when no checkout is
mounted there. A script that imports metamist's server modules (`api`, `db`) needs
`docker exec -w /build` (or the mount path); the installed `metamist` client
imports from anywhere.

## PostgreSQL

The database is `metamist`, with metamist's tables in the `main` schema and their
history in the `history` schema (maintained by `temporal_tables` triggers). It
listens on 5432 in the container. Access is set by `pg_hba.conf` rules that the
entrypoint writes on every boot, and no role has a password:

- `metamist_server`, the application login, is the only role that can connect over
  TCP, and only to the `metamist` database. It is not a superuser. This is the
  login the API uses, and the one to use from the host (with 5432 published on
  loopback). Its `search_path` is `public, main, history`, so tables need no schema
  prefix.
- `metamist_superuser`, which owns the schemas, can connect only over the local
  socket, that is from inside the container.

```bash
# from the host, with -p 127.0.0.1:5432:5432
psql "postgresql://metamist_server@localhost:5432/metamist"
# from inside the container, as the superuser
docker exec -it metamist_local_server psql -U metamist_superuser metamist
```

The role `metamist_server` is created by metamist's migrations without login (they
cannot hold a password). The entrypoint gives it `LOGIN` after the migrations
run, and grants `SM_LOCALONLY_DEFAULTUSER` the `project-creators` and
`members-admin` groups with the `INSERT INTO group_member ...` statement from
metamist's `docs/installation.md`.

## From tests (testcontainers)

A test suite can start the image itself with
[testcontainers](https://testcontainers-python.readthedocs.io/), waiting on the
`ready` log line (tested with testcontainers 4.15):

```python
import re

from testcontainers.core.container import DockerContainer
from testcontainers.core.wait_strategies import LogMessageWaitStrategy

metamist = (
    DockerContainer('australia-southeast1-docker.pkg.dev/cpg-common/images/metamist_local_server:<tag>')
    .with_exposed_ports(8000)
    .waiting_for(
        LogMessageWaitStrategy(re.compile(r'\[metamist_local_server\] ready'))
        .with_startup_timeout(300)
    )
)
with metamist:
    url = f'http://{metamist.get_container_host_ip()}:{metamist.get_exposed_port(8000)}'
    # point the client or tests at url (SM_URL=url, SM_ENVIRONMENT=local)
```

Exposed ports land on random host ports. That is fine for the API. If host code
uploads to the fake GCS, publish 4443 on a fixed port and pass it as
`HOST_GCS_PORT`, because the emulator returns upload URLs built from that port.

Unlike the `docker run` commands above, this publishes the port on every host
interface (`0.0.0.0` and `[::]`), not on loopback, and testcontainers 4.15 has no
typed way to choose the interface (`with_bind_ports` takes only a port). The API is
unauthenticated, so use this only on a machine whose ports other hosts cannot
reach, such as a laptop behind a firewall or a CI runner. To bind to loopback
anyway, set `metamist.ports['8000'] = ('127.0.0.1', 0)` before starting the
container. That works in 4.15 (docker-py accepts the tuple) but is outside
testcontainers' typed interface, which declares the values as ports only.

## Environment variables

Set these with `docker run -e`.

| Variable | Default | Meaning |
| --- | --- | --- |
| `SM_LOCALONLY_DEFAULTUSER` | `local-user` | The user every API call runs as. Added to the `project-creators` and `members-admin` groups at boot. Letters, digits, and `_ . @ -` only; anything else stops the container. |
| `SEED` | `0` | `1` runs `SEED_SCRIPT` once on first boot. Only `0` and `1` are accepted; anything else (such as `true`) stops the container. |
| `SEED_SCRIPT` | unset | The seed script to run. Required when `SEED=1`; the container refuses to start without it. |
| `HOST_API_PORT` | `8000` | The host port the API is published on. Used only in the boot banner. |
| `HOST_GCS_PORT` | `4443` | The host port fake GCS is published on. Used in the boot banner and as the external URL fake GCS puts in the links it hands out (such as resumable-upload sessions), so it must match what host clients connect to. |
| `METAMIST_DIR` | `/app/metamist` | Where to look for a mounted metamist checkout (see below). |

The image also sets metamist's own settings: `SM_ENVIRONMENT=local`,
`SM_URL=http://localhost:8000`, `SM_DEV_DB_NAME=metamist`,
`SM_DEV_DB_USER=metamist_server`, `SM_DEV_DB_HOST=127.0.0.1`,
`SM_DEV_DB_PORT=5432`, `STORAGE_EMULATOR_HOST=http://localhost:4443` and
`GOOGLE_CLOUD_PROJECT=metamist-local`. They match the services started inside the
container; leave them alone unless you know why you are changing them.
`SM_DEV_DB_NAME` and `SM_DEV_DB_USER` are the exception: metamist's migrations fix
those two names, so any other value stops the container.

Metamist's id prefixes (`SM_SAMPLEPREFIX`, `SM_SEQUENCINGGROUPPREFIX`,
`SM_COHORTPREFIX` and `SM_COHORTTEMPLATEPREFIX`, with the matching
`SM_*CHECKOFFSET` variables) are not set by the image, so metamist's defaults
apply. Pass them with `-e` to change them.

## Ports

| Port | Service |
| --- | --- |
| 8000 | metamist API: Swagger at `/docs`, GraphiQL at `/graphql`, REST at `/api` |
| 5432 | PostgreSQL: from the host only `metamist_server` (no password) can connect, and only to `metamist`; the superuser is local to the container |
| 4443 | fake GCS (`STORAGE_EMULATOR_HOST` for host-side clients) |

## Volumes and paths

The image declares no `VOLUME`: without `-v`, data lives in the container and goes
with it.

| Path | What |
| --- | --- |
| `/var/lib/postgresql` | PostgreSQL volume. The cluster is `PGDATA=/var/lib/postgresql/18/docker` (the layout of the official `postgres:18` image), and the seed markers `.metamist_seeded` and `.metamist_seed_in_progress` are inside it |
| `/data/fakegcs` | fake GCS objects (filesystem backend) |
| `/app/metamist` | `METAMIST_DIR`: optional metamist checkout to serve (see below) |
| `/opt/seed` | conventional place to mount or copy a seed script; nothing is baked in |
| `/var/log/metamist` | `api.log`, `postgres.log`, `initdb.log`, `fakegcs.log` |

Mount the database volume at `/var/lib/postgresql`, not at
`/var/lib/postgresql/data` (the path postgres images up to 17 used): nothing is
written there, so the container refuses to start if that path exists.

Keep `/var/lib/postgresql` and `/data/fakegcs` together: mount both or neither, so
the uploaded objects stay in step with the database that references them. The
container refuses to create a new database when `/data/fakegcs` already holds
objects.

Mount the database volume at `/var/lib/postgresql`, not at the data directory
inside it. On first boot the entrypoint runs `initdb` in a scratch directory
(`/var/lib/postgresql/18/docker.initdb`) and renames it to `PGDATA` when it has
finished, so a stop during `initdb` never leaves a half-written `PGDATA`: the next
boot removes the scratch directory and starts again. The container refuses to start
if `PGDATA` exists but is not a complete cluster (it lacks `PG_VERSION` or
`global/pg_control`, for example after a copy that was cut short), and never deletes
what is there; remove the volumes and start again.

## Seeding

Seeding is off by default; the stack boots empty. To load data on first boot:

```bash
docker run -d --name metamist_local_server \
    -p 127.0.0.1:8000:8000 -p 127.0.0.1:4443:4443 \
    -v metamist_local_server_db:/var/lib/postgresql \
    -v metamist_local_server_gcs:/data/fakegcs \
    -v "$PWD/my_seed:/opt/seed" \
    -e SEED=1 -e SEED_SCRIPT=/opt/seed/generate.py \
    australia-southeast1-docker.pkg.dev/cpg-common/images/metamist_local_server:<tag>
```

- The script runs once, in phase 6, after the API is up, as
  `cd $METAMIST_DIR && /opt/venv/bin/python <script>`. The environment is already
  set for the generated python client, so the script can talk to the local API.
- An absolute `SEED_SCRIPT` runs as is; a relative one is resolved against
  `METAMIST_DIR`. The resolved path is what the markers below record.
- If the script is not found, the container refuses to start before it starts any
  service (unless the volume is already seeded with that path, which skips the
  seed). If the script exits non-zero, the container exits non-zero before logging
  `ready`. Either way it never serves less data than was asked for.
- Metamist ships a synthetic seed: `-e SEED=1 -e SEED_SCRIPT=test/data/generate_data.py`
  creates the `greek-myth` project, with participants, samples, sequencing groups
  and analyses. It resolves to `/build/test/data/generate_data.py` in the baked source.

The seed state is kept as two marker files inside the PostgreSQL data directory
(`/var/lib/postgresql/18/docker`), each holding the resolved script path:

| Marker | Written | Meaning |
| --- | --- | --- |
| `.metamist_seed_in_progress` | before the script runs | the script started; removed when it succeeds |
| `.metamist_seeded` | after the script succeeds | the volume holds that script's data |

The container refuses to start, with a message that names the script, when:

- `.metamist_seed_in_progress` exists (whatever `SEED` is): the seed failed or was
  interrupted, so the data may be partial.
- `SEED=1` and `.metamist_seeded` holds a different path than the resolved
  `SEED_SCRIPT`.

The fix is to remove both volumes (`/var/lib/postgresql` and `/data/fakegcs`) and
start again, as under "Run". In the second case you can instead keep the data
already there by starting with `SEED=0`. If the marker matches the resolved
`SEED_SCRIPT`, a `SEED=1` boot skips the seed.

To seed an already-running container instead, leave `SEED=0` and run the script
with the image's python (in `/build`, as the entrypoint does):

```bash
docker cp ./my_seed metamist_local_server:/opt/seed
docker exec -w /build metamist_local_server python /opt/seed/generate.py
```

## Healthy

The `HEALTHCHECK` runs `curl` against `GET /docs` every 10 seconds (180 second
start period, 6 retries). Healthy means the API is serving. That happens before
the seed step, so a container that is seeding turns healthy while the seed is
still running. If you need the seed finished, wait for the `ready` log line instead.

```bash
docker inspect -f '{{.State.Health.Status}}' metamist_local_server
```

## Boot log

The entrypoint logs six numbered phases and then a final line. These lines are an
interface: host tooling waits on them.

```text
[metamist_local_server 1/6] starting PostgreSQL
[metamist_local_server 2/6] running dbmate migrations
[metamist_local_server 3/6] enabling metamist_server login and granting <user> project-creators / members-admin
[metamist_local_server 4/6] starting fake-gcs-server on :4443 (filesystem backend at /data/fakegcs)
[metamist_local_server 5/6] starting API on :8000 (uvicorn without reload)
[metamist_local_server 6/6] seeding: <script>
[metamist_local_server] ready
```

- On a successful boot each `[metamist_local_server N/6]` line appears once, in order, and
  `[metamist_local_server] ready` comes last, after phase 6 has finished. Only the prefix
  is the contract; the text after it can change (the 7.x images had seven phases,
  with MariaDB set-up spread over three of them).
- Phase 6 reads `seeding: <script>`, `seed skipped (SEED=0)` or
  `seed skipped (already seeded with <script>)`.
- Other lines with the plain `[metamist_local_server]` prefix are progress or errors. A
  fatal error is logged with that prefix and the container exits non-zero. An
  unexpected failure in a phase also logs `failed during phase N/6 ...`.
- A service that does not come up in time (PostgreSQL 60 s, fake GCS 30 s, API 60 s)
  logs `<service> did not come up`, followed by the last 50 lines of its log. A
  service whose process exits while the entrypoint is waiting for it logs
  `<service> exited before it came up (status N)` at once, with the same log
  tail, instead of waiting out the timeout.
- The first boot on an empty volume also logs `initialising PostgreSQL data
  directory` before phase 1.

## Shutdown

The entrypoint stays PID 1 and runs the command (by default `tail -F` of the API
log) as a child, so it can stop everything it started:

- On `docker stop` (SIGTERM) or Ctrl-C (SIGINT) it stops the command, the API,
  fake GCS and then PostgreSQL, which gets a clean (fast) shutdown. It exits with
  143 or 130. The next boot on the same volume needs no PostgreSQL crash recovery.
  This also applies to a stop in the middle of a boot.
- PostgreSQL gets 20 s for the fast shutdown. If it is still running after that, it
  gets an immediate shutdown (SIGQUIT, with no checkpoint, so the next boot runs
  crash recovery) and, 5 s later, SIGKILL. `docker stop` waits only 10 s before
  it kills the container, so a database with a long checkpoint to write needs
  `docker stop -t 30`.
- When the command exits, the services are stopped and the container exits with
  the command's status, so `docker run --rm <image> true` exits 0.
- When a service exits on its own, the container logs `<service> exited
  unexpectedly`, shows the end of that service's log, stops the rest and exits
  non-zero. The entrypoint checks every half second whether the command or a
  service is still running.
- A fatal error during boot stops whatever was started before it exits.

Ctrl-C under `docker run -it` stops the container; it does not cancel what an
interactive command is doing. The entrypoint starts the command and the services as
background jobs that start with SIGINT ignored and share its process group, so the
terminal's Ctrl-C reaches the entrypoint, which shuts everything down. To use an
interactive client such as `psql`, run it in the container with `docker exec -it`
instead, where Ctrl-C reaches only that process.

## Serving a metamist checkout (live reload)

To develop against metamist itself, mount a checkout over `METAMIST_DIR`:

```bash
docker run -d -v /path/to/metamist:/app/metamist ...
```

If `api/server.py` is found there, the migrations run from that tree's
`db/migrations` and `uvicorn --reload` serves its `api/`, so edits to the python
source are picked up live. Otherwise the baked-in copy at `/build` is served
without reload. The checkout must be metamist 8 or later, which today means the
`pg-main` branch (a fresh clone is on `dev`, which is 7.x): one without
`db/migrations` stops the container. The python dependencies and generated client
are the ones baked at build time, so a checkout with different dependencies, or a
changed API model, needs an image built from it (see
[Building against a local checkout](#building-against-a-local-checkout)).

## Building

From the repository root, with BuildKit (the default builder in current Docker;
the `Dockerfile` uses `RUN --mount`):

```bash
docker build -t metamist_local_server:dev images/metamist_local_server
```

This fetches metamist at the pinned commit, as CI does. The release workflow
(`deploy.yaml`, on a push to `main`) builds with no build args, so `ARG VERSION`
keeps the `Dockerfile`'s default (the image tag, `8.0.0-N`, comes from the
workflow). The manual `deploy_dev.yaml` and `deploy_container.yaml` workflows pass
`--build-arg VERSION=<the tag you give>`, so that tag must be `8.0.0` or
`8.0.0-<suffix>`: the build fails unless the part before the first `-` equals the
`version` in the tree's `pyproject.toml`.

## Building against a local checkout

The metamist source is a build stage named `metamist_src` whose root is the
metamist tree, so you can replace it with a checkout. Pass
`--build-arg VERSION=<that checkout's version>`: the same check applies.

```bash
mkdir /tmp/metamist-clean
git -C /path/to/metamist archive HEAD | tar -x -C /tmp/metamist-clean
docker build --build-context metamist_src=/tmp/metamist-clean \
    --build-arg VERSION=<version in /tmp/metamist-clean/pyproject.toml> \
    -t metamist_local_server:dev images/metamist_local_server
```

Pass a clean copy, as `git archive` makes, not a working tree: a `.venv`,
`node_modules` or `web/` build output would be copied into the image. The tree
must have `pyproject.toml` and `uv.lock` at its root.

The generated client comes from the code generators in metamist's locked `dev`
dependency group. They run in a build stage that has a JRE (the generator is a
Java program), and only the resulting python env is copied into the final image:
`/opt/venv` ends up with the locked runtime dependencies plus the client, and the
image has no JRE.

## Bumping

1. Pick a commit on metamist's `pg-main` (or whichever branch then carries
   releases) and set `ARG METAMIST_COMMIT` in the `Dockerfile`.
2. Set `ARG VERSION` to that commit's `version` in metamist's `pyproject.toml`. It
   must be the first `VERSION` line in the file: CI reads it to tag the image
   `<version>-N`. The build compares it (without the `-N`) with `/build/pyproject.toml`
   and fails on a mismatch, so a stale `ARG VERSION` cannot be published.
3. Follow metamist's `db/Dockerfile` for the database: when it moves to a new
   PostgreSQL major version, change `ARG PG_MAJOR` (this renames the data directory,
   so existing volumes need resetting), and when it changes the way `temporal_tables`
   or `dbmate` is installed, mirror that.
4. When bumping `temporal_tables`, `dbmate` or `fake-gcs-server`, update the
   matching `*_SHA256` build args in the `Dockerfile`. Cross-check each new hash
   against the publisher's own checksum, as the comments there describe: the
   SHA-1 in the release's `META.json` on PGXN for `temporal_tables`, the
   per-asset digest on the GitHub release page for `dbmate` and `fake-gcs-server`
   (for the latter it also matches the release's `checksums.txt`).
5. CI rebuilds an image only when a `Dockerfile` changes. A change anywhere else in
   this folder (the entrypoint, `sitecustomize.py`, `mmquery`) has no effect until the
   `Dockerfile` is touched in the same change.
6. Update this README to match: the [Versions](#versions) table, the `8.0.0` in the
   tag and [Building](#building) text, and, after a PostgreSQL major version
   change, the `18` in the data paths and in the summary at the top. Searching
   this file for the old version numbers finds them all. If releases have moved
   off `pg-main`, also update the branch names here, in the `Dockerfile` comment
   above `ARG METAMIST_COMMIT` and in the entrypoint's `db/migrations not found`
   message.
