# metamist_local_server

Metamist in local mode, in one container: an automation of metamist's own
development setup (`docs/installation.md`). It gives a working, empty metamist
with no GCP credentials:

- a MariaDB 11.7 server, with the `sm_dev` database, `sm_api` user and role
- metamist's liquibase migrations, run at boot
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
| metamist | 7.14.3, commit `ebfa456da4f9f47a3d0a8ca8cfb2feeb81abd242` (metamist has no tags) |
| MariaDB | 11.7 |
| liquibase | 4.26.0, with MariaDB JDBC driver 3.0.3 |
| fake-gcs-server | 1.54.0 |
| Python | 3.11, dependencies from metamist's `uv.lock` (uv 0.8.22) |

The image tag is `7.14.3-N`: metamist's version plus CI's build counter.

Pinned to a commit on metamist's `main` (MariaDB). Not `pg-dev` / `pg-main`:
the Postgres line is unreleased, and this image installs MariaDB.

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
    -v metamist_local_server_db:/var/lib/mysql \
    -v metamist_local_server_gcs:/data/fakegcs \
    australia-southeast1-docker.pkg.dev/cpg-common/images/metamist_local_server:<tag>
```

Then `http://localhost:8000/docs` is Swagger and `http://localhost:8000/graphql`
is GraphiQL. If you publish a container port on a different host port, tell the
container with `HOST_API_PORT` / `HOST_GCS_PORT` (below).

The API is unauthenticated and the database user `sm_api` has no password, so
publish the ports on loopback only, as above.

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
`SM_URL=http://localhost:8000`, `SM_DEV_DB_NAME=sm_dev`, `SM_DEV_DB_USER=sm_api`,
`SM_DEV_DB_HOST=127.0.0.1`, `SM_DEV_DB_PORT=3306`,
`STORAGE_EMULATOR_HOST=http://localhost:4443` and
`GOOGLE_CLOUD_PROJECT=metamist-local`. They match the services started inside the
container; leave them alone unless you know why you are changing them.

## Ports

| Port | Service |
| --- | --- |
| 8000 | metamist API: Swagger at `/docs`, GraphiQL at `/graphql`, REST at `/api` |
| 3306 | MariaDB: from the host only `sm_api` (no password) can connect; `root` is local to the container |
| 4443 | fake GCS (`STORAGE_EMULATOR_HOST` for host-side clients) |

## Volumes and paths

The image declares no `VOLUME`: without `-v`, data lives in the container and goes
with it.

| Path | What |
| --- | --- |
| `/var/lib/mysql` | MariaDB data, and the seed markers `.metamist_seeded` and `.metamist_seed_in_progress` |
| `/data/fakegcs` | fake GCS objects (filesystem backend) |
| `/app/metamist` | `METAMIST_DIR`: optional metamist checkout to serve (see below) |
| `/opt/seed` | conventional place to mount or copy a seed script; nothing is baked in |
| `/var/log/metamist` | `api.log`, `mariadb.log`, `fakegcs.log` |

Keep `/var/lib/mysql` and `/data/fakegcs` together: mount both or neither, so the
uploaded objects stay in step with the database that references them.

## Seeding

Seeding is off by default; the stack boots empty. To load data on first boot:

```bash
docker run -d --name metamist_local_server \
    -p 127.0.0.1:8000:8000 -p 127.0.0.1:4443:4443 \
    -v metamist_local_server_db:/var/lib/mysql \
    -v metamist_local_server_gcs:/data/fakegcs \
    -v "$PWD/my_seed:/opt/seed" \
    -e SEED=1 -e SEED_SCRIPT=/opt/seed/generate.py \
    australia-southeast1-docker.pkg.dev/cpg-common/images/metamist_local_server:<tag>
```

- The script runs once, in phase 7, after the API is up, as
  `cd $METAMIST_DIR && /opt/venv/bin/python <script>`. The environment is already
  set for the generated python client, so the script can talk to the local API.
- An absolute `SEED_SCRIPT` runs as is; a relative one is resolved against
  `METAMIST_DIR`. The resolved path is what the markers below record.
- If the script is not found or exits non-zero, the container exits non-zero
  before logging `ready`, so it never serves less data than was asked for.
- Metamist ships a synthetic seed: `-e SEED=1 -e SEED_SCRIPT=test/data/generate_data.py`
  creates the `greek-myth` project, with participants, samples, sequencing groups
  and analyses. It resolves to `/build/test/data/generate_data.py` in the baked source.

The seed state is kept as two marker files in `/var/lib/mysql`, each holding the
resolved script path:

| Marker | Written | Meaning |
| --- | --- | --- |
| `.metamist_seed_in_progress` | before the script runs | the script started; removed when it succeeds |
| `.metamist_seeded` | after the script succeeds | the volume holds that script's data |

The container refuses to start, with a message that names the script, when:

- `.metamist_seed_in_progress` exists (whatever `SEED` is): the seed failed or was
  interrupted, so the data may be partial.
- `SEED=1` and `.metamist_seeded` holds a different path than the resolved
  `SEED_SCRIPT`, including the empty marker that older images wrote.

The fix is to remove both volumes (`/var/lib/mysql` and `/data/fakegcs`) and start
again, as under "Run". In the second case you can instead keep the data already
there by starting with `SEED=0`. If the marker matches the resolved `SEED_SCRIPT`,
a `SEED=1` boot skips the seed.

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

The entrypoint logs seven numbered phases and then a final line. These lines are an
interface: host tooling waits on them.

```text
[metamist_local_server 1/7] starting MariaDB
[metamist_local_server 2/7] ensuring database, user and role exist
[metamist_local_server 3/7] running liquibase migrations
[metamist_local_server 4/7] granting <user> project-creators / members-admin
[metamist_local_server 5/7] starting fake-gcs-server on :4443 (filesystem backend at /data/fakegcs)
[metamist_local_server 6/7] starting API on :8000 (uvicorn without reload)
[metamist_local_server 7/7] seeding: <script>
[metamist_local_server] ready
```

- On a successful boot each `[metamist_local_server N/7]` line appears once, in order, and
  `[metamist_local_server] ready` comes last, after phase 7 has finished. Only the prefix
  is the contract; the text after it can change.
- Phase 7 reads `seeding: <script>`, `seed skipped (SEED=0)` or
  `seed skipped (already seeded with <script>)`.
- Other lines with the plain `[metamist_local_server]` prefix are progress or errors. A
  fatal error is logged with that prefix and the container exits non-zero. An
  unexpected failure in a phase also logs `failed during phase N/7 ...`.
- A service that does not come up in time (MariaDB 60 s, fake GCS 30 s, API 60 s)
  logs `<service> did not come up`, followed by the last 50 lines of its log.

## Shutdown

The entrypoint stays PID 1 and runs the command (by default `tail -F` of the API
log) as a child, so it can stop everything it started:

- On `docker stop` (SIGTERM) or Ctrl-C (SIGINT) it stops the command, the API,
  fake GCS and then MariaDB, which gets a clean shutdown. It exits with 143 or 130.
  The next boot on the same volume needs no MariaDB crash recovery. This also
  applies to a stop in the middle of a boot.
- When the command exits, the services are stopped and the container exits with
  the command's status, so `docker run --rm <image> true` exits 0.
- When a service exits on its own, the container logs `<service> exited
  unexpectedly`, shows the end of that service's log, stops the rest and exits
  non-zero.
- A fatal error during boot stops whatever was started before it exits.

## Serving a metamist checkout (live reload)

To develop against metamist itself, mount a checkout over `METAMIST_DIR`:

```bash
docker run -d -v /path/to/metamist:/app/metamist ...
```

If `api/server.py` is found there, the migrations run from that tree's `db/` and
`uvicorn --reload` serves its `api/`, so edits to the python source are picked up
live. Otherwise the baked-in copy at `/build` is served without reload. The python dependencies and
generated client are the ones baked at build time, so a checkout with different
dependencies, or a changed API model, needs an image built from it (next section).

## Building against a local checkout

CI builds with `--build-arg VERSION=<tag>` (`7.14.3-N`), fetching metamist at the
pinned commit. The metamist source is a build stage named `metamist_src` whose
root is the metamist tree, so you can replace it with a checkout. Pass
`--build-arg VERSION=<that checkout's version>`: the build fails unless it equals
the `version` in the tree's `pyproject.toml`.

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
dependency group. They are used during the build only: `/opt/venv` ends up with
the locked runtime dependencies plus the client.

## Bumping

1. Pick a commit on metamist's `main` and set `ARG METAMIST_COMMIT` in the
   `Dockerfile`.
2. Set `ARG VERSION` to that commit's `version` in metamist's `pyproject.toml`. It
   must be the first `VERSION` line in the file: CI reads it to tag the image
   `<version>-N`. The build compares it (without the `-N`) with `/build/pyproject.toml`
   and fails on a mismatch, so a stale `ARG VERSION` cannot be published.
3. When bumping liquibase, the JDBC driver or the MariaDB keyring, update the
   matching `*_SHA256` build arg in the `Dockerfile`. Cross-check the new hash
   against the publisher's own checksum, as the comment there describes.
4. CI rebuilds an image only when a `Dockerfile` changes. A change anywhere else in
   this folder (the entrypoint, `sitecustomize.py`, `mmquery`) has no effect until the
   `Dockerfile` is touched in the same change.
