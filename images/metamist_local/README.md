# metamist_local

Metamist in local mode, in one container, set up the way metamist's own
`docs/installation.md` describes. The image automates those steps so a script or
a CI job can have a working, empty metamist with no GCP credentials:

- a MariaDB 11.7 server, with the `sm_dev` database, `sm_api` user and role
- metamist's liquibase migrations, run at boot
- the metamist API (uvicorn, `--reload`), serving Swagger and GraphiQL
- the generated metamist python client, installed into the image's python
- fake-gcs-server, a local stand-in for Google Cloud Storage (the emulator
  metamist's own tests use)

It is credential-free: in local mode the API authenticates every request as
`SM_LOCALONLY_DEFAULTUSER` and the client sends no token. A `sitecustomize.py`
points the unmodified `google.cloud.storage.Client()` at the fake GCS with
anonymous credentials, and `hail` / `hailtop` resolve to import-only placeholders
(`/opt/hail-stubs`) so metamist scripts that import them still run.

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

### Which metamist branch

metamist's `dev` is its integration branch (deployed to its development
environment) and `main` receives release merges (deployed to production). Both
run on MariaDB, and their `docs/installation.md`, which this image automates, is
identical. The image pins a commit on **`main`**, so it matches what production
runs.

The move to Postgres lives on the separate `pg-dev` / `pg-main` branches and is
not released. Do not pin a commit from those branches: this image installs
MariaDB and runs its migrations. When Postgres is released to `main`, the
image's database layer and entrypoint change with it.

## Run

```bash
docker run -d --name metamist_local \
    -p 127.0.0.1:8000:8000 -p 127.0.0.1:4443:4443 \
    -v metamist_local_db:/var/lib/mysql \
    -v metamist_local_gcs:/data/fakegcs \
    australia-southeast1-docker.pkg.dev/cpg-common/images/metamist_local:<tag>
```

Then `http://localhost:8000/docs` is Swagger and `http://localhost:8000/graphql`
is GraphiQL. If you publish a container port on a different host port, tell the
container with `HOST_API_PORT` / `HOST_GCS_PORT` (below).

The API is unauthenticated and the database root user has no password, so
publish the ports on loopback only, as above.

`mmquery` is installed in the image: a small stdlib-only GraphQL client that
prints JSON and exits non-zero on GraphQL errors.

```bash
docker exec metamist_local mmquery '{ myProjects { name } }'
```

It also runs on the host (copy it out with `docker cp`); point it at the
published port with `SM_URL=http://localhost:8000`.

## Environment variables

Set these with `docker run -e`.

| Variable | Default | Meaning |
| --- | --- | --- |
| `SM_LOCALONLY_DEFAULTUSER` | `local-user` | The user every API call runs as. Added to the `project-creators` and `members-admin` groups at boot. Letters, digits, and `_ . @ -` only; anything else stops the container. |
| `SEED` | `0` | `1` runs `SEED_SCRIPT` once on first boot. |
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
| 3306 | MariaDB (`root` and `sm_api`, no password) |
| 4443 | fake GCS (`STORAGE_EMULATOR_HOST` for host-side clients) |

## Volumes and paths

The image declares no `VOLUME`: without `-v`, data lives in the container and goes
with it.

| Path | What |
| --- | --- |
| `/var/lib/mysql` | MariaDB data, and the `.metamist_seeded` marker file |
| `/data/fakegcs` | fake GCS objects (filesystem backend) |
| `/app/metamist` | `METAMIST_DIR`: optional metamist checkout to serve (see below) |
| `/opt/seed` | conventional place to mount or copy a seed script; nothing is baked in |
| `/var/log/metamist` | `api.log`, `mariadb.log`, `fakegcs.log` |

Keep `/var/lib/mysql` and `/data/fakegcs` together: mount both or neither, so the
uploaded objects stay in step with the database that references them.

## Seeding

Seeding is off by default; the stack boots empty. To load data on first boot:

```bash
docker run -d --name metamist_local \
    -v "$PWD/my_seed:/opt/seed" \
    -e SEED=1 -e SEED_SCRIPT=/opt/seed/generate.py \
    ...
```

- The script runs once, in phase 7, after the API is up, as
  `cd $METAMIST_DIR && /opt/venv/bin/python <script>`. The environment is already
  set for the generated python client, so the script can talk to the local API.
- An absolute `SEED_SCRIPT` runs as is; a relative one is resolved against
  `METAMIST_DIR`.
- On success the marker `/var/lib/mysql/.metamist_seeded` is written, and the seed
  does not run again while that volume persists.
- If the script is not found or exits non-zero, the container exits non-zero
  before logging `ready`, so it never serves less data than was asked for. No
  marker is written, so the next boot tries again.

To seed an already-running container instead, leave `SEED=0` and run the script
with the image's python:

```bash
docker cp ./my_seed metamist_local:/opt/seed
docker exec metamist_local python /opt/seed/generate.py
```

## Healthy

The `HEALTHCHECK` polls `GET /docs` every 10 seconds (180 second start period, 6
retries). Healthy means the API is serving. That happens before the seed step, so
a container that is seeding turns healthy while the seed is still running. If you
need the seed finished, wait for the `ready` log line instead.

```bash
docker inspect -f '{{.State.Health.Status}}' metamist_local
```

## Boot log

The entrypoint logs seven numbered phases and then a final line. These lines are an
interface: host tooling waits on them.

```text
[metamist_local 1/7] starting MariaDB
[metamist_local 2/7] ensuring database, user and role exist
[metamist_local 3/7] running liquibase migrations
[metamist_local 4/7] granting <user> project-creators / members-admin
[metamist_local 5/7] starting fake-gcs-server on :4443 (filesystem backend at /data/fakegcs)
[metamist_local 6/7] starting API on :8000 (uvicorn --reload)
[metamist_local 7/7] seeding: <script>
[metamist_local] ready
```

- On a successful boot each `[metamist_local N/7]` line appears once, in order, and
  `[metamist_local] ready` comes last, after phase 7 has finished. Only the prefix
  is the contract; the text after it can change.
- Phase 7 reads `seeding: <script>` or `seed skipped (already seeded, or SEED=0)`.
- Other lines with the plain `[metamist_local]` prefix are progress or errors. A
  fatal error is logged with that prefix and the container exits non-zero.

## Serving a metamist checkout (live reload)

To develop against metamist itself, mount a checkout over `METAMIST_DIR`:

```bash
docker run -d -v /path/to/metamist:/app/metamist ...
```

If `api/server.py` is found there, the migrations run from that tree's `db/` and
`uvicorn --reload` serves its `api/`, so edits to the python source are picked up
live. Otherwise the baked-in copy at `/build` is used. The python dependencies and
generated client are the ones baked at build time, so a checkout with different
dependencies, or a changed API model, needs an image built from it (next section).

## Building against a local checkout

CI builds with no build arguments, fetching metamist at the pinned commit. The
metamist source is a build stage named `metamist_src` whose root is the metamist
tree, so you can replace it with a checkout:

```bash
mkdir /tmp/metamist-clean
git -C /path/to/metamist archive HEAD | tar -x -C /tmp/metamist-clean
docker build --build-context metamist_src=/tmp/metamist-clean \
    -t metamist_local:dev images/metamist_local
```

Pass a clean copy, as `git archive` makes, not a working tree: a `.venv`,
`node_modules` or `web/` build output would be copied into the image. The tree
must have `pyproject.toml` and `uv.lock` at its root.

## Bumping

1. Pick a commit on metamist's `main` and set `ARG METAMIST_COMMIT` in the
   `Dockerfile`.
2. Set `ARG VERSION` to that commit's `version` in metamist's `pyproject.toml`. It
   must be the first `VERSION` line in the file: CI reads it to tag the image
   `<version>-N`.
3. CI rebuilds an image only when a `Dockerfile` changes. A change anywhere else in
   this folder (the entrypoint, the stubs, `mmquery`) has no effect until the
   `Dockerfile` is touched in the same change.
