# syntax=docker/dockerfile:1.7
#
# metamist_local: metamist in local mode, set up as its docs/installation.md
# describes, in one container:
#   - MariaDB 11.7 server (data on a volume)
#   - liquibase migrations (Java + MariaDB JDBC driver)
#   - fake-gcs-server, a local stand-in for Google Cloud Storage
#   - the FastAPI/GraphQL API (uvicorn), serving Swagger and GraphiQL (no web UI)
#   - the generated metamist python client, installed into the venv
#
# The metamist source is baked in at /build (deps + client). At runtime you can
# bind-mount a metamist working tree over /app/metamist (METAMIST_DIR) so edits
# are picked up live by `uvicorn --reload`. Nothing here needs GCP credentials:
# in local mode the server falls back to SM_LOCALONLY_DEFAULTUSER and the client
# sends no token.
#
# See README.md in this folder for the runtime contract.

# CI tags the image from this version (the first such line in this file) plus an
# auto-incremented suffix. It is metamist's pyproject version at METAMIST_COMMIT;
# bump the two together.
ARG VERSION=7.14.3
# metamist has no git tags or releases, so the source is pinned by a commit on
# its main (production, MariaDB) branch; not pg-dev/pg-main (Postgres, unreleased).
ARG METAMIST_COMMIT=ebfa456da4f9f47a3d0a8ca8cfb2feeb81abd242

# ---------------------------------------------------------------------------
# Fetch the metamist source at the pinned commit, without .git and without the
# React web UI (web/), which is intentionally not built: the API serves Swagger
# and GraphiQL directly, which is all that query/validation testing needs.
# ---------------------------------------------------------------------------
FROM alpine:3.21 AS metamist_fetch
ARG METAMIST_COMMIT
WORKDIR /src
# hadolint ignore=DL3018
RUN apk add --no-cache git \
    && git init --quiet . \
    && git fetch --quiet --depth 1 https://github.com/populationgenomics/metamist.git "${METAMIST_COMMIT}" \
    && git checkout --quiet FETCH_HEAD \
    && test "$(git rev-parse HEAD)" = "${METAMIST_COMMIT}" \
    && rm -rf .git web

# The metamist source tree as the filesystem root. The main stage copies only
# from this stage, so a local checkout replaces it without editing this file:
#   docker build --build-context metamist_src=/path/to/clean/metamist images/metamist_local
FROM scratch AS metamist_src
COPY --from=metamist_fetch /src/ /

# ---------------------------------------------------------------------------
# The image
# ---------------------------------------------------------------------------
FROM python:3.11-slim-bookworm

# pipefail so a failed `curl | bash` or checksum pipeline fails the build instead
# of carrying on (with the MariaDB repo missing, apt would install the wrong server).
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    UV_PROJECT_ENVIRONMENT=/opt/venv \
    UV_LINK_MODE=copy \
    UV_COMPILE_BYTECODE=1 \
    PATH="/opt/venv/bin:/usr/local/bin:${PATH}" \
    LIQUIBASE_VERSION=4.26.0 \
    MARIADB_JDBC_VERSION=3.0.3

# System packages: MariaDB 11.7 server/client, a JRE (for liquibase and
# the openapi generator), and the usual build/runtime helpers.
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl wget gnupg git unzip jq gosu procps \
        default-jre-headless \
    && curl -fsSL https://r.mariadb.com/downloads/mariadb_repo_setup \
        | bash -s -- --mariadb-server-version="mariadb-11.7" \
    && apt-get update && apt-get install -y --no-install-recommends \
        mariadb-server mariadb-client \
    && rm -rf /var/lib/apt/lists/*

# liquibase + the MariaDB JDBC driver it needs to run the migrations
RUN wget -q "https://github.com/liquibase/liquibase/releases/download/v${LIQUIBASE_VERSION}/liquibase-${LIQUIBASE_VERSION}.zip" \
    && unzip -q "liquibase-${LIQUIBASE_VERSION}.zip" -d /opt/liquibase \
    && rm "liquibase-${LIQUIBASE_VERSION}.zip" \
    && ln -s /opt/liquibase/liquibase /usr/local/bin/liquibase \
    && wget -q "https://repo1.maven.org/maven2/org/mariadb/jdbc/mariadb-java-client/${MARIADB_JDBC_VERSION}/mariadb-java-client-${MARIADB_JDBC_VERSION}.jar" \
        -O /opt/mariadb-java-client.jar

# fake-gcs-server: a local stand-in for Google Cloud Storage, the same emulator
# metamist's own test suite uses. metamist builds a storage.Client() whenever an
# analysis has gs:// outputs; this lets that work locally with no GCP credentials.
# Pinned and SHA-256 verified against the upstream checksums file.
ARG FAKE_GCS_VERSION=1.54.0
RUN arch="$(dpkg --print-architecture)" \
    && asset="fake-gcs-server_${FAKE_GCS_VERSION}_Linux_${arch}.tar.gz" \
    && base="https://github.com/fsouza/fake-gcs-server/releases/download/v${FAKE_GCS_VERSION}" \
    && wget -q "${base}/${asset}" -O /tmp/fake-gcs.tar.gz \
    && wget -q "${base}/checksums.txt" -O /tmp/fake-gcs-checksums.txt \
    && (cd /tmp && grep " ${asset}\$" fake-gcs-checksums.txt | sed "s| ${asset}| /tmp/fake-gcs.tar.gz|" | sha256sum -c -) \
    && tar -xzf /tmp/fake-gcs.tar.gz -C /usr/local/bin fake-gcs-server \
    && rm /tmp/fake-gcs.tar.gz /tmp/fake-gcs-checksums.txt \
    && chmod +x /usr/local/bin/fake-gcs-server

# uv (for the python env)
COPY --from=ghcr.io/astral-sh/uv:0.8.22 /uv /uvx /bin/

# ---------------------------------------------------------------------------
# Bake the metamist source, the python env and the generated client.
# ---------------------------------------------------------------------------
WORKDIR /build

# Dependency manifests first, so the (slow) dependency layer is reused when only
# source changed; the full tree lands in the next layer.
COPY --from=metamist_src /pyproject.toml /uv.lock /build/
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen --no-dev --no-install-project

COPY --from=metamist_src / /build/
RUN --mount=type=cache,target=/root/.cache/uv uv sync --frozen --no-dev

# Generate the python client (openapi-generator-cli is a build-only tool), then
# install the client into site-packages as a real (non-editable) package so
# `import metamist` works regardless of what is mounted over the source tree.
# openapi-generator-cli runs the generator; strawberry's CLI is used by
# regenerate_api.py to export the GraphQL schema. Both are build-only.
RUN --mount=type=cache,target=/root/.cache/uv \
    uv pip install --python /opt/venv openapi-generator-cli==5.3.0 'strawberry-graphql[cli]==0.268.1' 'libcst==1.5.1' \
    && OPENAPI_COMMAND=openapi-generator-cli python regenerate_api.py \
    && uv pip install --python /opt/venv --no-deps /build/packages/metamist

# Make the unmodified storage.Client() use anonymous creds against the emulator,
# and (local mode) resolve hail/hailtop to placeholders so hail-importing scripts
# run without the heavy dependency.
COPY sitecustomize.py /opt/venv/lib/python3.11/site-packages/sitecustomize.py
COPY hail-stubs/ /opt/hail-stubs/

# No seed is baked in: seeding is opt-in and project-specific. With SEED=1 the
# entrypoint runs SEED_SCRIPT once on first boot (mount or copy the script in;
# /opt/seed is the conventional place). Without it the stack boots empty.
COPY docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
COPY mmquery /usr/local/bin/mmquery
RUN chmod +x /usr/local/bin/docker-entrypoint.sh /usr/local/bin/mmquery \
    && mkdir -p /var/log/metamist /app/metamist \
    && chown -R mysql:mysql /var/lib/mysql

# Local-mode defaults. The default user is granted project-creators by the
# entrypoint, so the same value is used for auth and for the grant.
ENV SM_ENVIRONMENT=local \
    SM_URL=http://localhost:8000 \
    SM_LOCALONLY_DEFAULTUSER=local-user \
    SM_DEV_DB_NAME=sm_dev \
    SM_DEV_DB_USER=sm_api \
    SM_DEV_DB_HOST=127.0.0.1 \
    SM_DEV_DB_PORT=3306 \
    METAMIST_DIR=/app/metamist \
    STORAGE_EMULATOR_HOST=http://localhost:4443 \
    GOOGLE_CLOUD_PROJECT=metamist-local \
    SEED=0

WORKDIR /app/metamist
EXPOSE 8000 3306 4443

# Healthy = the API answers on :8000. That happens before the seed step runs, so
# a seeding container turns healthy while the seed is still loading; wait for the
# "[metamist_local] ready" log line if you need the seed finished. CI can wait on
# `docker inspect --format '{{.State.Health.Status}}'` / `--wait` /
# `condition: service_healthy` instead of scraping the boot log.
HEALTHCHECK --interval=10s --timeout=5s --start-period=180s --retries=6 \
    CMD ["python", "-c", "import urllib.request; urllib.request.urlopen('http://localhost:8000/docs', timeout=4)"]

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
# Keep the container alive and stream the API log when no command is given.
CMD ["tail", "-F", "/var/log/metamist/api.log"]
