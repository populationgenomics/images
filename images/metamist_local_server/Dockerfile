# syntax=docker/dockerfile:1.7
#
# metamist_local_server: metamist in local mode, set up as its docs/installation.md
# describes, in one container:
#   - PostgreSQL 18 server with the temporal_tables extension (data on a volume)
#   - dbmate migrations (metamist's db/migrations, run at boot)
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
# bump the two together (the build fails if they disagree).
ARG VERSION=8.0.0
# metamist has no git tags or releases, so the source is pinned by a commit. This
# one is the 8.0.0 release, on metamist's pg-main (PostgreSQL) branch; the 7.x
# line on dev and main used MariaDB, which this image no longer installs.
ARG METAMIST_COMMIT=785e3b02337850d736ad9864cc49d8f1d1897704
# The PostgreSQL major version: server, client, extension build and data path.
ARG PG_MAJOR=18

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

# The metamist source tree as the filesystem root. The build stages copy only
# from this stage, so a local checkout replaces it without editing this file:
#   docker build --build-context metamist_src=/path/to/clean/metamist images/metamist_local_server
FROM scratch AS metamist_src
COPY --from=metamist_fetch /src/ /

# ---------------------------------------------------------------------------
# The python base with the PostgreSQL apt repository (PGDG) configured. The
# PGDG signing key is the one Debian's own postgresql-common package ships, so
# it is verified through Debian's apt signatures and there is no key download to
# pin. create_main_cluster is turned off before any server package is installed:
# otherwise the package initialises a cluster under /var/lib/postgresql, which a
# fresh named volume would then inherit. The entrypoint runs initdb on first boot.
# ---------------------------------------------------------------------------
FROM python:3.14-slim-trixie AS pg_apt

# pipefail so a failed download or checksum pipeline fails the build.
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ENV DEBIAN_FRONTEND=noninteractive

# hadolint ignore=DL3008
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates postgresql-common \
    && sed -ri 's/^#?(create_main_cluster) .*$/\1 = false/' /etc/postgresql-common/createcluster.conf \
    && grep -q '^create_main_cluster = false$' /etc/postgresql-common/createcluster.conf \
    && echo "deb [signed-by=/usr/share/postgresql-common/pgdg/apt.postgresql.org.asc] https://apt.postgresql.org/pub/repos/apt trixie-pgdg main" \
        >/etc/apt/sources.list.d/pgdg.list \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# Build the temporal_tables extension (metamist's migrations use its
# versioning() trigger for the history tables) against the server headers, so
# the compiler and headers stay out of the final image. PGXN, the extension's
# own distribution channel, is what metamist's db/Dockerfile installs from; here
# the release is pinned and its SHA-256 recorded. The hash was cross-checked
# against the SHA-1 in the release's META.json on PGXN. Bump a hash together
# with its version.
# ---------------------------------------------------------------------------
FROM pg_apt AS temporal_tables_build
ARG PG_MAJOR
ARG TEMPORAL_TABLES_VERSION=1.2.2
ARG TEMPORAL_TABLES_SHA256=81f99ee6d5ad1a253c18e3891da0c3bdf699eaabae26e50b4686d8dfcb2f18ba

# hadolint ignore=DL3008
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl unzip gcc make libc6-dev "postgresql-server-dev-${PG_MAJOR}" \
    && curl -fsSL "https://api.pgxn.org/dist/temporal_tables/${TEMPORAL_TABLES_VERSION}/temporal_tables-${TEMPORAL_TABLES_VERSION}.zip" \
        -o /tmp/temporal_tables.zip \
    && echo "${TEMPORAL_TABLES_SHA256}  /tmp/temporal_tables.zip" | sha256sum -c - \
    && unzip -q /tmp/temporal_tables.zip -d /tmp \
    && make -C "/tmp/temporal_tables-${TEMPORAL_TABLES_VERSION}" \
    && make -C "/tmp/temporal_tables-${TEMPORAL_TABLES_VERSION}" install DESTDIR=/out

# ---------------------------------------------------------------------------
# Build the python env and the generated client. The client generator
# (openapi-generator-cli, from metamist's dev dependency group) runs on Java, so
# the JRE lives only in this stage: the final image gets just /opt/venv and /build.
# ---------------------------------------------------------------------------
FROM python:3.14-slim-trixie AS client_build

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# UV_PYTHON_DOWNLOADS=never: use this image's python (metamist pins 3.14) or fail,
# rather than have uv fetch another one.
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    UV_PROJECT_ENVIRONMENT=/opt/venv \
    UV_LINK_MODE=copy \
    UV_COMPILE_BYTECODE=1 \
    UV_PYTHON_DOWNLOADS=never \
    PATH="/opt/venv/bin:${PATH}"

# hadolint ignore=DL3008
RUN apt-get update && apt-get install -y --no-install-recommends \
        default-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# uv (for the python env)
COPY --from=ghcr.io/astral-sh/uv:0.8.22 /uv /uvx /bin/

WORKDIR /build

# Dependency manifests first, so the (slow) dependency layer is reused when only
# source changed; the full tree lands in the next layer.
COPY --from=metamist_src /pyproject.toml /uv.lock /build/
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen --no-dev --no-install-project

COPY --from=metamist_src / /build/

# The image tag (VERSION: `<version>` or CI's `<version>-N`) must name the metamist
# source baked in, so a bumped METAMIST_COMMIT or a --build-context checkout cannot
# leave a stale tag.
ARG VERSION
RUN src_version="$(python -c 'import tomllib; print(tomllib.load(open("/build/pyproject.toml", "rb"))["project"]["version"])')" \
    && if [ "${VERSION%%-*}" != "${src_version}" ]; then \
        echo "VERSION ${VERSION} does not match metamist ${src_version} in /build/pyproject.toml" >&2; \
        exit 1; \
    fi

# Generate the python client with the code generators from metamist's own dev
# group (openapi-generator-cli runs the generator; strawberry's CLI exports the
# GraphQL schema), then sync back to the runtime set so /opt/venv holds exactly
# the locked runtime dependencies. Install the client into site-packages as a real
# (non-editable) package so `import metamist` works regardless of what is mounted
# over the source tree. Last, remove /build/web: the source has no web/ (see
# metamist_fetch), so it holds only what regenerate_api.py recreated there (the
# generated TypeScript client and a copied image), which nothing serves.
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen \
    && OPENAPI_COMMAND=openapi-generator-cli python regenerate_api.py \
    && uv sync --frozen --no-dev \
    && uv pip install --python /opt/venv --no-deps /build/packages/metamist \
    && rm -rf /build/web

# ---------------------------------------------------------------------------
# The image
# ---------------------------------------------------------------------------
FROM pg_apt
ARG PG_MAJOR

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    LANG=en_US.utf8 \
    PG_MAJOR=${PG_MAJOR} \
    PGDATA=/var/lib/postgresql/${PG_MAJOR}/docker \
    PATH="/opt/venv/bin:/usr/local/bin:${PATH}:/usr/lib/postgresql/${PG_MAJOR}/bin"

# dbmate, the migration runner metamist uses (its docs run `dbmate up`). Pinned to
# a release; each binary's SHA-256 is the digest GitHub publishes for that release
# asset (checked against the downloaded files). Bump a hash together with the
# version. (An ARG, not an ENV: dbmate reads DBMATE_* environment variables as options.)
ARG DBMATE_VERSION=2.38.1
ARG DBMATE_SHA256_AMD64=3e2c62735e0e63f3e5e4de6a3fa5c8a29ac9a34dcb8631a08d9589ab943ad9bd
ARG DBMATE_SHA256_ARM64=bb7910069907621feb96e75d26676c55a3bf3514dad7d9c4ca69dfd8465c0543

# System packages: the PostgreSQL server and client (the newest point release of
# the pinned major: the PGDG repo drops old ones), and the usual runtime
# helpers. curl also downloads the binaries below, so wget is not installed.
# en_US.UTF-8 is generated and made the default because that is the locale
# metamist's own db image (postgres:18) initialises its cluster with, so text
# ordering matches what metamist's tests and queries expect.
# hadolint ignore=DL3008
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl jq gosu procps locales \
        "postgresql-${PG_MAJOR}" "postgresql-client-${PG_MAJOR}" \
    && echo 'en_US.UTF-8 UTF-8' >/etc/locale.gen \
    && locale-gen \
    && locale -a | grep -qi '^en_US.utf8$' \
    && rm -rf /var/lib/apt/lists/* \
    && arch="$(dpkg --print-architecture)" \
    && case "${arch}" in \
        amd64) dbmate_sha256="${DBMATE_SHA256_AMD64}" ;; \
        arm64) dbmate_sha256="${DBMATE_SHA256_ARM64}" ;; \
        *) echo "no dbmate checksum for ${arch}" >&2; exit 1 ;; \
    esac \
    && curl -fsSL "https://github.com/amacneil/dbmate/releases/download/v${DBMATE_VERSION}/dbmate-linux-${arch}" \
        -o /usr/local/bin/dbmate \
    && echo "${dbmate_sha256}  /usr/local/bin/dbmate" | sha256sum -c - \
    && chmod +x /usr/local/bin/dbmate

# fake-gcs-server: a local stand-in for Google Cloud Storage, the same emulator
# metamist's own test suite uses. metamist builds a storage.Client() whenever an
# analysis has gs:// outputs; this lets that work locally with no GCP credentials.
# Pinned to a release; each archive's SHA-256 is the digest GitHub publishes for
# that release asset, which matches the release's checksums.txt (checked against
# the downloaded files). Bump a hash together with the version.
ARG FAKE_GCS_VERSION=1.54.0
ARG FAKE_GCS_SHA256_AMD64=3001b2da1fd135eac3a2b38aa2e50c0dbcf89b9359380d4d2b9809c58ca9710b
ARG FAKE_GCS_SHA256_ARM64=b48053046d5e3b569452ad3dfdde141c7059fafb072b73603947c6af59ba524c
RUN arch="$(dpkg --print-architecture)" \
    && case "${arch}" in \
        amd64) fake_gcs_sha256="${FAKE_GCS_SHA256_AMD64}" ;; \
        arm64) fake_gcs_sha256="${FAKE_GCS_SHA256_ARM64}" ;; \
        *) echo "no fake-gcs-server checksum for ${arch}" >&2; exit 1 ;; \
    esac \
    && curl -fsSL "https://github.com/fsouza/fake-gcs-server/releases/download/v${FAKE_GCS_VERSION}/fake-gcs-server_${FAKE_GCS_VERSION}_Linux_${arch}.tar.gz" \
        -o /tmp/fake-gcs.tar.gz \
    && echo "${fake_gcs_sha256}  /tmp/fake-gcs.tar.gz" | sha256sum -c - \
    && tar -xzf /tmp/fake-gcs.tar.gz -C /usr/local/bin fake-gcs-server \
    && rm /tmp/fake-gcs.tar.gz \
    && chmod +x /usr/local/bin/fake-gcs-server

# The temporal_tables extension files (built against the same PG_MAJOR above),
# the python env and the metamist source with its generated client.
COPY --from=temporal_tables_build /out/ /
COPY --from=client_build /opt/venv /opt/venv
COPY --from=client_build /build/ /build/

WORKDIR /build

# Make the unmodified storage.Client() use anonymous creds against the emulator,
# and (local mode) serve hail/hailtop imports as placeholders from an import hook,
# so hail-importing scripts run without the heavy dependency.
COPY sitecustomize.py /opt/venv/lib/python3.14/site-packages/sitecustomize.py

# No seed is baked in: seeding is opt-in and project-specific. With SEED=1 the
# entrypoint runs SEED_SCRIPT once on first boot (mount or copy the script in;
# /opt/seed is the conventional place). Without it the stack boots empty.
COPY docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
COPY mmquery /usr/local/bin/mmquery
RUN chmod +x /usr/local/bin/docker-entrypoint.sh /usr/local/bin/mmquery \
    && mkdir -p /var/log/metamist /app/metamist /etc/metamist

# Local-mode defaults. The default user is granted project-creators by the
# entrypoint, so the same value is used for auth and for the grant. The database
# and login names are fixed by metamist's migrations (they create the
# metamist_server role and grant CONNECT on the metamist database); the
# entrypoint refuses other values.
ENV SM_ENVIRONMENT=local \
    SM_URL=http://localhost:8000 \
    SM_LOCALONLY_DEFAULTUSER=local-user \
    SM_DEV_DB_NAME=metamist \
    SM_DEV_DB_USER=metamist_server \
    SM_DEV_DB_HOST=127.0.0.1 \
    SM_DEV_DB_PORT=5432 \
    METAMIST_DIR=/app/metamist \
    STORAGE_EMULATOR_HOST=http://localhost:4443 \
    GOOGLE_CLOUD_PROJECT=metamist-local \
    SEED=0

WORKDIR /app/metamist
EXPOSE 8000 5432 4443

# Healthy = curl gets /docs from the API on :8000. That happens before the seed step runs, so
# a seeding container turns healthy while the seed is still loading; wait for the
# "[metamist_local_server] ready" log line if you need the seed finished. CI can wait on
# `docker inspect --format '{{.State.Health.Status}}'` / `--wait` /
# `condition: service_healthy` instead of scraping the boot log.
HEALTHCHECK --interval=10s --timeout=5s --start-period=180s --retries=6 \
    CMD curl -fsS http://localhost:8000/docs >/dev/null || exit 1

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
# Keep the container alive and stream the API log when no command is given.
CMD ["tail", "-F", "/var/log/metamist/api.log"]
