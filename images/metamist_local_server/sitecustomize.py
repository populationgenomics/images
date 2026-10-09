"""Local mode only: point google-cloud-storage at the fake GCS emulator with
anonymous credentials, and make `hail` / `hailtop` importable as placeholders.

metamist builds an unmodified ``storage.Client()`` (no args) whenever an
analysis has outputs. Off-GCP that calls ``google.auth.default()`` and fails.
When ``STORAGE_EMULATOR_HOST`` is set we instead default the client to
``AnonymousCredentials`` and a dummy project, so the client talks to the local
fake-gcs-server with no real GCP credentials. Inert when the env var is unset,
and it never overrides credentials/project a caller passes explicitly (so it
does not change metamist's own test client).

Patching ``storage.Client.__init__`` is deliberate: the client is constructed
inside metamist code that this image runs unmodified, so there is no call site
to pass credentials at.
"""

from __future__ import annotations

import os
import sys

_HAIL_STUBS = '/opt/hail-stubs'


def _add_hail_stubs() -> None:
    """Make `import hail` / `import hailtop` resolve to the local placeholders.

    Scripts like create_test_subset.py import cpg_utils.hail_batch but never call
    into hail when SM_ENVIRONMENT=local. The stubs are appended (lowest priority)
    so a real hail install, if ever added, takes precedence.
    """
    if os.path.isdir(_HAIL_STUBS) and _HAIL_STUBS not in sys.path:
        sys.path.append(_HAIL_STUBS)


def _patch_storage_client() -> None:
    """Default storage.Client() to anonymous credentials and a dummy project.

    Does nothing if google-cloud-storage isn't installed in this interpreter.
    """
    try:
        from google.auth import credentials
        from google.cloud import storage
    except ImportError:
        return

    original_init = storage.Client.__init__

    def local_init(self: storage.Client, *args: object, **kwargs: object) -> None:
        kwargs.setdefault('credentials', credentials.AnonymousCredentials())
        kwargs.setdefault(
            'project', os.environ.get('GOOGLE_CLOUD_PROJECT', 'metamist-local')
        )
        original_init(self, *args, **kwargs)

    storage.Client.__init__ = local_init


if os.environ.get('SM_ENVIRONMENT', '').lower() == 'local':
    _add_hail_stubs()

if os.environ.get('STORAGE_EMULATOR_HOST'):
    _patch_storage_client()
