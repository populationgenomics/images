"""Local mode only: point google-cloud-storage at the fake GCS emulator with
anonymous credentials, and make `hail` / `hailtop` importable as placeholders.

metamist builds an unmodified ``storage.Client()`` (no args) whenever an
analysis has outputs. Off-GCP that calls ``google.auth.default()`` and fails.
When ``STORAGE_EMULATOR_HOST`` is set we instead default the client to
``AnonymousCredentials`` and a dummy project, so the client talks to the local
fake-gcs-server with no real GCP credentials. Inert when the env var is unset,
and it never overrides credentials/project a caller passes explicitly (so it
does not change metamist's own test client).
"""

import os
import sys

# Local-mode only: make `import hail` / `import hailtop` resolve to minimal
# placeholders so scripts like create_test_subset.py (which import
# cpg_utils.hail_batch but never call into hail when SM_ENVIRONMENT=local) can
# run without the heavy hail dependency. Appended (lowest priority) so a real
# hail install, if ever added, takes precedence.
if os.environ.get('SM_ENVIRONMENT', '').lower() == 'local':
    _hail_stubs = '/opt/hail-stubs'
    if os.path.isdir(_hail_stubs) and _hail_stubs not in sys.path:
        sys.path.append(_hail_stubs)

if os.environ.get('STORAGE_EMULATOR_HOST'):
    try:
        from google.auth.credentials import AnonymousCredentials
        from google.cloud import storage
    except ImportError:
        # google-cloud-storage isn't installed in this interpreter: nothing to patch.
        pass
    else:
        _orig_init = storage.Client.__init__

        def _local_init(self: storage.Client, *args: object, **kwargs: object) -> None:
            kwargs.setdefault('credentials', AnonymousCredentials())
            kwargs.setdefault(
                'project', os.environ.get('GOOGLE_CLOUD_PROJECT', 'metamist-local')
            )
            _orig_init(self, *args, **kwargs)

        storage.Client.__init__ = _local_init
