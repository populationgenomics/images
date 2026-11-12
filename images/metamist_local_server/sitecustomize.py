"""Local mode only: point google-cloud-storage at the fake GCS emulator with
anonymous credentials, and make `hail` / `hailtop` importable as placeholders.

metamist builds an unmodified ``storage.Client()`` (no args) whenever an
analysis has outputs. Off-GCP that calls ``google.auth.default()`` and fails.
When ``STORAGE_EMULATOR_HOST`` is set we instead default the client to
``AnonymousCredentials`` and a dummy project, so the client talks to the local
fake-gcs-server with no real GCP credentials. Inert when the env var is unset,
and it never overrides credentials/project a caller passes explicitly, by
keyword or positionally (so it does not change metamist's own test client).

Patching ``storage.Client.__init__`` is deliberate: the client is constructed
inside metamist code that this image runs unmodified, so there is no call site
to pass credentials at.

Scripts such as create_test_subset.py import cpg_utils.hail_batch, which imports
hail and hailtop at module level, but never run hail code when
SM_ENVIRONMENT=local. A meta path finder serves those packages, and any
submodule of them, as placeholders so the imports succeed without the heavy
dependency. Anything that would actually run hail raises a RuntimeError.
"""

from __future__ import annotations

import importlib.abc
import importlib.machinery
import os
import sys
import types
from collections.abc import Sequence

_HAIL_ROOTS = frozenset({'hail', 'hailtop'})
_HAIL_UNAVAILABLE = 'hail is not available in metamist_local_server (local mode)'


def _is_dunder(name: str) -> bool:
    return name.startswith('__') and name.endswith('__')


class _HailPlaceholderMeta(type):
    """Class attribute access on a placeholder returns the placeholder class."""

    def __getattr__(cls, name: str) -> type[HailPlaceholder]:
        if _is_dunder(name):
            raise AttributeError(name)
        return HailPlaceholder


class HailPlaceholder(metaclass=_HailPlaceholderMeta):
    """Stand-in for any hail object, enough to import and subclass.

    Attribute access on the class or an instance returns the placeholder class
    (dunder names raise AttributeError, so protocols such as copy and pickle
    see nothing), and it can be subclassed: cpg_utils subclasses hb.Batch and
    a ServiceBackend at import time, and evaluates annotations such as
    hb.batch.job.Job, so any attribute chain must resolve.

    Using one raises RuntimeError: instantiating it, calling an instance, or
    calling a subclass that does not override ``__init__``.
    """

    def __init__(self, *_args: object, **_kwargs: object) -> None:
        raise RuntimeError(_HAIL_UNAVAILABLE)

    def __getattr__(self, name: str) -> type[HailPlaceholder]:
        if _is_dunder(name):
            raise AttributeError(name)
        return HailPlaceholder

    def __call__(self, *_args: object, **_kwargs: object) -> None:
        raise RuntimeError(_HAIL_UNAVAILABLE)


class _HailPlaceholderModule(types.ModuleType):
    """A package whose every public attribute is the placeholder class."""

    def __getattr__(self, name: str) -> type[HailPlaceholder]:
        if _is_dunder(name):
            raise AttributeError(f'module {self.__name__!r} has no attribute {name!r}')
        return HailPlaceholder


class _HailPlaceholderFinder(importlib.abc.MetaPathFinder, importlib.abc.Loader):
    """Serve hail, hailtop and their dotted submodules as placeholder packages.

    Appended to sys.meta_path, so the normal finders run first and a real hail
    install wins. Submodules are served only under a placeholder parent: with a
    real hail installed, a missing submodule stays an ImportError.
    """

    def find_spec(
        self,
        fullname: str,
        path: Sequence[str] | None = None,  # noqa: ARG002
        target: types.ModuleType | None = None,  # noqa: ARG002
    ) -> importlib.machinery.ModuleSpec | None:
        if fullname.partition('.')[0] not in _HAIL_ROOTS:
            return None
        parent = fullname.rpartition('.')[0]
        if parent and not isinstance(sys.modules.get(parent), _HailPlaceholderModule):
            return None
        # is_package makes __path__ = [], so submodule imports come back here.
        return importlib.machinery.ModuleSpec(fullname, self, is_package=True)

    def create_module(
        self, spec: importlib.machinery.ModuleSpec
    ) -> types.ModuleType | None:
        return _HailPlaceholderModule(spec.name)

    def exec_module(self, module: types.ModuleType) -> None:
        """Nothing to run: the placeholder module has no code."""


def _add_hail_placeholders() -> None:
    """Make `import hail` / `import hailtop` resolve to the local placeholders."""
    if not any(isinstance(f, _HailPlaceholderFinder) for f in sys.meta_path):
        sys.meta_path.append(_HailPlaceholderFinder())


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
        # Client(project, credentials, ...): a positional argument already
        # supplies the parameter, and a keyword for it would collide.
        if len(args) < 1:
            kwargs.setdefault(
                'project', os.environ.get('GOOGLE_CLOUD_PROJECT', 'metamist-local')
            )
        if len(args) < 2:  # noqa: PLR2004
            kwargs.setdefault('credentials', credentials.AnonymousCredentials())
        original_init(self, *args, **kwargs)

    storage.Client.__init__ = local_init


if os.environ.get('SM_ENVIRONMENT', '').lower() == 'local':
    _add_hail_placeholders()

if os.environ.get('STORAGE_EMULATOR_HOST'):
    _patch_storage_client()
