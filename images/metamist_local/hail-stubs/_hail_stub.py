"""Universal placeholder for the local-mode hail/hailtop stubs.

The stubs exist only so metamist scripts that import cpg_utils.hail_batch can be
imported without hail installed. They are never executed in local mode
(get_batch/init_batch are not called when SM_ENVIRONMENT=local), so behaviour
does not matter, only importability.

The metaclass and ``__getattr__`` are deliberate: cpg_utils subclasses hb.Batch
and a hail ServiceBackend at import time and evaluates annotations such as
hb.batch.job.Job, so any attribute chain must resolve. Hand-written stubs naming
each attribute would break whenever cpg_utils starts touching a new one.
"""

from __future__ import annotations


class _HailPlaceholderMeta(type):
    def __getattr__(cls, name: str) -> type:
        return HailPlaceholder


class HailPlaceholder(metaclass=_HailPlaceholderMeta):
    """Stand-in for any hail object.

    Every attribute access or call, on the class or an instance, returns
    HailPlaceholder, and it can be subclassed.
    """

    def __getattr__(self, name: str) -> type:
        return HailPlaceholder

    def __call__(self, *_args: object, **_kwargs: object) -> type:
        return HailPlaceholder
