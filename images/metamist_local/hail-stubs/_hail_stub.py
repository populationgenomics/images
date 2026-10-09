"""Universal placeholder for the local-mode hail/hailtop stubs.

`Any` is a class so it can be used as a base class (cpg_utils.hail_batch
subclasses hb.Batch and a hail ServiceBackend at import time), and any attribute
access or call on it returns itself, so eager annotations like hb.batch.job.Job
resolve. It is never executed in local mode (get_batch/init_batch are not called
when SM_ENVIRONMENT=local), so behaviour does not matter, only importability.
"""

from typing import Any as _TypingAny


class _AnyMeta(type):
    def __getattr__(cls, name: str) -> type:
        return Any


class Any(metaclass=_AnyMeta):
    def __getattr__(self, name: str) -> type:
        return Any

    def __call__(self, *_args: _TypingAny, **_kwargs: _TypingAny) -> type:  # noqa: ANN401
        return Any
