"""Universal placeholder for the local-mode hail/hailtop stubs.

`HailPlaceholder` is a class so it can be used as a base class (cpg_utils.hail_batch
subclasses hb.Batch and a hail ServiceBackend at import time), and any attribute
access or call on it returns the class itself, so eager annotations like
hb.batch.job.Job resolve. It is never executed in local mode (get_batch/init_batch
are not called when SM_ENVIRONMENT=local), so behaviour does not matter, only
importability.
"""


class _HailPlaceholderMeta(type):
    def __getattr__(cls, name: str) -> type:
        return HailPlaceholder


class HailPlaceholder(metaclass=_HailPlaceholderMeta):
    def __getattr__(self, name: str) -> type:
        return HailPlaceholder

    def __call__(self, *_args: object, **_kwargs: object) -> type:
        return HailPlaceholder
