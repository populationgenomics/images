"""Local-mode placeholder: `import hailtop.batch as hb` (hb.Batch is subclassed)."""

from _hail_stub import Any as Batch


def __getattr__(_name: str) -> type:
    return Batch
