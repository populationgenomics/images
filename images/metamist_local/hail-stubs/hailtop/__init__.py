"""Local-mode placeholder for `hailtop` (see _hail_stub)."""

from _hail_stub import Any as _Any


def __getattr__(_name: str) -> type:
    return _Any
