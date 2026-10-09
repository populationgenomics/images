"""Local-mode placeholder for `hail` (see _hail_stub)."""

from _hail_stub import HailPlaceholder


def __getattr__(_name: str) -> type:
    return HailPlaceholder
