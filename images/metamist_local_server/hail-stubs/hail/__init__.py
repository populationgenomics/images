"""Local-mode placeholder for `hail` (see _hail_stub)."""

from __future__ import annotations

import _hail_stub


def __getattr__(_name: str) -> type:
    return _hail_stub.HailPlaceholder
