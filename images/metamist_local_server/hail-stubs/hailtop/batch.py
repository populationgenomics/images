"""Local-mode placeholder: `import hailtop.batch as hb` (hb.Batch is subclassed)."""

from __future__ import annotations

import _hail_stub

Batch = _hail_stub.HailPlaceholder


def __getattr__(_name: str) -> type:
    return _hail_stub.HailPlaceholder
