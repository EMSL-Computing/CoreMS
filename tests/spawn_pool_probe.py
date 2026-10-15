"""Return from a spawn pool after pythonnet has started in the parent.

Run this file as a script. A spawn child imports it again, so the pool is
started only from ``__main__``.
"""

import sys
from pathlib import Path

# This file lives in tests/. A script run puts that directory on sys.path,
# and the editable install can point at a different checkout. Import the
# CoreMS next to this file.
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))


def _load_runtime():
    from pythonnet import load

    try:
        load()
    except RuntimeError:
        # This interpreter already loaded a runtime.
        pass
    import clr  # noqa: F401


def touch_runtime(value):
    """Load the .NET runtime and return ``value`` unchanged."""
    _load_runtime()
    return value


def main():
    _load_runtime()
    from corems.encapsulation.spawn_pool import spawn_pool

    with spawn_pool(2) as pool:
        result = pool.map(touch_runtime, [1, 2])
    if result != [1, 2]:
        print(result)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
