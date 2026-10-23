"""Spawn pools return their results and do not leave workers behind."""

import subprocess
import sys
from pathlib import Path

import pytest

from corems.encapsulation.spawn_pool import spawn_pool


def _square(value):
    return value * value


def _boom(_value):
    raise RuntimeError("worker failed")


def _echo(value):
    return value


def _workers(pool):
    return list(pool._pool)


def test_map_returns_results_and_reaps_workers():
    with spawn_pool(2) as pool:
        workers = _workers(pool)
        result = pool.map(_square, [1, 2, 3])
    assert result == [1, 4, 9]
    assert workers
    assert all(not worker.is_alive() for worker in workers)


def test_worker_exception_reaps_workers():
    workers = []
    with pytest.raises(RuntimeError, match="worker failed"):
        with spawn_pool(2) as pool:
            workers.extend(_workers(pool))
            pool.map(_boom, [1])
    assert workers
    assert all(not worker.is_alive() for worker in workers)


def test_imap_unordered_returns_every_item():
    with spawn_pool(2) as pool:
        workers = _workers(pool)
        got = list(pool.imap_unordered(_echo, [1, 2, 3, 4]))
    assert sorted(got) == [1, 2, 3, 4]
    assert all(not worker.is_alive() for worker in workers)


def test_one_worker_reaps():
    with spawn_pool(1) as pool:
        workers = _workers(pool)
        result = pool.map(_square, [3])
    assert result == [9]
    assert all(not worker.is_alive() for worker in workers)


def test_empty_map_reaps():
    with spawn_pool(2) as pool:
        workers = _workers(pool)
        result = pool.map(_square, [])
    assert result == []
    assert all(not worker.is_alive() for worker in workers)


def test_spawn_pool_returns_after_pythonnet():
    check = subprocess.run(
        [sys.executable, "-c", "from pythonnet import load; load()"],
        capture_output=True,
        text=True,
        timeout=60,
    )
    if check.returncode != 0:
        message = check.stderr.strip() or check.stdout.strip()
        pytest.skip(message or "pythonnet runtime is not available")
    probe = Path(__file__).resolve().with_name("spawn_pool_probe.py")
    subprocess.run([sys.executable, str(probe)], check=True, timeout=60)
