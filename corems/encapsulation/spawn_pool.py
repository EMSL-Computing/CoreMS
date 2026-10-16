"""Process pools that still return after pythonnet has started."""

from contextlib import contextmanager
from multiprocessing import get_context
from multiprocessing.pool import TERMINATE


@contextmanager
def spawn_pool(processes):
    """Yield a spawn process pool and kill its workers on exit.

    Parameters
    ----------
    processes : int
        Number of worker processes. The caller chooses this count.

    Yields
    ------
    multiprocessing.pool.Pool
        Pool started with the ``spawn`` start method.

    Notes
    -----
    Linux ``multiprocessing.Pool`` forks. A forked child deadlocks after
    pythonnet has started CoreCLR or Mono, so the pool never finishes.
    Spawn starts a new interpreter. Those workers ignore SIGTERM, and
    ``Pool.terminate`` waits forever in its ``join``. The caller reads
    results inside the ``with`` block. On the way out, including when the
    body raises, the workers are killed and reaped.
    """
    pool = get_context("spawn").Pool(processes)
    try:
        yield pool
    finally:
        _kill_spawn_pool(pool)


def _kill_spawn_pool(pool):
    """Stop a spawn pool without ``Pool.terminate``.

    ``terminate`` signals the workers with SIGTERM and then joins them.
    A worker that has loaded CoreCLR or Mono ignores that signal, so the
    join never returns. The results are already in the parent. Kill the
    workers and cancel the pool finalizer, which would call ``terminate``
    again at garbage collection.
    """
    pool._state = TERMINATE
    for handler in (pool._worker_handler, pool._task_handler, pool._result_handler):
        handler._state = TERMINATE
    pool._change_notifier.put(None)
    for process in list(pool._pool):
        if process.is_alive():
            try:
                process.kill()
            except ProcessLookupError:
                pass
        process.join()
    pool._terminate.cancel()
