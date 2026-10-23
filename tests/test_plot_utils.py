"""Unit tests for shared plot finalize helper."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from corems.encapsulation.plot_utils import _finalize_plot


def test_finalize_plot_return_fig_leaves_open(tmp_path):
    fig = plt.figure()
    plt.plot([0, 1], [0, 1])
    out = _finalize_plot(fig, return_fig=True)
    assert out is fig
    assert isinstance(out, Figure)
    assert plt.fignum_exists(fig.number)
    path = tmp_path / "open.png"
    fig.savefig(path)
    assert path.stat().st_size > 0
    plt.close(fig)


def test_finalize_plot_path_without_return_saves_and_closes(tmp_path):
    fig = plt.figure()
    plt.plot([0, 1], [0, 1])
    path = tmp_path / "batch.png"
    out = _finalize_plot(fig, path=path)
    assert out is None
    assert path.is_file() and path.stat().st_size > 0
    assert not plt.fignum_exists(fig.number)


def test_finalize_plot_return_fig_and_path(tmp_path):
    fig = plt.figure()
    plt.plot([0, 1], [0, 1])
    path = tmp_path / "both.png"
    out = _finalize_plot(fig, return_fig=True, path=path)
    assert out is fig
    assert path.is_file() and path.stat().st_size > 0
    assert plt.fignum_exists(fig.number)
    plt.close(fig)


def test_finalize_plot_show_path_returns_none():
    fig = plt.figure()
    plt.plot([0, 1], [0, 1])
    # Agg backend may warn; should not raise
    out = _finalize_plot(fig, return_fig=False, path=None)
    assert out is None
    plt.close("all")
