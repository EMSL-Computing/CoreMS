"""Shared plotting helpers for CoreMS figure lifecycle management."""

from __future__ import annotations

import matplotlib.pyplot as plt


def _finalize_plot(fig, return_fig=False, path=None, **savefig_kwargs):
    """Shared exit path for CoreMS plot methods that support ``return_fig``.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to finalize.
    return_fig : bool, optional
        If True, leave the figure open and return it. The caller owns the
        figure lifecycle (e.g. further customization or ``plt.close(fig)``).
        Default is False.
    path : str or path-like, optional
        If set, save the figure to this path via ``fig.savefig`` before
        showing or returning.
    **savefig_kwargs
        Forwarded to ``fig.savefig`` when ``path`` is set.

    Returns
    -------
    matplotlib.figure.Figure or None
        The open figure if ``return_fig`` is True; otherwise None.

    Notes
    -----
    Behavior matrix:

    - ``return_fig=False``, ``path=None``: call ``plt.show()`` (interactive /
      notebook default).
    - ``return_fig=False``, ``path`` set: save, close the figure, do **not**
      call ``plt.show()`` (batch / headless friendly).
    - ``return_fig=True``: optionally save if ``path`` is set; return the open
      figure without showing or closing it.
    """
    if path is not None:
        fig.savefig(path, **savefig_kwargs)

    if return_fig:
        return fig

    if path is None:
        plt.show()
    else:
        plt.close(fig)
    return None
