"""
corems.molecular_networking
===========================

Molecular networking for comparing experimental and library MS2 spectra
using FlashEntropy similarity (and optional cosine) in ``open``,
``identity``, or ``neutral_loss`` modes.

Quick start
-----------

.. code-block:: python

    from corems.molecular_networking import MolecularNetwork

    # Optional: build FE library via MSPInterface / other CoreMS DB helpers
    mn = MolecularNetwork(
        fe_lib=fe_lib,           # or None for query-only
        search_type="open",      # default; try "identity" / "neutral_loss"
        similarity_thresholds={"entropy_similarity": 0.3, "cosine": 0.3},
    )

    # From an LCMSBase object (mass features with best_ms2):
    spectra, ids, pmzs = MolecularNetwork.prepare_query_spectra_from_lcms_object(lcms)

    # All-in-one (stages 1–2; optional stage 3 with hydrate_library_similarities=True)
    mn.query_vs_library(spectra, ids, query_precursor_mzs=pmzs)

    edges = mn.get_network_edges(metric="entropy_similarity")
    mn.save_edge_list("edges.csv")
    # Plotting/clustering require: pip install "corems[networking]"  (networkx + ipysigma)
    mn.plot_network(path="network.png", return_fig=True)  # static
    # mn.plot_interactive_network(out_path="network.html")  # interactive HTML

Stages
------
1. **Query–query** — temporary FE index over experimental spectra.
2. **Query–library** — search against a reference FlashEntropy library.
3. **Library–library** (optional) — all-vs-all among library hits above a threshold.

Internal library node IDs use the ``lib:<index>`` form so they never collide with
numeric mass-feature IDs.  CSV export maps those nodes to FlashEntropy
``spectra_id`` (or ``id``) when available.

Spectrum protocol
-----------------
Query spectra are duck-typed: objects must expose ``.mz_exp`` and ``.abundance``
array-like attributes (e.g. CoreMS MS2 spectra).

Classes
-------
MolecularNetwork
    Main interface for building and querying molecular networks.
SimilarityMatrix
    Sparse symmetric matrix for one similarity metric.
SimilarityEngine
    Computes pairwise similarities (FlashEntropy + optional cosine).
NetworkVisualizeMixin
    Static and interactive plotting (mixed into MolecularNetwork).

See also
--------
examples/notebooks/LCMS_Tutorial.ipynb
    End-to-end single-sample networking section.
examples/notebooks/LCMS_Collection_Tutorial.ipynb
    Multi-sample / consensus representatives.
tests/test_molecular_networking.py
    Unit tests and usage fixtures.
"""

from corems.molecular_networking.similarity_matrix import SimilarityMatrix
from corems.molecular_networking.similarity_engine import SimilarityEngine
from corems.molecular_networking.network_builder import MolecularNetwork
from corems.molecular_networking.network_visualize import NetworkVisualizeMixin

__all__ = [
    "MolecularNetwork",
    "SimilarityMatrix",
    "SimilarityEngine",
    "NetworkVisualizeMixin",
]
