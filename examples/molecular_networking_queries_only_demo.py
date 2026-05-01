"""
Queries-Only Molecular Networking Demo (no reference FE library)
================================================================

Demonstrates using MolecularNetwork with experimental spectra only:

- No pre-built FlashEntropy library (fe_lib=None).
- Stage 1: query-vs-query similarities only.
- No query-vs-library or library-vs-library stages.

Run from repo root:
    python examples/molecular_networking_queries_only_demo.py
"""

import sys
from pathlib import Path

import numpy as np

# ── Ensure repo root on path ─────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from corems.molecular_networking import MolecularNetwork  # noqa: E402


class MockSpectrum:
    """Minimal spectrum object compatible with SimilarityEngine.

    Attributes
    ----------
    mz_exp : np.ndarray
        Fragment m/z values.
    abundance : np.ndarray
        Corresponding intensities.
    name : str | None
        Optional human-readable label.
    """

    def __init__(self, mz_exp, abundance, name=None):
        self.mz_exp = np.asarray(mz_exp, dtype=float)
        self.abundance = np.asarray(abundance, dtype=float)
        self.name = name


# Base spectrum: 5 fragment peaks well below precursor
_BASE_MZ = np.array([100.0, 150.0, 200.0, 250.0, 280.0], dtype=float)
_BASE_ABUN = np.array([1.0, 0.8, 0.6, 0.4, 0.2], dtype=float)
_BASE_PMZ = 300.0
_SHIFT = 3  # uniform shift for EXP_4


def build_experimental_spectra():
    """Build controlled mock spectra with known expected similarity outcomes.

    Spectrum design
    ---------------
    EXP_0 : base spectrum (pmz=300, fragments=[100,150,200,250,280])
    EXP_1 : EXP_0 + 0.0001 Da offset on all mz (pmz same)
             → open: entropy≈1, cosine≈1 | neutral_loss: entropy≈1, cosine≈1
    EXP_2 : same mz as EXP_0, different abundance (pmz same)
             → open: high but <1 | neutral_loss: high but <1
    EXP_3 : same mz/abun as EXP_0, different precursor (pmz=400)
             → open: entropy=1, cosine=1 | neutral_loss: entropy<1, cosine<1
               (neutral losses differ because precursor differs)
    EXP_4 : EXP_0 shifted +50 Da on all mz AND precursor
             → open: entropy≈0, cosine≈0 (no peak overlap in fragment space)
             → neutral_loss: entropy≈1, cosine≈1 (same neutral losses)

    Returns
    -------
    spectra : list[MockSpectrum]
    ids : list[str]
    precursor_mzs : list[float]
    """
    spectra: list[MockSpectrum] = []
    ids: list[str] = []
    precursor_mzs: list[float] = []

    # EXP_0 – base
    spectra.append(MockSpectrum(_BASE_MZ.copy(), _BASE_ABUN.copy(), name="exp_0"))
    ids.append("EXP_0")
    precursor_mzs.append(_BASE_PMZ)

    # EXP_1 – tiny offset (+0.0001 Da) on all fragment mz, same precursor
    spectra.append(MockSpectrum(_BASE_MZ + 0.0001, _BASE_ABUN.copy(), name="exp_1"))
    ids.append("EXP_1")
    precursor_mzs.append(_BASE_PMZ)

    # EXP_2 – same mz as EXP_0, different abundance
    spectra.append(MockSpectrum(_BASE_MZ.copy(), _BASE_ABUN * np.array([0.5, 1.2, 0.9, 1.5, 0.7]), name="exp_2"))
    ids.append("EXP_2")
    precursor_mzs.append(_BASE_PMZ)

    # EXP_3 – same mz/abun as EXP_0, different precursor (400 Da)
    spectra.append(MockSpectrum(_BASE_MZ.copy(), _BASE_ABUN.copy(), name="exp_3"))
    ids.append("EXP_3")
    precursor_mzs.append(301.0)

    # EXP_4 – EXP_0 shifted +50 Da on all fragments AND precursor
    spectra.append(MockSpectrum(_BASE_MZ + _SHIFT, _BASE_ABUN.copy(), name="exp_4"))
    ids.append("EXP_4")
    precursor_mzs.append(_BASE_PMZ + _SHIFT)

    return spectra, ids, precursor_mzs


# Expected similarity outcomes per search type
# (id1, id2, open_entropy, open_cosine, nl_entropy, nl_cosine, note)
_EXPECTED = [
    ("EXP_0", "EXP_1", "≈1.0", "≈1.0", "≈1.0", "≈1.0", "tiny offset → same NL"),
    ("EXP_0", "EXP_2", "high", "high", "high", "high", "same mz, diff abun"),
    ("EXP_0", "EXP_3", "≈1.0", "≈1.0", "<1",   "<1",   "same frags, diff pmz → diff NL"),
    ("EXP_0", "EXP_4", "≈0",   "≈0",   "≈1.0", "≈1.0", "+50 Da shift → same NL"),
]


def run_search(search_type: str, q_spectra, q_ids, q_precursor_mzs):
    """Run query-vs-query for a given search_type and print results."""
    for metric in ["entropy_similarity", "cosine"]:
        print(f"\n{'-' * 60}")
        print(f"Search type: {search_type.upper()} | Similarity metric: {metric}")
        print(f"{'-' * 60}")


        network = MolecularNetwork(
            fe_lib=None,
            search_type=search_type,
            additional_similarities=["cosine"],
            similarity_thresholds={
                "entropy_similarity": 0.5,
                "cosine": 0.5,
            },
            use_parallel=False,
            n_jobs=1,
        )

        network.run_query_vs_query_only(
            query_spectra=q_spectra,
            query_ids=q_ids,
            query_precursor_mzs=q_precursor_mzs,
        )

        print(f"  {network}")

        for qid in ['EXP_0']:
            neighbors_entropy = network.get_spectrum_neighbors(qid, metric="entropy_similarity")
            neighbors_cosine = network.get_spectrum_neighbors(qid, metric="cosine")
            print(f"\n  Neighbors of {qid} (entropy_similarity):")
            for nid, score in neighbors_entropy:
                print(f"    {qid} ↔ {nid}: {score:.4f}")
            print(f"  Neighbors of {qid} (cosine):")
            for nid, score in neighbors_cosine:
                print(f"    {qid} ↔ {nid}: {score:.4f}")

        for metric in ["entropy_similarity", "cosine"]:
            stats = network.get_network_stats(metric=metric)
            print(f"  [{metric}] n_edges={stats['n_edges']}, density={stats['density']:.2f}")



def main():
    print("=" * 65)
    print("QUERIES-ONLY DEMO – No FE library")
    print("=" * 65)

    q_spectra, q_ids, q_precursor_mzs = build_experimental_spectra()
    print(f"  Built {len(q_spectra)} experimental spectra")
    print(f"  EXP_4 is exact +{_SHIFT:.0f} Da shift of EXP_0 (neutral-loss test)")
    print(f"  EXP_3 is same fragments as EXP_0 but different precursor (303 Da)")

    for search_type in ("open", "neutral_loss"):
        run_search(search_type, q_spectra, q_ids, q_precursor_mzs)


if __name__ == "__main__":
    main()
