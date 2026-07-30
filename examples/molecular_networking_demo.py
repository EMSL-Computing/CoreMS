"""
Molecular Networking Demo Script
=================================

Demonstrates the corems.molecular_networking module end-to-end:
    1. Parse an MSP file and build a FlashEntropy library
    2. Build mock query spectra (with noise) from the library entries
    3. Create MolecularNetwork objects for two search modes:
        - open search + cosine
        - neutral_loss search + cosine
    4. For each network, run a tiered query:
         Stage 1 – Query-vs-Query  (all query pairs)
         Stage 2 – Query-vs-Library (each query vs full library)
         Stage 3 – Library-vs-Library (only library spectra that matched
                   a query above library_similarity_threshold)
    5. Query edges, neighbors, and network statistics per stage
    6. Save outputs for each mode:
        edge list CSV, similarity matrix CSV, SimilarityMatrix NPZ, PNG plots

Run from the repo root:
    python examples/molecular_networking_demo.py
"""

import sys
import numpy as np
from pathlib import Path

# ── Make sure the repo root is on the path ───────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

# ── Output directory ─────────────────────────────────────────────────────────
OUT_DIR = REPO_ROOT / "temp.corems" / "molecular_networking_demo"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────────
# 1. Parse the MSP file and build a FlashEntropy library
# ─────────────────────────────────────────────────────────────────────────────
from corems.molecular_id.search.database_interfaces import MSPInterface
from corems.molecular_networking import MolecularNetwork, SimilarityMatrix

# Public fixture shipped with the repo (no private tmp_data / large MSP required).
# Optional override for a larger local library:
#   export COREMS_NETWORKING_MSP=/path/to/library.msp
import os

_default_msp = REPO_ROOT / "tests/tests_data/lcms/test_db.msp"
_env_msp = os.environ.get("COREMS_NETWORKING_MSP")
MSP_FILE = Path(_env_msp) if _env_msp else _default_msp

# Keep query set small so demo finishes fast and query-vs-query stays readable.
DEMO_QUERY_COUNT = 6

# Threshold for including a library spectrum in Stage 3 (library-vs-library).
LIBRARY_SIMILARITY_THRESHOLD = 0.3

print("=" * 65)
print("STEP 1 – Load MSP library and build FlashEntropy index")
print("=" * 65)

if not MSP_FILE.is_file():
    raise FileNotFoundError(
        f"MSP library not found: {MSP_FILE}\n"
        f"Use the repo fixture or set COREMS_NETWORKING_MSP."
    )

msp = MSPInterface(file_path=str(MSP_FILE))
df = msp._data_frame
print(f"  Parsed {len(df)} spectra from {MSP_FILE.name}")

fe_lib = msp._to_flashentropy(
    input_dataframe=df,
    normalize=True,
    fe_kwargs={
        "normalize_intensity": True,
        "min_ms2_difference_in_da": 0.02,   # must be exactly 2x max_ms2_tolerance_in_da
        "max_ms2_tolerance_in_da": 0.01,
        "max_indexed_mz": 3000,
        "precursor_ions_removal_da": None,
        "noise_threshold": 0,
    },
)
print(f"  FlashEntropy library built ({len(df)} entries indexed)")

# ─────────────────────────────────────────────────────────────────────────────
# 2. Build mock query spectra from the library entries
#    Each spectrum gets small random noise to simulate experimental data.
#    We also record the library index so the engine can do exact lookups.
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("STEP 2 – Build mock query spectra")
print("=" * 65)


class MockSpectrum:
    """Minimal spectrum object compatible with SimilarityEngine."""

    def __init__(self, mz_exp, abundance, name=None):
        self.mz_exp = np.asarray(mz_exp, dtype=float)
        self.abundance = np.asarray(abundance, dtype=float)
        self.name = name


all_spectra = []
all_ids = []
all_precursor_mzs = []

n_lib = len(df)
idxs_of_interest = list(range(min(DEMO_QUERY_COUNT, n_lib)))
for idx in idxs_of_interest:
    row = df.iloc[idx]
    peaks = np.array(row.peaks, dtype=float)
    if len(peaks) == 0:
        continue

    mz_vals = peaks[:, 0]
    abun_vals = peaks[:, 1]

    rng = np.random.default_rng(seed=idx)
    noisy_mz = mz_vals + rng.normal(0, 0.0005, size=mz_vals.shape)  # Reduced from 0.001
    noisy_abun = abun_vals * (1 + rng.normal(0, 0.02, size=abun_vals.shape))  # Reduced from 0.05
    noisy_abun = np.clip(noisy_abun, 0, None)

    name = getattr(row, "refmet_name", None) or f"spectrum_{idx}"
    name = "MockSpec_" + str(name)
    spec_id = getattr(row, "spectra_id", None) or f"spec_{idx:04d}"
    spec_id = "MockSpec_" + str(spec_id)
    pmz = float(getattr(row, "precursormz", 0.0) or 0.0)

    all_spectra.append(MockSpectrum(noisy_mz, noisy_abun, name=name))
    all_ids.append(str(spec_id))
    all_precursor_mzs.append(pmz)

print(f"\n  Created {len(all_spectra)} query spectra from library rows {idxs_of_interest}")
for i in range(min(3, len(all_spectra))):
    print(f"    [{i}] id={all_ids[i]!r}  precursor_mz={all_precursor_mzs[i]:.4f}  "
          f"n_peaks={len(all_spectra[i].mz_exp)}")

print(f"\n  Will run tiered query with {len(all_spectra)} query spectra")
print(f"  Library-vs-library threshold: {LIBRARY_SIMILARITY_THRESHOLD}")

def run_network_demo(search_type: str, label: str):
    print("\n" + "=" * 65)
    print(f"RUN – {label}")
    print("=" * 65)

    network = MolecularNetwork(
        fe_lib=fe_lib,
        search_type=search_type,
        additional_similarities=["cosine"],
        similarity_thresholds={
            "entropy_similarity": 0.6,
            "cosine": 0.6,
        },
        use_parallel=False,
        n_jobs=1,
    )

    run_stage3 = True

    print("\n  Stage 1 + 2: query-vs-query and query-vs-library …")
    network.query_vs_library(
        query_spectra=all_spectra,
        query_ids=all_ids,
        query_precursor_mzs=all_precursor_mzs,
        hydrate_library_similarities=run_stage3,
        library_similarity_threshold=LIBRARY_SIMILARITY_THRESHOLD,
    )
    print(f"\n  After tiered query: {network}")

    print("\n" + "=" * 65)
    print(f"VALIDATION – {label}: Top library matches per query")
    print("=" * 65)

    positive_control_source_idx = 0
    query_id_set = set(all_ids)

    cosine_mat = network.similarity_matrices["cosine"]

    for qid in all_ids:
        neighbors = network.get_spectrum_neighbors(qid, metric="entropy_similarity")
        lib_neighbors = [(nid, score) for nid, score in neighbors if nid.isdigit()]
        if lib_neighbors:
            top_match = lib_neighbors[0]
            status = "✓ MATCH" if top_match[1] >= 0.7 else "~ WEAK"
            print(f"  {qid:25s} → library[{top_match[0]:5s}]: {top_match[1]:.4f} {status}")
        else:
            print(f"  {qid:25s} → NO LIBRARY MATCHES")

    print("\n" + "=" * 65)
    print(f"STATS – {label}")
    print("=" * 65)

    for metric in ["entropy_similarity", "cosine"]:
        mat = network.similarity_matrices[metric]
        all_pairs = mat.get_pairs_above_threshold(0.0)

        qq_pairs = [(a, b, s) for a, b, s in all_pairs if a in query_id_set and b in query_id_set]
        ql_pairs = [(a, b, s) for a, b, s in all_pairs if (a in query_id_set) != (b in query_id_set)]
        ll_pairs = [(a, b, s) for a, b, s in all_pairs if a not in query_id_set and b not in query_id_set]

        threshold = network._threshold_for(metric)
        qq_edges = [(a, b, s) for a, b, s in qq_pairs if s >= threshold]
        ql_edges = [(a, b, s) for a, b, s in ql_pairs if s >= threshold]
        ll_edges = [(a, b, s) for a, b, s in ll_pairs if s >= threshold]

        print(f"\n  [{metric}]  (edge threshold={threshold})")
        print(f"    Stage 1 – Query-vs-Query  : {len(qq_pairs):4d} stored pairs, {len(qq_edges):4d} edges above threshold")
        print(f"    Stage 2 – Query-vs-Library: {len(ql_pairs):4d} stored pairs, {len(ql_edges):4d} edges above threshold")
        print(f"    Stage 3 – Library-vs-Lib  : {len(ll_pairs):4d} stored pairs, {len(ll_edges):4d} edges above threshold")
        print(f"    Total nodes in matrix     : {mat.n_spectra}")

    print("\n" + "=" * 65)
    print(f"SAVE OUTPUTS – {label}")
    print("=" * 65)

    for metric in ["entropy_similarity", "cosine"]:
        csv_path = str(OUT_DIR / f"{search_type}_edges_{metric}.csv")
        network.save_edge_list(csv_path, metric=metric)

    for metric in ["entropy_similarity", "cosine"]:
        mat_path = str(OUT_DIR / f"{search_type}_matrix_{metric}.csv")
        network.save_similarity_matrix(mat_path, metric=metric, threshold=0.0)

    npz_path = str(OUT_DIR / f"{search_type}_entropy_similarity_matrix.npz")
    network.similarity_matrices["entropy_similarity"].save(npz_path)
    reloaded = SimilarityMatrix.load(npz_path)
    print(f"\n  Reloaded matrix: {reloaded}")
    assert reloaded.n_spectra == network.similarity_matrices["entropy_similarity"].n_spectra
    print("  ✓ Save/load round-trip OK")

    # Clustering / HTML plot need optional deps: pip install "corems[networking]"
    try:
        import networkx  # noqa: F401
        import ipysigma  # noqa: F401
        has_viz = True
    except ImportError:
        has_viz = False
        print(
            "  Skipping cluster/plot steps (install optional viz deps with: "
            'pip install "corems[networking]")'
        )

    if has_viz:
        for metric in ["entropy_similarity", "cosine"]:
            cluster_summary = network.compute_network_clusters(
                metric=metric,
                include_queries_only=True,
                max_edges=500,
                cluster_method="weighted_greedy_modularity",
                cluster_super_threshold=400,
                cluster_sparsify_top_k=8,
                cluster_recursive_split=True,
                layout_seed=42,
            )
            print(
                f"  ✓ [{metric}] clusters: {cluster_summary['n_clusters']} clusters "
                f"across {cluster_summary['n_nodes']} nodes"
            )

            cluster_paths = network.save_network_clusters(
                str(OUT_DIR),
                metric=metric,
                run_id=search_type,
            )
            print(f"  ✓ [{metric}] cluster artifacts: {cluster_paths['manifest']}")

            html_path = OUT_DIR / f"{search_type}_network_{metric}.html"
            network.plot_network(
                metric=metric,
                out_path=str(html_path),
                max_edges=500,
                library_label_field=("compound_name", "name", "spectra_id"),
                library_node_attrs=(
                    "compound_name",
                    "name",
                    "spectra_id",
                    "precursor_mz",
                    "precursortype",
                    "inchikey",
                ),
                bypass_clustering=False,
            )
            print(f"  ✓ [{metric}] interactive network HTML: {html_path}")

    return network

print("\n" + "=" * 65)
print("STEP 3+ – Run demos for open and neutral_loss")
print("=" * 65)

networks = {}
networks["open"] = run_network_demo(search_type="open", label="Open Search")

networks["neutral_loss"] = run_network_demo(search_type="neutral_loss", label="Neutral Loss Search")

print("\n" + "=" * 65)
print("DONE – Both networks completed")
print("=" * 65)

print(f"  Built networks: {list(networks.keys())}")
