"""
Molecular Networking Demo Script
=================================

Demonstrates the corems.molecular_networking module end-to-end:
    1. Parse an MSP file and build a FlashEntropy library
    2. Build mock query spectra (with noise) from the library entries
    3. Create a MolecularNetwork (open search + cosine)
    4. Run a tiered query:
         Stage 1 – Query-vs-Query  (all query pairs)
         Stage 2 – Query-vs-Library (each query vs full library)
         Stage 3 – Library-vs-Library (only library spectra that matched
                   a query above library_similarity_threshold)
    5. Query edges, neighbors, and network statistics per stage
    6. Save outputs: edge list CSV, similarity matrix CSV, GraphML, PNG plots

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

MSP_FILE = REPO_ROOT / "tests/tests_data/lcms/test_db.msp"

# Larger lipid library (51 k spectra) for library search.
LARGE_MSP = REPO_ROOT / "tmp_data" / "20250407_database.msp"

# Keep query set small so demo finishes fast and query-vs-query stays readable.
DEMO_QUERY_COUNT = 6

# Threshold for including a library spectrum in Stage 3 (library-vs-library).
LIBRARY_SIMILARITY_THRESHOLD = 0.3

print("=" * 65)
print("STEP 1 – Load MSP library and build FlashEntropy index")
print("=" * 65)

msp = MSPInterface(file_path=str(LARGE_MSP))
df = msp._data_frame
print(f"  Parsed {len(df)} spectra from {LARGE_MSP.name}")

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

for lib_idx, row in enumerate(df.itertuples(index=False)):
    peaks = np.array(row.peaks, dtype=float)
    if len(peaks) == 0:
        continue

    mz_vals = peaks[:, 0]
    abun_vals = peaks[:, 1]

    rng = np.random.default_rng(seed=lib_idx)
    noisy_mz = mz_vals + rng.normal(0, 0.0005, size=mz_vals.shape)  # Reduced from 0.001
    noisy_abun = abun_vals * (1 + rng.normal(0, 0.02, size=abun_vals.shape))  # Reduced from 0.05
    noisy_abun = np.clip(noisy_abun, 0, None)

    name = getattr(row, "compound_name", None) or f"spectrum_{lib_idx}"
    spec_id = getattr(row, "spectra_id", None) or f"spec_{lib_idx:04d}"
    pmz = float(getattr(row, "precursormz", 0.0) or 0.0)

    all_spectra.append(MockSpectrum(noisy_mz, noisy_abun, name=name))
    all_ids.append(str(spec_id))
    all_precursor_mzs.append(pmz)

    if len(all_spectra) >= DEMO_QUERY_COUNT:
        break

# Add positive control: exact library match with different ID
if len(df) > 0:
    control_idx = 0  # Use first library entry
    control_row = df.iloc[control_idx]
    control_peaks = np.array(control_row.peaks, dtype=float)
    
    if len(control_peaks) > 0:
        control_mz = control_peaks[:, 0]
        control_abun = control_peaks[:, 1]
        control_name = f"{control_row.compound_name}_EXACT_COPY"
        control_id = "POSITIVE_CONTROL"
        control_pmz = float(control_row.precursormz or 0.0)
        control_spectra_id = control_row.spectra_id if hasattr(control_row, 'spectra_id') else f"spec_{control_idx:04d}"
        
        all_spectra.append(MockSpectrum(control_mz, control_abun, name=control_name))
        all_ids.append(control_id)
        all_precursor_mzs.append(control_pmz)
        
        print(f"\n  Added positive control: {control_id}")
        print(f"    Exact copy of dataframe row[{control_idx}]:")
        print(f"      - compound_name: {control_row.compound_name}")
        print(f"      - spectra_id: {control_spectra_id}")
        print(f"      - precursor_mz: {control_pmz:.4f}")
        print(f"      - n_peaks: {len(control_peaks)}")
        print(f"    Note: FE library indices differ from dataframe indices (sorted by precursor m/z)")

print(f"\n  Created {len(all_spectra)} query spectra ({len(all_spectra) - 1} noisy + 1 exact)")
for i in range(min(3, len(all_spectra))):
    print(f"    [{i}] id={all_ids[i]!r}  precursor_mz={all_precursor_mzs[i]:.4f}  "
          f"n_peaks={len(all_spectra[i].mz_exp)}")

print(f"\n  Will run tiered query with {len(all_spectra)} query spectra")
print(f"  Library-vs-library threshold: {LIBRARY_SIMILARITY_THRESHOLD}")

# ─────────────────────────────────────────────────────────────────────────────
# 3. Create the MolecularNetwork (open search with cosine)
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("STEP 3 – Create MolecularNetwork (open search)")
print("=" * 65)

# Convert dataframe to list of dicts for library_spectra parameter
library_spectra = df.to_dict(orient="records")

network = MolecularNetwork(
    fe_lib=fe_lib,
    search_type="open",             # open search
    additional_similarities=["cosine"],
    similarity_thresholds={
        "entropy_similarity": 0.3,      # lower threshold for demo data
        "cosine": 0.3,
    },
    peak_sep_da=0.02,
    ms1_tolerance_da=0.05,             # generous tolerance for demo
    ms2_tolerance_da=0.01,             # must be <= peak_sep_da / 2
    entropy_threshold_low=0.05,        # trigger cosine for any non-trivial match
    use_parallel=False,                # keep demo single-threaded
    n_jobs=1,
    library_spectra=library_spectra,   # Pass library spectra for cosine computation
)

# ─────────────────────────────────────────────────────────────────────────────
# 4. Tiered query
#    Stage 1: Query-vs-Query
#    Stage 2: Query-vs-Library
#    Stage 3: Library-vs-Library (filtered by library_similarity_threshold)
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("STEP 4 – Tiered query (3 stages)")
print("=" * 65)

print("\n  Stage 1 + 2: query-vs-query and query-vs-library …")
network.query_vs_library(
    query_spectra=all_spectra,
    query_ids=all_ids,
    query_precursor_mzs=all_precursor_mzs,
    hydrate_library_similarities=True,
    library_similarity_threshold=LIBRARY_SIMILARITY_THRESHOLD,
)
print(f"\n  After tiered query: {network}")

# ─────────────────────────────────────────────────────────────────────────────
# 4.5. Validation: Check query-to-library matches
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("VALIDATION – Top library matches per query")
print("=" * 65)

# Record which library spectrum POSITIVE_CONTROL is copied from
positive_control_source_idx = 0  # Set when creating POSITIVE_CONTROL above

entropy_mat = network.similarity_matrices["entropy_similarity"]
cosine_mat = network.similarity_matrices["cosine"]

for qid in all_ids:
    neighbors = network.get_spectrum_neighbors(qid, metric="entropy_similarity")
    # Filter to library IDs only (numeric strings)
    lib_neighbors = [(nid, score) for nid, score in neighbors if nid.isdigit()]
    if lib_neighbors:
        top_match = lib_neighbors[0]
        status = "✓ MATCH" if top_match[1] >= 0.7 else "~ WEAK"
        print(f"  {qid:25s} → library[{top_match[0]:5s}]: {top_match[1]:.4f} {status}")
    else:
        print(f"  {qid:25s} → NO LIBRARY MATCHES")

# Special validation for POSITIVE_CONTROL
print("\n" + "=" * 65)
print("POSITIVE_CONTROL Validation")
print("=" * 65)
print(f"  POSITIVE_CONTROL is an exact copy of dataframe row[{positive_control_source_idx}]")
print(f"  Note: FE library indices may differ from dataframe indices (sorted by precursor m/z)")

# Get top entropy match
pc_entropy_neighbors = network.get_spectrum_neighbors("POSITIVE_CONTROL", metric="entropy_similarity")
pc_lib_entropy = [(nid, score) for nid, score in pc_entropy_neighbors if nid.isdigit()]
if pc_lib_entropy:
    top_entropy_lib_id, top_entropy_score = pc_lib_entropy[0]
    
    # Try to get library spectrum info
    try:
        lib_spec = fe_lib[int(top_entropy_lib_id)]
        lib_pmz = lib_spec.get("precursor_mz", "unknown")
        lib_peaks_count = len(lib_spec.get("peaks", []))
        print(f"\n  Top entropy match: library[{top_entropy_lib_id}]")
        print(f"    - precursor_mz: {lib_pmz:.4f}")
        print(f"    - n_peaks (cleaned): {lib_peaks_count}")
        print(f"    - entropy similarity: {top_entropy_score:.6f}")
    except:
        print(f"\n  Top entropy match: library[{top_entropy_lib_id}] = {top_entropy_score:.6f}")
    
    if top_entropy_score >= 0.95:
        print(f"    ✓ EXCELLENT (≥0.95)")
    elif top_entropy_score >= 0.80:
        print(f"    ✓ GOOD (≥0.80)")
    else:
        print(f"    ~ MODERATE (<0.80) - Expected ~1.0 for exact copy!")
    
    # Check cosine for the same library spectrum
    cosine_score = cosine_mat.get_similarity("POSITIVE_CONTROL", top_entropy_lib_id)
    print(f"\n  Cosine similarity to library[{top_entropy_lib_id}]: {cosine_score:.6f}")
    if cosine_score >= 0.95:
        print(f"    ✓ EXCELLENT (≥0.95)")
    elif cosine_score >= 0.80:
        print(f"    ✓ GOOD (≥0.80)")
    else:
        print(f"    ~ MODERATE (<0.80) - Expected ~1.0 for exact copy!")
    
    # Check consistency between metrics
    if abs(top_entropy_score - cosine_score) < 0.1:
        print(f"\n  ✓ Metrics are consistent (diff = {abs(top_entropy_score - cosine_score):.4f})")
    else:
        print(f"\n  ⚠ Metrics differ significantly (diff = {abs(top_entropy_score - cosine_score):.4f})")
        print(f"     This suggests the metrics are using different cleaned spectra!")
else:
    print(f"\n  ✗ No library matches found!")

# ─────────────────────────────────────────────────────────────────────────────
# 5. Inspect per-stage results (renumbered from original step 5)
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("STEP 5 – Per-stage network statistics")
print("=" * 65)

query_id_set = set(all_ids)

for metric in ["entropy_similarity", "cosine"]:
    mat = network.similarity_matrices[metric]
    all_pairs = mat.get_pairs_above_threshold(0.0)

    # Classify each pair by stage
    qq_pairs = [(a, b, s) for a, b, s in all_pairs
                if a in query_id_set and b in query_id_set]
    ql_pairs = [(a, b, s) for a, b, s in all_pairs
                if (a in query_id_set) != (b in query_id_set)]
    ll_pairs = [(a, b, s) for a, b, s in all_pairs
                if a not in query_id_set and b not in query_id_set]

    threshold = network._threshold_for(metric)
    qq_edges = [(a, b, s) for a, b, s in qq_pairs if s >= threshold]
    ql_edges = [(a, b, s) for a, b, s in ql_pairs if s >= threshold]
    ll_edges = [(a, b, s) for a, b, s in ll_pairs if s >= threshold]

    print(f"\n  [{metric}]  (edge threshold={threshold})")
    print(f"    Stage 1 – Query-vs-Query  : {len(qq_pairs):4d} stored pairs, "
          f"{len(qq_edges):4d} edges above threshold")
    print(f"    Stage 2 – Query-vs-Library: {len(ql_pairs):4d} stored pairs, "
          f"{len(ql_edges):4d} edges above threshold")
    print(f"    Stage 3 – Library-vs-Lib  : {len(ll_pairs):4d} stored pairs, "
          f"{len(ll_edges):4d} edges above threshold")
    print(f"    Total nodes in matrix     : {mat.n_spectra}")

# ─────────────────────────────────────────────────────────────────────────────
# 6. Query the network (edges, neighbors, stats)
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("STEP 6 – Query the network")
print("=" * 65)

for metric in ["entropy_similarity", "cosine"]:
    edges = network.get_network_edges(metric=metric)
    stats = network.get_network_stats(metric=metric)
    print(f"\n  [{metric}]")
    print(f"    Stats: {stats}")
    if edges:
        print(f"    Edges ({len(edges)} total, showing up to 5):")
        for id1, id2, score in edges[:5]:
            print(f"      {id1!r} ↔ {id2!r}  score={score:.4f}")
    else:
        print("    No edges above threshold.")

# Neighbors for the first spectrum
first_id = all_ids[0]
print(f"\n  Neighbors of {first_id!r} (entropy_similarity):")
neighbors = network.get_spectrum_neighbors(first_id, metric="entropy_similarity")
if neighbors:
    for nid, score in neighbors:
        node_type = "query" if nid in query_id_set else "library"
        print(f"    {nid!r} [{node_type}]  score={score:.4f}")
else:
    print("    None above threshold.")

# ─────────────────────────────────────────────────────────────────────────────
# 7. Inspect the SimilarityMatrix objects directly
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("STEP 7 – Inspect SimilarityMatrix objects")
print("=" * 65)

for metric, mat in network.similarity_matrices.items():
    print(f"\n  {mat}")
    df_pairs = mat.to_dataframe(threshold=0.0)
    print(f"    All stored pairs: {len(df_pairs)}")
    if not df_pairs.empty:
        print(df_pairs.head(10).to_string(index=False))

# ─────────────────────────────────────────────────────────────────────────────
# 8. Save outputs
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("STEP 8 – Save outputs")
print("=" * 65)

# Edge list CSV
for metric in ["entropy_similarity", "cosine"]:
    csv_path = str(OUT_DIR / f"edges_{metric}.csv")
    network.save_edge_list(csv_path, metric=metric)

# Similarity matrix CSV (all stored pairs, threshold=0)
for metric in ["entropy_similarity", "cosine"]:
    mat_path = str(OUT_DIR / f"matrix_{metric}.csv")
    network.save_similarity_matrix(mat_path, metric=metric, threshold=0.0)

# Save / reload SimilarityMatrix
npz_path = str(OUT_DIR / "entropy_similarity_matrix.npz")
network.similarity_matrices["entropy_similarity"].save(npz_path)
reloaded = SimilarityMatrix.load(npz_path)
print(f"\n  Reloaded matrix: {reloaded}")
assert reloaded.n_spectra == network.similarity_matrices["entropy_similarity"].n_spectra
print("  ✓ Save/load round-trip OK")
