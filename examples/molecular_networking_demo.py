"""
Molecular Networking Demo Script
=================================

Demonstrates the corems.molecular_networking module end-to-end:
  1. Parse an MSP file and build a FlashEntropy library
  2. Build mock query spectra (with noise) from the library entries
  3. Create a MolecularNetwork (identity search + cosine)
  4. Add spectra in two batches (demonstrating incremental updates)
  5. Query edges, neighbors, and network statistics
  6. Save outputs: edge list CSV, similarity matrix CSV, GraphML, PNG plots

Run from the repo root:
    python examples/molecular_networking_demo.py
"""

import sys
import os
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

# Larger lipid library (51 k spectra) – used in STEP 10
LARGE_MSP = REPO_ROOT / "tmp_data" / "20250407_database.msp"
# Number of query spectra to draw from the library for the all-vs-all demo
LARGE_MSP_N_QUERY = 200

print("=" * 65)
print("STEP 1 – Load MSP library and build FlashEntropy index")
print("=" * 65)

msp = MSPInterface(file_path=str(MSP_FILE))
df = msp._data_frame
print(f"  Parsed {len(df)} spectra from {MSP_FILE.name}")
print(f"  Columns: {list(df.columns)}")

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
all_lib_indices = []   # position of each spectrum in the FE library

for lib_idx, row in enumerate(df.itertuples(index=False)):
    peaks = np.array(row.peaks, dtype=float)
    if len(peaks) == 0:
        continue

    mz_vals = peaks[:, 0]
    abun_vals = peaks[:, 1]

    rng = np.random.default_rng(seed=lib_idx)
    noisy_mz = mz_vals + rng.normal(0, 0.001, size=mz_vals.shape)
    noisy_abun = abun_vals * (1 + rng.normal(0, 0.05, size=abun_vals.shape))
    noisy_abun = np.clip(noisy_abun, 0, None)

    name = getattr(row, "compound_name", None) or f"spectrum_{lib_idx}"
    spec_id = getattr(row, "spectra_id", None) or f"spec_{lib_idx:04d}"
    pmz = float(getattr(row, "precursormz", 0.0) or 0.0)

    all_spectra.append(MockSpectrum(noisy_mz, noisy_abun, name=name))
    all_ids.append(str(spec_id))
    all_precursor_mzs.append(pmz)
    all_lib_indices.append(lib_idx)

print(f"  Created {len(all_spectra)} mock spectra")
for i in range(min(3, len(all_spectra))):
    print(f"    [{i}] id={all_ids[i]!r}  precursor_mz={all_precursor_mzs[i]:.4f}  "
          f"n_peaks={len(all_spectra[i].mz_exp)}")

# Split into two batches
batch1_spectra     = all_spectra[:3]
batch1_ids         = all_ids[:3]
batch1_pmzs        = all_precursor_mzs[:3]
batch1_lib_indices = all_lib_indices[:3]

batch2_spectra     = all_spectra[3:]
batch2_ids         = all_ids[3:]
batch2_pmzs        = all_precursor_mzs[3:]
batch2_lib_indices = all_lib_indices[3:]

print(f"\n  Batch 1: {len(batch1_spectra)} spectra  ids={batch1_ids}")
print(f"  Batch 2: {len(batch2_spectra)} spectra  ids={batch2_ids}")

# ─────────────────────────────────────────────────────────────────────────────
# 3. Create the MolecularNetwork (identity search + cosine)
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("STEP 3 – Create MolecularNetwork (identity search)")
print("=" * 65)

network = MolecularNetwork(
    fe_lib=fe_lib,
    search_type="identity",             # precursor-matched
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
)
print(f"  {network}")

# ─────────────────────────────────────────────────────────────────────────────
# 4. Add spectra in two batches
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("STEP 4 – Add spectra in two batches (incremental update)")
print("=" * 65)

print("\n  Adding Batch 1 …")
network.add_spectra(
    spectra=batch1_spectra,
    spectrum_ids=batch1_ids,
    precursor_mzs=batch1_pmzs,
    lib_indices=batch1_lib_indices,
)
print(f"  After Batch 1: {network}")

print("\n  Adding Batch 2 …")
network.add_spectra(
    spectra=batch2_spectra,
    spectrum_ids=batch2_ids,
    precursor_mzs=batch2_pmzs,
    lib_indices=batch2_lib_indices,
)
print(f"  After Batch 2: {network}")

# ─────────────────────────────────────────────────────────────────────────────
# 5. Query the network
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("STEP 5 – Query the network")
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
        print(f"    {nid!r}  score={score:.4f}")
else:
    print("    None above threshold.")

# ─────────────────────────────────────────────────────────────────────────────
# 6. Inspect the SimilarityMatrix objects directly
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("STEP 6 – Inspect SimilarityMatrix objects")
print("=" * 65)

for metric, mat in network.similarity_matrices.items():
    print(f"\n  {mat}")
    df_pairs = mat.to_dataframe(threshold=0.0)
    print(f"    All stored pairs: {len(df_pairs)}")
    if not df_pairs.empty:
        print(df_pairs.to_string(index=False))

# ─────────────────────────────────────────────────────────────────────────────
# 7. Save outputs
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("STEP 7 – Save outputs")
print("=" * 65)

# Edge list CSV
for metric in ["entropy_similarity", "cosine"]:
    csv_path = str(OUT_DIR / f"edges_{metric}.csv")
    network.save_edge_list(csv_path, metric=metric)

# Similarity matrix CSV (all stored pairs, threshold=0)
for metric in ["entropy_similarity", "cosine"]:
    mat_path = str(OUT_DIR / f"matrix_{metric}.csv")
    network.save_similarity_matrix(mat_path, metric=metric, threshold=0.0)

# GraphML (requires networkx)
try:
    for metric in ["entropy_similarity", "cosine"]:
        gml_path = str(OUT_DIR / f"network_{metric}.graphml")
        network.save_graphml(gml_path, metric=metric)
except ImportError as e:
    print(f"  Skipping GraphML (networkx not installed): {e}")

# Network plot PNG (requires matplotlib + networkx)
try:
    for metric in ["entropy_similarity", "cosine"]:
        png_path = str(OUT_DIR / f"network_{metric}.png")
        network.plot_network(
            metric=metric,
            layout="spring",
            output_file=png_path,
            figsize=(8, 6),
            dpi=100,
        )
except ImportError as e:
    print(f"  Skipping network plot (matplotlib/networkx not installed): {e}")

# Heatmap PNG
try:
    for metric in ["entropy_similarity", "cosine"]:
        hm_path = str(OUT_DIR / f"heatmap_{metric}.png")
        network.plot_similarity_heatmap(
            metric=metric,
            output_file=hm_path,
            figsize=(6, 5),
            dpi=100,
        )
except ImportError as e:
    print(f"  Skipping heatmap (matplotlib not installed): {e}")

# Save / reload SimilarityMatrix
npz_path = str(OUT_DIR / "entropy_similarity_matrix.npz")
network.similarity_matrices["entropy_similarity"].save(npz_path)
reloaded = SimilarityMatrix.load(npz_path)
print(f"\n  Reloaded matrix: {reloaded}")
assert reloaded.n_spectra == network.similarity_matrices["entropy_similarity"].n_spectra
print("  ✓ Save/load round-trip OK")

# ─────────────────────────────────────────────────────────────────────────────
# 8. Demonstrate open search (no precursor_mzs required)
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("STEP 8 – Open search (no precursor m/z required)")
print("=" * 65)

network_open = MolecularNetwork(
    fe_lib=fe_lib,
    search_type="open",
    additional_similarities=["cosine"],
    similarity_thresholds={"entropy_similarity": 0.3, "cosine": 0.3},
    peak_sep_da=0.02,
    ms2_tolerance_da=0.01,             # must be <= peak_sep_da / 2
    entropy_threshold_low=0.05,
    use_parallel=False,
    n_jobs=1,
)

network_open.add_spectra(
    spectra=all_spectra,
    spectrum_ids=all_ids,
    lib_indices=all_lib_indices,
    # precursor_mzs intentionally omitted for open search
)
print(f"  {network_open}")
stats_open = network_open.get_network_stats(metric="entropy_similarity")
print(f"  Open search stats: {stats_open}")

# ─────────────────────────────────────────────────────────────────────────────
# 9. Validate precursor_mzs enforcement
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("STEP 9 – Validate precursor_mzs enforcement")
print("=" * 65)

try:
    bad_network = MolecularNetwork(fe_lib=fe_lib, search_type="identity")
    bad_network.add_spectra(
        spectra=batch1_spectra,
        spectrum_ids=batch1_ids,
        # precursor_mzs intentionally omitted
    )
    print("  ERROR: Should have raised ValueError!")
except ValueError as e:
    print(f"  ✓ Correctly raised ValueError: {e}")

# ─────────────────────────────────────────────────────────────────────────────
# 10. Larger lipid library demo (FAMLS – 482 spectra)
#     Uses the FAMLS MSP file if available; skips gracefully if not found.
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("STEP 10 – Larger lipid library (FAMLS, ~482 spectra)")
print("=" * 65)

if not FAMLS_MSP.exists():
    print(f"  FAMLS MSP not found at {FAMLS_MSP} – skipping STEP 10.")
else:
    import time

    # ── Parse the FAMLS MSP ───────────────────────────────────────────────────
    msp_famls = MSPInterface(file_path=str(FAMLS_MSP))
    df_famls = msp_famls._data_frame
    print(f"  Parsed {len(df_famls)} spectra from {FAMLS_MSP.name}")

    # ── Build FlashEntropy library ────────────────────────────────────────────
    fe_famls = msp_famls._to_flashentropy(
        input_dataframe=df_famls,
        normalize=True,
        fe_kwargs={
            "normalize_intensity": True,
            "min_ms2_difference_in_da": 0.02,
            "max_ms2_tolerance_in_da": 0.01,
            "max_indexed_mz": 3000,
            "precursor_ions_removal_da": None,
            "noise_threshold": 0,
        },
    )
    print(f"  FlashEntropy library built ({len(df_famls)} entries indexed)")

    # ── Build mock query spectra from the library ─────────────────────────────
    famls_spectra, famls_ids, famls_pmzs, famls_lib_idx = [], [], [], []
    for lib_idx, row in enumerate(df_famls.itertuples(index=False)):
        peaks = np.array(row.peaks, dtype=float)
        if len(peaks) == 0:
            continue
        mz_vals, abun_vals = peaks[:, 0], peaks[:, 1]
        rng = np.random.default_rng(seed=lib_idx + 1000)
        noisy_mz = mz_vals + rng.normal(0, 0.001, size=mz_vals.shape)
        noisy_abun = np.clip(abun_vals * (1 + rng.normal(0, 0.05, size=abun_vals.shape)), 0, None)

        # Try common column names for spectrum ID and precursor m/z
        spec_id = (
            getattr(row, "spectra_id", None)
            or getattr(row, "name", None)
            or f"famls_{lib_idx:04d}"
        )
        pmz = float(
            getattr(row, "precursormz", None)
            or getattr(row, "precursor_mz", None)
            or 0.0
        )
        famls_spectra.append(MockSpectrum(noisy_mz, noisy_abun))
        famls_ids.append(str(spec_id))
        famls_pmzs.append(pmz)
        famls_lib_idx.append(lib_idx)

    print(f"  Built {len(famls_spectra)} mock query spectra")

    # ── Create network (open search – no precursor required) ──────────────────
    t0 = time.time()
    net_famls = MolecularNetwork(
        fe_lib=fe_famls,
        search_type="open",
        additional_similarities=["cosine"],
        similarity_thresholds={"entropy_similarity": 0.5, "cosine": 0.5},
        peak_sep_da=0.02,
        ms2_tolerance_da=0.01,
        entropy_threshold_low=0.3,
        use_parallel=True,
        n_jobs=-1,
    )

    # Add in two batches to demonstrate incremental update
    half = len(famls_spectra) // 2
    net_famls.add_spectra(
        spectra=famls_spectra[:half],
        spectrum_ids=famls_ids[:half],
        lib_indices=famls_lib_idx[:half],
    )
    net_famls.add_spectra(
        spectra=famls_spectra[half:],
        spectrum_ids=famls_ids[half:],
        lib_indices=famls_lib_idx[half:],
    )
    elapsed = time.time() - t0
    print(f"  Network built in {elapsed:.1f}s  →  {net_famls}")

    # ── Report stats ──────────────────────────────────────────────────────────
    for metric in ["entropy_similarity", "cosine"]:
        stats = net_famls.get_network_stats(metric=metric)
        edges = net_famls.get_network_edges(metric=metric)
        print(f"\n  [{metric}]  stats={stats}")
        if edges:
            print(f"    Top 5 edges:")
            for id1, id2, score in edges[:5]:
                print(f"      {id1!r} ↔ {id2!r}  score={score:.4f}")

    # ── Save outputs ──────────────────────────────────────────────────────────
    famls_out = OUT_DIR / "famls"
    famls_out.mkdir(exist_ok=True)

    for metric in ["entropy_similarity", "cosine"]:
        net_famls.save_edge_list(str(famls_out / f"edges_{metric}.csv"), metric=metric)
        net_famls.save_similarity_matrix(
            str(famls_out / f"matrix_{metric}.csv"), metric=metric, threshold=0.0
        )

    try:
        for metric in ["entropy_similarity", "cosine"]:
            net_famls.save_graphml(str(famls_out / f"network_{metric}.graphml"), metric=metric)
    except ImportError as e:
        print(f"  Skipping GraphML: {e}")

    try:
        for metric in ["entropy_similarity", "cosine"]:
            net_famls.plot_network(
                metric=metric,
                layout="spring",
                output_file=str(famls_out / f"network_{metric}.png"),
                figsize=(12, 10),
                dpi=100,
            )
    except ImportError as e:
        print(f"  Skipping network plot: {e}")

    try:
        for metric in ["entropy_similarity", "cosine"]:
            net_famls.plot_similarity_heatmap(
                metric=metric,
                output_file=str(famls_out / f"heatmap_{metric}.png"),
                figsize=(10, 9),
                dpi=100,
            )
    except ImportError as e:
        print(f"  Skipping heatmap: {e}")

    print(f"\n  FAMLS outputs written to: {famls_out}")

print("\n" + "=" * 65)
print("DEMO COMPLETE")
print(f"Outputs written to: {OUT_DIR}")
print("=" * 65)
