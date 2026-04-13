"""
Molecular Networking Demo Script
=================================

Demonstrates the corems.molecular_networking module end-to-end:
    1. Parse an MSP file and build a FlashEntropy library
    2. Build mock query spectra (with noise) from the library entries
    3. Create a MolecularNetwork (identity search + cosine)
    4. Add spectra in a single search
    5. Query edges, neighbors, and network statistics
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

# Larger lipid library (51 k spectra) – used in STEP 10
LARGE_MSP = REPO_ROOT / "tmp_data" / "20250407_database.msp"

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

print(f"\n  Created {len(all_spectra)} mock spectra")
for i in range(min(3, len(all_spectra))):
    print(f"    [{i}] id={all_ids[i]!r}  precursor_mz={all_precursor_mzs[i]:.4f}  "
          f"n_peaks={len(all_spectra[i].mz_exp)}")

print(f"\n  Will add all mock spectra in a single search: {len(all_spectra)} spectra")

# ─────────────────────────────────────────────────────────────────────────────
# 3. Create the MolecularNetwork (open search with cosine)
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("STEP 3 – Create MolecularNetwork (open search)")
print("=" * 65)

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
)

# ─────────────────────────────────────────────────────────────────────────────
# 4. Add spectra (single run)
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("STEP 4 – Add spectra (single run)")
print("=" * 65)

print("\n  Adding all mock spectra …")
network.query_vs_library(
    query_spectra=all_spectra,
    query_ids=all_ids,
    query_precursor_mzs=all_precursor_mzs,
    query_lib_indices=all_lib_indices,
)
print(f"  After adding: {network}")

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

# Save / reload SimilarityMatrix
npz_path = str(OUT_DIR / "entropy_similarity_matrix.npz")
network.similarity_matrices["entropy_similarity"].save(npz_path)
reloaded = SimilarityMatrix.load(npz_path)
print(f"\n  Reloaded matrix: {reloaded}")
assert reloaded.n_spectra == network.similarity_matrices["entropy_similarity"].n_spectra
print("  ✓ Save/load round-trip OK")