"""
Molecular Networking Demo Script
=================================

This script demonstrates how the molecular networking module will work in CoreMS.
It shows the intended API and workflow for:
  1. Parsing an MSP file and building a FlashEntropy library
  2. Creating a MolecularNetwork from a list of MassSpectrumBase-like objects
  3. Adding spectra incrementally (only computing new similarities)
  4. Querying the network for edges and neighbors
  5. Saving and outputting network diagrams

This script is intended to be run as a standalone demo and will later be
converted into a formal test.

Dependencies (to be added to requirements.txt):
  - scipy (for sparse matrix storage)
  - ms_entropy (already in CoreMS for FlashEntropy)
  - networkx (for graph analysis and layout)
  - matplotlib (for network visualization)

Usage:
    python examples/molecular_networking_demo.py

For larger-scale testing, use:
    MSP_FILE = "/Users/heal742/LOCAL/05_NMDC/02_MetaMS/metams/test_data/test_lcms_metab_data/20250407_database.msp"
"""

import numpy as np
from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
# 1. Parse the MSP file and build a FlashEntropy library
#    MSPInterface reads a local .msp file and can convert it to FlashEntropy
# ─────────────────────────────────────────────────────────────────────────────

from corems.molecular_id.search.database_interfaces import MSPInterface

# Use the test MSP file shipped with CoreMS
MSP_FILE = Path("tests/tests_data/lcms/test_db.msp")

# For larger-scale testing, uncomment:
# MSP_FILE = Path("/Users/heal742/LOCAL/05_NMDC/02_MetaMS/metams/test_data/test_lcms_metab_data/20250407_database.msp")

print(f"Loading MSP library from: {MSP_FILE}")
msp = MSPInterface(file_path=str(MSP_FILE))

print(f"  → Parsed {len(msp._data_frame)} spectra from MSP file")
print(f"  → Columns: {list(msp._data_frame.columns)}")

# Build a FlashEntropy index from the MSP library
fe_lib = msp._to_flashentropy(
    input_dataframe=msp._data_frame,
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

print(f"  → FlashEntropy library built with {len(fe_lib)} spectra")

# ─────────────────────────────────────────────────────────────────────────────
# 2. Simulate query spectra from the MSP library entries
#    In real usage, these would be MassSpectrumBase objects from an LCMS run.
#    Here we create lightweight mock objects that mimic the interface.
# ─────────────────────────────────────────────────────────────────────────────

class MockSpectrum:
    """
    Lightweight mock of MassSpectrumBase for demo purposes.
    In real usage, these would be actual MassSpectrumBase objects from
    an LCMSBase._ms dictionary (e.g., lcms_obj._ms[scan_number]).

    The molecular networking module will accept any object with:
      - .mz_exp     : array-like of m/z values
      - .abundance  : array-like of abundance values
    """
    def __init__(self, mz_exp, abundance, name=None):
        self.mz_exp = np.array(mz_exp, dtype=float)
        self.abundance = np.array(abundance, dtype=float)
        self.name = name  # optional label for visualization

    def to_peaks_array(self):
        """Return spectrum as (N, 2) array of [mz, abundance] pairs."""
        return np.column_stack((self.mz_exp, self.abundance))


# Build mock spectra from the MSP library entries (simulating experimental data)
all_spectra = []
all_ids = []
all_precursor_mzs = []

for i, row in msp._data_frame.iterrows():
    peaks = np.array(row["peaks"], dtype=float)
    if len(peaks) == 0:
        continue

    mz_vals = peaks[:, 0]
    abun_vals = peaks[:, 1]

    # Add small noise to simulate experimental spectra
    rng = np.random.default_rng(seed=i)
    noisy_mz = mz_vals + rng.normal(0, 0.001, size=mz_vals.shape)
    noisy_abun = abun_vals * (1 + rng.normal(0, 0.05, size=abun_vals.shape))
    noisy_abun = np.clip(noisy_abun, 0, None)

    spec = MockSpectrum(
        mz_exp=noisy_mz,
        abundance=noisy_abun,
        name=row.get("compound_name", f"spectrum_{i}"),
    )
    all_spectra.append(spec)
    all_ids.append(row.get("spectra_id", f"spectrum_{i:06d}"))
    # Precursor m/z is stored separately - required for identity/neutral_loss search
    all_precursor_mzs.append(float(row.get("precursormz", 0.0)))

print(f"\nCreated {len(all_spectra)} mock query spectra from MSP library")

# Split into two batches to demonstrate incremental updates
batch1_spectra = all_spectra[:3]
batch1_ids = all_ids[:3]
batch1_precursor_mzs = all_precursor_mzs[:3]

batch2_spectra = all_spectra[3:]
batch2_ids = all_ids[3:]
batch2_precursor_mzs = all_precursor_mzs[3:]

print(f"  Batch 1: {len(batch1_spectra)} spectra")
print(f"    IDs:           {batch1_ids}")
print(f"    Precursor m/z: {batch1_precursor_mzs}")
print(f"  Batch 2: {len(batch2_spectra)} spectra")
print(f"    IDs:           {batch2_ids}")
print(f"    Precursor m/z: {batch2_precursor_mzs}")

# ─────────────────────────────────────────────────────────────────────────────
# 3. Demonstrate the intended MolecularNetwork API
#    (This is the API we will implement in corems/molecular_networking/)
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("INTENDED API DEMONSTRATION")
print("=" * 70)

print("""
# ── Intended usage (once module is implemented) ──────────────────────────

from corems.molecular_networking import MolecularNetwork

# ── FlashEntropy search types ─────────────────────────────────────────────
#
# The 'search_type' parameter controls how FlashEntropy matches spectra:
#
#   "identity"     - Precursor-matched search (ms_entropy method="identity")
#                    Requires precursor_mzs in add_spectra().
#                    Most stringent: only matches spectra with similar precursor m/z.
#
#   "open"         - No precursor matching (ms_entropy method="open")
#                    precursor_mzs not required in add_spectra().
#                    Most permissive: matches any spectra with similar fragment patterns.
#
#   "neutral_loss" - Neutral loss search (ms_entropy method="neutral_loss")
#                    Requires precursor_mzs in add_spectra().
#                    Matches spectra with similar neutral loss patterns.
#
# ── Additional similarity metrics ────────────────────────────────────────
#
# The 'additional_similarities' parameter specifies which extra similarity
# metrics to compute for pairs that pass the entropy similarity threshold.
# Currently supported: ["cosine"]
# Each metric gets its own separate SimilarityMatrix and network.
#
# ── Similarity thresholds ────────────────────────────────────────────────
#
# 'similarity_thresholds' is a dict mapping metric name → threshold value.
# Edges are created in each metric's network when score >= threshold.
# If a metric is not in the dict, a default threshold of 0.5 is used.

# Initialize the network with a pre-built FlashEntropy library
network = MolecularNetwork(
    fe_lib=fe_lib,
    search_type="identity",             # "identity", "open", or "neutral_loss"
    additional_similarities=["cosine"], # Extra metrics to compute alongside entropy
    similarity_thresholds={
        "entropy_similarity": 0.5,      # Threshold for entropy similarity network
        "cosine": 0.6,                  # Threshold for cosine similarity network
    },
    peak_sep_da=0.01,                   # Peak separation for FE search
    ms1_tolerance_da=0.01,             # Precursor m/z tolerance (identity/neutral_loss)
    ms2_tolerance_da=0.005,            # Fragment m/z tolerance for FE search
    use_parallel=True,                  # Enable multiprocessing
    n_jobs=-1,                          # Use all available cores
)

# ── Add first batch of spectra ────────────────────────────────────────────
#
# For "identity" and "neutral_loss" search types, precursor_mzs is REQUIRED.
# Each element of precursor_mzs corresponds to the spectrum at the same index.
# For "open" search type, precursor_mzs is ignored (can be omitted or None).

network.add_spectra(
    spectra=batch1_spectra,
    spectrum_ids=batch1_ids,
    precursor_mzs=batch1_precursor_mzs,   # Required for identity/neutral_loss
)
# → Computes 3*(3-1)/2 = 3 unique pairwise entropy similarities
# → For pairs with entropy_similarity > entropy_threshold_low, also computes cosine
# → Stores each metric in its own separate SimilarityMatrix

# ── Add second batch - only computes NEW similarities ─────────────────────
network.add_spectra(
    spectra=batch2_spectra,
    spectrum_ids=batch2_ids,
    precursor_mzs=batch2_precursor_mzs,   # Required for identity/neutral_loss
)
# → Computes:
#     len(batch2) * len(batch1) = cross-batch pairs
#     len(batch2)*(len(batch2)-1)/2 = within-batch2 pairs
# → Does NOT recompute batch1 vs batch1 (already stored)
# → Each metric's SimilarityMatrix updated independently

# ── Open search (no precursor m/z needed) ────────────────────────────────
network_open = MolecularNetwork(
    fe_lib=fe_lib,
    search_type="open",                 # No precursor matching
    additional_similarities=["cosine"],
    similarity_thresholds={"entropy_similarity": 0.5, "cosine": 0.6},
)

network_open.add_spectra(
    spectra=batch1_spectra,
    spectrum_ids=batch1_ids,
    # precursor_mzs not required for "open" search
)

# ── Query the network ─────────────────────────────────────────────────────

# Get all edges above the similarity threshold for a given metric
edges_entropy = network.get_network_edges(metric="entropy_similarity")
# Returns: list of (id1, id2, score)

edges_cosine = network.get_network_edges(metric="cosine")
# Returns: list of (id1, id2, score)

# Get neighbors for a specific spectrum (default metric: entropy_similarity)
neighbors = network.get_spectrum_neighbors(batch1_ids[0])
# Returns: list of (neighbor_id, score)

neighbors_cosine = network.get_spectrum_neighbors(batch1_ids[0], metric="cosine")
# Returns: list of (neighbor_id, score)

# Get basic network statistics (parameterized by metric)
stats_entropy = network.get_network_stats(metric="entropy_similarity")
# Returns: dict with n_nodes, n_edges, avg_degree, density, etc.

stats_cosine = network.get_network_stats(metric="cosine")

# Access the underlying similarity matrices (one per metric)
sim_matrix_entropy = network.similarity_matrices["entropy_similarity"]
sim_matrix_cosine  = network.similarity_matrices["cosine"]

# Each matrix is a separate SimilarityMatrix object
score = sim_matrix_entropy.get_similarity(batch1_ids[0], batch1_ids[1])
""")

# ─────────────────────────────────────────────────────────────────────────────
# 4. Demonstrate the SimilarityMatrix API (one per metric)
# ─────────────────────────────────────────────────────────────────────────────

print("""
# ── SimilarityMatrix API (one instance per metric) ───────────────────────

from corems.molecular_networking import SimilarityMatrix

# Each metric has its own SimilarityMatrix instance
# (accessed via network.similarity_matrices["entropy_similarity"], etc.)

# Create a standalone SimilarityMatrix
sim_matrix = SimilarityMatrix(metric_name="entropy_similarity")

# Register spectra with user-provided IDs
sim_matrix.register_spectra(spectrum_ids=['lipid_001', 'lipid_002', 'lipid_003'])

# Store a similarity score
sim_matrix.set_similarity('lipid_001', 'lipid_002', score=0.75)

# Retrieve a similarity score
score = sim_matrix.get_similarity('lipid_001', 'lipid_002')

# Get all pairs above a threshold
pairs = sim_matrix.get_pairs_above_threshold(threshold=0.5)
# Returns: list of (id1, id2, score)

# Convert to dense numpy array (for small matrices)
dense = sim_matrix.to_dense()

# Convert to pandas DataFrame
df = sim_matrix.to_dataframe(threshold=0.5)
# Returns: DataFrame with columns ['id1', 'id2', 'score']

# Save/load the matrix
sim_matrix.save('entropy_similarity_matrix.npz')
sim_matrix_loaded = SimilarityMatrix.load('entropy_similarity_matrix.npz')
""")

# ─────────────────────────────────────────────────────────────────────────────
# 5. Demonstrate the SimilarityEngine API
# ─────────────────────────────────────────────────────────────────────────────

print("""
# ── SimilarityEngine API ─────────────────────────────────────────────────

from corems.molecular_networking import SimilarityEngine

# Initialize with a pre-built FlashEntropy library
engine = SimilarityEngine(
    fe_lib=fe_lib,
    search_type="identity",             # "identity", "open", or "neutral_loss"
    additional_similarities=["cosine"], # Extra metrics to compute
    peak_sep_da=0.01,
    ms1_tolerance_da=0.01,
    ms2_tolerance_da=0.005,
    entropy_threshold_low=0.1,          # Min entropy score to trigger additional metrics
    use_parallel=True,
    n_jobs=-1,
)

# Compute all-vs-all similarities for a list of spectra
# Returns: dict of {metric_name: {(id1, id2): score}}
similarities = engine.compute_all_vs_all(
    spectra=batch1_spectra,
    spectrum_ids=batch1_ids,
    precursor_mzs=batch1_precursor_mzs,   # Required for identity/neutral_loss
)
# similarities["entropy_similarity"][(id1, id2)] = 0.75
# similarities["cosine"][(id1, id2)] = 0.82

# Compute similarities between new spectra and existing spectra
# (incremental update - only computes new pairs)
new_similarities = engine.compute_new_vs_existing(
    new_spectra=batch2_spectra,
    new_ids=batch2_ids,
    new_precursor_mzs=batch2_precursor_mzs,       # Required for identity/neutral_loss
    existing_spectra=batch1_spectra,
    existing_ids=batch1_ids,
    existing_precursor_mzs=batch1_precursor_mzs,  # Required for identity/neutral_loss
)
""")

# ─────────────────────────────────────────────────────────────────────────────
# 6. Demonstrate the network visualization and output API
# ─────────────────────────────────────────────────────────────────────────────

print("""
# ── Network Visualization and Output API ─────────────────────────────────

# Get a networkx graph for a specific metric
G_entropy = network.to_networkx(metric="entropy_similarity")
G_cosine  = network.to_networkx(metric="cosine")
# Returns: networkx.Graph with:
#   - nodes: spectrum IDs (with 'name', 'precursor_mz' attributes if available)
#   - edges: (id1, id2, {'score': ...})

# Save the network as a GraphML file (compatible with Cytoscape)
network.save_graphml('my_network_entropy.graphml', metric="entropy_similarity")
network.save_graphml('my_network_cosine.graphml', metric="cosine")

# Save the network as a CSV edge list
network.save_edge_list('my_network_entropy_edges.csv', metric="entropy_similarity")

# Save the similarity matrix as a CSV
network.save_similarity_matrix('my_network_entropy_matrix.csv', metric="entropy_similarity")

# Plot the network diagram using matplotlib
network.plot_network(
    metric="entropy_similarity",
    node_label='name',              # Label nodes with compound name
    edge_weight='score',            # Scale edge width by score
    layout='spring',                # networkx layout algorithm
    output_file='my_network_entropy.png',
    figsize=(12, 10),
    dpi=150,
)

# Plot a heatmap of the similarity matrix
network.plot_similarity_heatmap(
    metric="cosine",
    output_file='my_network_cosine_heatmap.png',
    figsize=(10, 8),
)
""")

# ─────────────────────────────────────────────────────────────────────────────
# 7. Demonstrate the two-stage similarity calculation logic
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("TWO-STAGE SIMILARITY CALCULATION LOGIC")
print("=" * 70)

print("""
Stage 1: FlashEntropy (fast, vectorized)
  - Uses the pre-built FlashEntropy index
  - Computes entropy-based similarity for all pairs using the chosen search_type:
      "identity"     → fe_lib.search(..., method={"identity"})["identity_search"]
                        Requires precursor_mzs; matches spectra with similar precursor m/z
      "open"         → fe_lib.search(..., method={"open"})["open_search"]
                        No precursor m/z required; matches any similar fragment patterns
      "neutral_loss" → fe_lib.search(..., method={"neutral_loss"})["neutral_loss_search"]
                        Requires precursor_mzs; matches spectra with similar neutral losses
  - Filters out pairs with entropy_similarity < entropy_threshold_low (e.g., 0.1)
  - Stores results in the "entropy_similarity" SimilarityMatrix

Stage 2: Additional metrics (only for pairs passing Stage 1)
  - For each pair with entropy_similarity > entropy_threshold_low:
      - Computes cosine similarity using CoreMS's SpectralSimilarity class
      - Stores results in the "cosine" SimilarityMatrix
  - Each metric has its own separate SimilarityMatrix

This approach is efficient because:
  - FlashEntropy is highly optimized for fast filtering
  - Additional metrics are only computed for a small fraction of pairs
  - For N spectra, Stage 1 is O(N²) but fast; Stage 2 is O(k) where k << N²
  - Each metric's network can be queried independently
""")

# ─────────────────────────────────────────────────────────────────────────────
# 8. Demonstrate the incremental update logic
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("INCREMENTAL UPDATE LOGIC")
print("=" * 70)

n_existing = len(batch1_spectra)
n_new = len(batch2_spectra)
n_total = n_existing + n_new

pairs_without_incremental = n_total * (n_total - 1) // 2
pairs_with_incremental = n_new * n_existing + n_new * (n_new - 1) // 2
savings_pct = 100 * (1 - pairs_with_incremental / max(pairs_without_incremental, 1))

print(f"""
When add_spectra() is called with new spectra:

  Existing spectra: {n_existing} (already computed: {n_existing*(n_existing-1)//2} pairs)
  New spectra:      {n_new}

  Required computations:
    - {n_new} × {n_existing} = {n_new * n_existing} cross-batch pairs (new vs existing)
    - {n_new}×({n_new}-1)/2 = {n_new*(n_new-1)//2} within-batch pairs (new vs new)
    Total: {pairs_with_incremental} new pairs

  Skipped computations:
    - {n_existing*(n_existing-1)//2} pairs (batch1 vs batch1, already stored)

  Without incremental: {pairs_without_incremental} total pairs
  With incremental:    {pairs_with_incremental} new pairs
  Savings: {savings_pct:.0f}% fewer computations

  For large datasets (e.g., 1000 existing + 100 new):
    Without incremental: 1100*1099/2 = 604,450 pairs
    With incremental:    100*1000 + 100*99/2 = 104,950 pairs
    Savings: ~83% fewer computations!

  Note: Each metric's SimilarityMatrix is updated independently.
  Note: precursor_mzs must be provided for each batch when using
        "identity" or "neutral_loss" search types.
""")

# ─────────────────────────────────────────────────────────────────────────────
# 9. Demonstrate sparse matrix storage benefits
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("SPARSE MATRIX STORAGE BENEFITS")
print("=" * 70)

print("""
Each SimilarityMatrix uses scipy.sparse for memory-efficient storage:

  - Stores only non-zero similarities (above entropy_threshold_low)
  - For typical molecular networks, most pairs have zero similarity
  - Memory usage scales with number of edges, not N²

  Example memory comparison for 10,000 spectra:
    Dense matrix:  10,000 × 10,000 × 8 bytes = 800 MB
    Sparse matrix: ~1% non-zero → ~8 MB (100x reduction)

  network.similarity_matrices is a dict of separate SimilarityMatrix objects:
    {
        "entropy_similarity": SimilarityMatrix(...),  # entropy scores
        "cosine":             SimilarityMatrix(...),  # cosine scores
    }

  Each SimilarityMatrix stores its data as scipy.sparse.csr_matrix.
  Matrices can be saved/loaded as .npz files for persistence.
""")

print("\n" + "=" * 70)
print("DEMO COMPLETE")
print("=" * 70)
print("\nNext steps:")
print("  1. Implement corems/molecular_networking/similarity_matrix.py")
print("  2. Implement corems/molecular_networking/similarity_engine.py")
print("  3. Implement corems/molecular_networking/network_builder.py")
print("  4. Implement corems/molecular_networking/__init__.py")
print("  5. Convert this demo into a formal test in tests/test_molecular_networking.py")
