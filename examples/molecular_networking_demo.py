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
    spec_id = "MockSpec_" + str(spec_id)
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
        control_spectra_id = "MockSpec_" + str(control_spectra_id)
        
        all_spectra.append(MockSpectrum(control_mz, control_abun, name=control_name))
        all_ids.append(control_id)
        all_precursor_mzs.append(control_pmz)
        
        print(f"\n  Added positive control: {control_id}")
        print(f"    Exact copy of dataframe row[{control_idx}]:")
        print(f"      - compound_name: {control_row.compound_name}")
        print(f"      - spectra_id: {control_spectra_id}")
        print(f"      - precursor_mz: {control_pmz:.4f}")
        print(f"      - n_peaks: {len(control_peaks)}")
        print("    Note: FE library indices differ from dataframe indices (sorted by precursor m/z)")

print(f"\n  Created {len(all_spectra)} query spectra ({len(all_spectra) - 1} noisy + 1 exact)")
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

    run_stage3 = False

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
    print(f"POSITIVE_CONTROL Validation – {label}")
    print("=" * 65)
    print(f"  POSITIVE_CONTROL is an exact copy of dataframe row[{positive_control_source_idx}]")
    print("  Note: FE library indices may differ from dataframe indices (sorted by precursor m/z)")

    pc_entropy_neighbors = network.get_spectrum_neighbors("POSITIVE_CONTROL", metric="entropy_similarity")
    pc_lib_entropy = [(nid, score) for nid, score in pc_entropy_neighbors if nid.isdigit()]
    if pc_lib_entropy:
        top_entropy_lib_id, top_entropy_score = pc_lib_entropy[0]

        try:
            lib_spec = fe_lib[int(top_entropy_lib_id)]
            lib_pmz = lib_spec.get("precursor_mz", "unknown")
            lib_peaks_count = len(lib_spec.get("peaks", []))
            print(f"\n  Top entropy match: library[{top_entropy_lib_id}]")
            print(f"    - precursor_mz: {lib_pmz:.4f}")
            print(f"    - n_peaks (cleaned): {lib_peaks_count}")
            print(f"    - entropy similarity: {top_entropy_score:.6f}")
        except Exception:
            print(f"\n  Top entropy match: library[{top_entropy_lib_id}] = {top_entropy_score:.6f}")

        if top_entropy_score >= 0.95:
            print("    ✓ EXCELLENT (≥0.95)")
        elif top_entropy_score >= 0.80:
            print("    ✓ GOOD (≥0.80)")
        else:
            print("    ~ MODERATE (<0.80) - Expected ~1.0 for exact copy!")

        cosine_score = cosine_mat.get_similarity("POSITIVE_CONTROL", top_entropy_lib_id)
        print(f"\n  Cosine similarity to library[{top_entropy_lib_id}]: {cosine_score:.6f}")
        if cosine_score >= 0.95:
            print("    ✓ EXCELLENT (≥0.95)")
        elif cosine_score >= 0.80:
            print("    ✓ GOOD (≥0.80)")
        else:
            print("    ~ MODERATE (<0.80) - Expected ~1.0 for exact copy!")

        if abs(top_entropy_score - cosine_score) < 0.1:
            print(f"\n  ✓ Metrics are consistent (diff = {abs(top_entropy_score - cosine_score):.4f})")
        else:
            print(f"\n  ⚠ Metrics differ significantly (diff = {abs(top_entropy_score - cosine_score):.4f})")
            print("     This suggests the metrics are using different cleaned spectra!")
    else:
        print("\n  ✗ No library matches found!")

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

    print("\n" + "=" * 65)
    print(f"VISUALIZE – {label}")
    print("=" * 65)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        for metric in ["entropy_similarity", "cosine"]:
            print(f"\n  Creating visualization for {metric}...")
            edges = network.get_network_edges(metric=metric)

            if not edges:
                print(f"    No edges above threshold for {metric}, skipping visualization.")
                continue

            adjacency_all = {}
            for id1, id2, _ in edges:
                if id1 not in adjacency_all:
                    adjacency_all[id1] = set()
                if id2 not in adjacency_all:
                    adjacency_all[id2] = set()
                adjacency_all[id1].add(id2)
                adjacency_all[id2].add(id1)

            nodes_to_keep = set()
            for query_id in query_id_set:
                if query_id not in adjacency_all:
                    continue
                queue = [query_id]
                visited = {query_id}
                while queue:
                    current = queue.pop(0)
                    nodes_to_keep.add(current)
                    for neighbor in adjacency_all.get(current, []):
                        if neighbor not in visited:
                            visited.add(neighbor)
                            queue.append(neighbor)

            edges = [(id1, id2, score) for id1, id2, score in edges if id1 in nodes_to_keep and id2 in nodes_to_keep]

            if not edges:
                print(f"    No edges connected to queries for {metric}, skipping visualization.")
                continue

            nodes = set()
            for id1, id2, _ in edges:
                nodes.add(id1)
                nodes.add(id2)
            nodes = sorted(nodes)
            n_nodes = len(nodes)

            print(f"    Network has {n_nodes} nodes and {len(edges)} edges (query-connected only)")

            rng = np.random.default_rng(seed=42)
            pos = {node: rng.uniform(-1, 1, size=2) for node in nodes}

            k = 1.0 / np.sqrt(n_nodes)
            iterations = 50
            for iteration in range(iterations):
                forces = {node: np.array([0.0, 0.0]) for node in nodes}

                for i, node1 in enumerate(nodes):
                    for node2 in nodes[i + 1 :]:
                        delta = pos[node1] - pos[node2]
                        dist = np.linalg.norm(delta)
                        if dist > 0:
                            force = k * k / dist
                            direction = delta / dist
                            forces[node1] += direction * force
                            forces[node2] -= direction * force

                for id1, id2, score in edges:
                    delta = pos[id1] - pos[id2]
                    dist = np.linalg.norm(delta)
                    if dist > 0:
                        force = dist * dist / k * score
                        direction = delta / dist
                        forces[id1] -= direction * force
                        forces[id2] += direction * force

                temp = 0.1 * (1.0 - iteration / iterations)
                for node in nodes:
                    force_mag = np.linalg.norm(forces[node])
                    if force_mag > 0:
                        displacement = forces[node] / force_mag * min(force_mag, temp)
                        pos[node] += displacement

            all_pos = np.array(list(pos.values()))
            pos_min = all_pos.min(axis=0)
            pos_max = all_pos.max(axis=0)
            pos_range = pos_max - pos_min
            pos_range[pos_range == 0] = 1
            for node in nodes:
                pos[node] = 2 * (pos[node] - pos_min) / pos_range - 1

            fig, ax = plt.subplots(figsize=(10, 8))

            for id1, id2, score in edges:
                x1, y1 = pos[id1]
                x2, y2 = pos[id2]
                alpha = min(1.0, score)
                width = 0.5 + 2.5 * score
                ax.plot([x1, x2], [y1, y2], "gray", alpha=alpha, linewidth=width, zorder=1)

            query_nodes = [node for node in nodes if node in query_id_set]
            library_nodes = [node for node in nodes if node not in query_id_set]

            if query_nodes:
                query_x = [pos[node][0] for node in query_nodes]
                query_y = [pos[node][1] for node in query_nodes]
                ax.scatter(
                    query_x,
                    query_y,
                    c="#FF6B6B",
                    s=80,
                    alpha=0.9,
                    edgecolors="black",
                    linewidths=1.0,
                    label="Query",
                    zorder=2,
                )

            if library_nodes:
                lib_x = [pos[node][0] for node in library_nodes]
                lib_y = [pos[node][1] for node in library_nodes]
                ax.scatter(
                    lib_x,
                    lib_y,
                    c="#4ECDC4",
                    s=30,
                    alpha=0.7,
                    edgecolors="black",
                    linewidths=0.5,
                    label="Library",
                    zorder=2,
                )

            for node in query_nodes:
                x, y = pos[node]
                label_text = node if len(node) <= 15 else node[:12] + "..."
                ax.text(x, y, label_text, fontsize=6, ha="center", va="center", weight="bold", zorder=3)

            all_x = [pos[node][0] for node in nodes]
            all_y = [pos[node][1] for node in nodes]
            x_margin = (max(all_x) - min(all_x)) * 0.05
            y_margin = (max(all_y) - min(all_y)) * 0.05
            ax.set_xlim(min(all_x) - x_margin, max(all_x) + x_margin)
            ax.set_ylim(min(all_y) - y_margin, max(all_y) + y_margin)
            ax.set_aspect("equal")
            ax.axis("off")

            threshold = network._threshold_for(metric)
            stats = network.get_network_stats(metric=metric)
            title = f"{label}: {metric.replace('_', ' ').title()} Network\n"
            title += f"Threshold: {threshold:.2f} | Nodes: {stats['n_nodes']} | Edges: {stats['n_edges']}"
            ax.set_title(title, fontsize=14, weight="bold", pad=20)
            ax.legend(loc="upper right", fontsize=10, framealpha=0.9)

            fig_path = OUT_DIR / f"{search_type}_network_{metric}.png"
            plt.tight_layout()
            plt.savefig(fig_path, dpi=150, bbox_inches="tight", facecolor="white")
            plt.close(fig)
            print(f"    ✓ Saved visualization to {fig_path}")

        print("\n  ✓ Network visualizations complete")

    except ImportError:
        print("\n  ⚠ matplotlib not available, skipping visualization")
        print("    Install with: pip install matplotlib")
    except Exception as e:
        print(f"\n  ⚠ Visualization failed: {e}")

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
