"""
Single-File LC-MS with Molecular Networking (research / debug)
==============================================================

**Not the primary public exemplar.** Paths below are placeholders for local
research data. For a fixture-based demo that runs out of the box, use:

- ``examples/molecular_networking_demo.py`` (MSP fixture + mock queries)
- ``examples/molecular_networking_queries_only_demo.py`` (query-only mocks)

End-to-end debug script for a single Thermo RAW DDA file:

    STEP 1  – Config / paths  (RAW_FILE, MSP_FILE, RUN_MOLECULAR_NETWORKING flag)
    STEP 2  – Load raw file → LCMSBase object
    STEP 3  – Set parameters inline
    STEP 4  – Peak picking  (find_mass_features → integrate → peak_metrics)
    STEP 5  – Add associated MS1 spectra
    STEP 6  – Remove unprocessed data (free memory)
    STEP 7  – Molecular formula search (SearchMolecularFormulasLC)
    STEP 7b – Identify Fe-containing mass features (putative siderophores)
    STEP 7c – Comprehensive export (HDF5 + report CSV via LCMSMetabolomicsExport)
    STEP 8  – Add MS2 spectra (add_associated_ms2_dda, auto-detect centroid/profile)
    STEP 10 – Build FlashEntropy library from MSP (for networking)    [if RUN_MOLECULAR_NETWORKING]
    STEP 11 – Save EIC/MS1/MS2 plots for Fe-containing mass features  [if RUN_MOLECULAR_NETWORKING and PLOT_FE_MASS_FEATURES]
    STEP 12 – Build MolecularNetwork (open + neutral_loss)            [if RUN_MOLECULAR_NETWORKING]
    STEP 13 – Summary

Run from the repo root (after editing RAW_FILE / MSP_FILE to local paths):
    python examples/single_file_lcms_with_molecular_networking.py

# Some hits we are expecting to see in the LCMS data (based on publication):
rt 16.1 min, m/zs at 670.151, 721.0677, 723.0630
molecular formulas C30H27N3O15 
"""

import os
import sys
import numpy as np
from pathlib import Path

# =============================================================================
# STEP 1 – Config / paths
# =============================================================================
print("=" * 65)
print("STEP 1 – Config / paths")
print("=" * 65)

# ── Input file (set COREMS_NETWORKING_RAW or edit this path) ─────────────────
RAW_FILE = Path(
    os.environ.get(
        "COREMS_NETWORKING_RAW",
        "/path/to/your_file.raw",
    )
)

# ── MSP spectral library ──────────────────────────────────────────────────────
# Prefer env override; else optional large local library; else repo fixture.
# Falls back to fe_lib=None (query-vs-query only) if file is missing.
_msp_env = os.environ.get("COREMS_NETWORKING_MSP")
if _msp_env:
    MSP_FILE = Path(_msp_env)
elif (Path("tmp_data") / "20250407_database.msp").is_file():
    MSP_FILE = Path("tmp_data") / "20250407_database.msp"
else:
    MSP_FILE = Path("tests/tests_data/lcms/test_db.msp")

# ── Feature flags ─────────────────────────────────────────────────────────────
# Set to False to skip the molecular networking steps (Steps 9–11) while
# debugging the upstream LCMS processing pipeline.
RUN_MOLECULAR_NETWORKING = True

# Set to False to skip saving EIC/MS1/MS2 plots for Fe-containing mass features.
PLOT_FE_MASS_FEATURES = True

# ── Output directory ──────────────────────────────────────────────────────────
OUT_DIR = Path("temp.corems") / "single_file_lcms_networking"
OUT_DIR.mkdir(parents=True, exist_ok=True)

print(f"  RAW file : {RAW_FILE}")
print(f"  MSP file : {MSP_FILE}  (exists={MSP_FILE.exists()})")
print(f"  Output   : {OUT_DIR}")

if not RAW_FILE.exists():
    raise FileNotFoundError(f"RAW file not found: {RAW_FILE}")

# =============================================================================
# STEP 2 – Load raw file → LCMSBase object
# =============================================================================
print("\n" + "=" * 65)
print("STEP 2 – Load raw file → LCMSBase object")
print("=" * 65)

from corems.mass_spectra.input.rawFileReader import ImportMassSpectraThermoMSFileReader  # noqa: E402

parser = ImportMassSpectraThermoMSFileReader(str(RAW_FILE))

# Instantiate LCMSBase with MS1 spectra only
myLCMSobj = parser.get_lcms_obj(spectra="ms1")

# =============================================================================
# STEP 3 – Set parameters inline
# =============================================================================
# NOTE: For a production workflow, we recommend setting parameters via YAML config files
# TODO KRH: add YAML to config loading and use that
print("\n" + "=" * 65)
print("STEP 3 – Set parameters inline")
print("=" * 65)

# ── Persistent-homology peak-picking thresholds ───────────────────────────────
# Raised from defaults (0.001) to speed up demo; lower for production.
myLCMSobj.parameters.lc_ms.ph_inten_min_rel = 0.005
myLCMSobj.parameters.lc_ms.ph_persis_min_rel = 0.005
myLCMSobj.parameters.lc_ms.ph_smooth_it = 0          # no smoothing

# ── MS1 mass-spectrum parameters ──────────────────────────────────────────────
ms1_params = myLCMSobj.parameters.mass_spectrum["ms1"]
ms1_params.mass_spectrum.noise_threshold_method = "relative_abundance"
ms1_params.mass_spectrum.noise_threshold_min_relative_abundance = 0.1
ms1_params.mass_spectrum.noise_min_mz = 0
ms1_params.mass_spectrum.min_picking_mz = 0
ms1_params.mass_spectrum.noise_max_mz = np.inf
ms1_params.mass_spectrum.max_picking_mz = np.inf
ms1_params.ms_peak.legacy_resolving_power = False

# ── Molecular formula search constraints ─────────────────────────────────────
# Expanded to cover expected range (e.g. C30H27N3O15 from publication notes)
ms1_params.molecular_search.url_database = ""   # use local sqlite
ms1_params.molecular_search.usedAtoms = {
    "C": (15, 50),
    "H": (10, 200),
    "O": (0, 25),
    "N": (0, 6),
    "Fe": (0, 2),
}
# Tighten ppm error for these examples
ms1_params.molecular_search.min_ppm_error = -1.0
ms1_params.molecular_search.max_ppm_error = 1.0
#TODO KRH: expand the usedAtoms ranges to encapsulate all expected sideraphore formulas from database

# ── MS2 parameters (copy from MS1, same resolution for Orbitrap DDA) ─────────
ms2_params = ms1_params.copy()
myLCMSobj.parameters.mass_spectrum["ms2"] = ms2_params

print("  Parameters set.")
print(f"  usedAtoms: {ms1_params.molecular_search.usedAtoms}")

# ── Detect if MS1 data are centroided ────────────────────────────────────────
ms1_scan_df = myLCMSobj.scan_df[myLCMSobj.scan_df.ms_level == 1]
ms1_is_centroid = all(x == "centroid" for x in ms1_scan_df.ms_format.tolist())
if ms1_is_centroid:
    print("  MS1 data are centroided → switching peak-picking method")
    myLCMSobj.parameters.lc_ms.peak_picking_method = "centroided_persistent_homology"
    ms1_params.mass_spectrum.noise_threshold_method = "relative_abundance"

# =============================================================================
# STEP 4 – Peak picking
# =============================================================================
print("\n" + "=" * 65)
print("STEP 4 – Peak picking")
print("=" * 65)

myLCMSobj.find_mass_features()
myLCMSobj.integrate_mass_features(drop_if_fail=True)
myLCMSobj.add_peak_metrics()

# =============================================================================
# STEP 5 – Add associated MS1 spectra
# =============================================================================
print("\n" + "=" * 65)
print("STEP 5 – Add associated MS1 spectra")
print("=" * 65)

if ms1_is_centroid:
    # Centroided Thermo RAW: use parser for centroid spectra which will give resovling power
    myLCMSobj.add_associated_ms1(
        auto_process=True, use_parser=True, spectrum_mode="centroid"
    )
else:
    # Profile mode: reconstruct from raw data stored in memory and do the centroiding
    myLCMSobj.add_associated_ms1(
        auto_process=True, use_parser=False, spectrum_mode="profile"
    )

# Remove unprocessed data to free memory before formula search and networking steps.
myLCMSobj.remove_unprocessed_data()

# =============================================================================
# STEP 6 – Molecular formula search
# =============================================================================
print("\n" + "=" * 65)
print("STEP 7 – Molecular formula search (SearchMolecularFormulasLC) on MS1 spectra associated to mass features")
print("=" * 65)

from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulasLC  # noqa: E402

# Instantiate and run formula search on MS1 spectra associated to mass features. 
# Only performs a molecular formula search on MS1 spectra that are associated to mass features.
mol_search = SearchMolecularFormulasLC(myLCMSobj)
mol_search.run_mass_feature_search()

# =============================================================================
# STEP 7 – Identify mass features to those with Fe in any assigned formula
# =============================================================================
print("\n" + "=" * 65)
print("STEP 7 – Identify Fe-containing mass features")
print("=" * 65)

import warnings  # noqa: E402


def _has_fe(mf) -> bool:
    """Return True if any molecular formula assigned to this mass feature contains Fe."""
    ms1_peak = mf.ms1_peak
    if ms1_peak is None:
        return False
    for formula in ms1_peak.molecular_formulas:
        if formula.get("Fe") and formula.get("Fe") > 0:
            return True
    return False


fe_mf_ids = {mf_id for mf_id, mf in myLCMSobj.mass_features.items() if _has_fe(mf)}
print(f"  Fe-containing mass features: {len(fe_mf_ids)} / {len(myLCMSobj.mass_features)}")

# =============================================================================
# STEP 8 - Save intermediate export with MS1 formula annotations (HDF5 + report CSV)
# =============================================================================
print("\n" + "=" * 65)
print("STEP 8 – Save intermediate export with MS1 formula annotations (HDF5 + report CSV)")
print("=" * 65)

from corems.mass_spectra.output.export import LCMSMetabolomicsExport  # noqa: E402

# Stem name for output files (mirrors workflow convention)
_out_stem = OUT_DIR / RAW_FILE.stem

# HDF5 snapshot + full report CSV (no MS2 library metadata yet at this stage)
exporter = LCMSMetabolomicsExport(str(_out_stem), myLCMSobj)
exporter.to_hdf(overwrite=True)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    exporter.report_to_csv()
print(f"  HDF5 + report CSV written to: {_out_stem}.corems/")

# =============================================================================
# STEP 9 – Add MS2 spectra to mass features (add_associated_ms2_dda, auto-detect centroid/profile)
# =============================================================================
print("\n" + "=" * 65)
print("STEP 9 – Add MS2 spectra (add_associated_ms2_dda)")
print("=" * 65)

# Detect whether MS2 scans are centroided or profile.
# spectrum_mode=None lets the parser check each scan individually (handles mixed files).
ms2_scan_df_check = myLCMSobj.scan_df[myLCMSobj.scan_df.ms_level == 2]
if len(ms2_scan_df_check) > 0:
    ms2_formats = ms2_scan_df_check.ms_format.tolist()
    all_centroid = all(x == "centroid" for x in ms2_formats)
    all_profile  = all(x == "profile"  for x in ms2_formats)
    if all_centroid:
        _ms2_mode = "centroid"
        print("  MS2 data: all centroided → spectrum_mode='centroid'")
    elif all_profile:
        _ms2_mode = "profile"
        print("  MS2 data: all profile → spectrum_mode='profile'")
    else:
        _ms2_mode = None   # mixed – let parser decide per scan
        print("  MS2 data: mixed centroid/profile → spectrum_mode=None (auto per scan)")
else:
    _ms2_mode = None
    print("  No MS2 scans found in scan_df.")

myLCMSobj.add_associated_ms2_dda(spectrum_mode=_ms2_mode)

n_with_ms2 = sum(
    1 for mf in myLCMSobj.mass_features.values()
    if len(mf.ms2_scan_numbers) > 0
)
print(f"  Mass features with ≥1 MS2 scan: {n_with_ms2}")
print(f"  Total MS2 spectra loaded: {sum(1 for k in myLCMSobj._ms if myLCMSobj.scan_df.loc[myLCMSobj.scan_df.scan == k, 'ms_level'].values[0] == 2) if len(myLCMSobj._ms) > 0 else 0}")

if RUN_MOLECULAR_NETWORKING:
    # =============================================================================
    # STEP 10 – Build FlashEntropy library from MSP (for networking)
    # =============================================================================
    print("\n" + "=" * 65)
    print("STEP 10 – Build FlashEntropy library from MSP")
    print("=" * 65)

    fe_lib = None           # will stay None if MSP file is missing

    if MSP_FILE.exists():
        from corems.molecular_id.search.database_interfaces import MSPInterface  # noqa: E402

        print(f"  Parsing MSP file: {MSP_FILE.name} …")
        my_msp = MSPInterface(file_path=str(MSP_FILE))

        fe_lib, _ = my_msp.get_metabolomics_spectra_library(
            polarity=myLCMSobj.polarity,
            format="flashentropy",
            normalize=True,
            fe_kwargs={
                "normalize_intensity": True,
                "min_ms2_difference_in_da": 0.02,   # 2× max_ms2_tolerance_in_da
                "max_ms2_tolerance_in_da": 0.01,
                "max_indexed_mz": 3000,
                "precursor_ions_removal_da": None,
                "noise_threshold": 0,
            },
        )
        print(f"  FlashEntropy library built ({len(my_msp._data_frame)} entries, polarity={myLCMSobj.polarity})")
    else:
        print(f"  MSP file not found ({MSP_FILE}). Skipping library build.")
        print("  Molecular networking will run query-vs-query only (fe_lib=None).")

    # =============================================================================
    # STEP 11 – Save EIC/MS1/MS2 plots for Fe-containing mass features
    # =============================================================================
    print("\n" + "=" * 65)
    print("STEP 11 – Plot Fe-containing mass features")
    print("=" * 65)

    if PLOT_FE_MASS_FEATURES:
        import matplotlib  # noqa: E402
        matplotlib.use("Agg")   # non-interactive backend – safe for scripts
        import matplotlib.pyplot as plt  # noqa: E402

        fe_plot_dir = OUT_DIR / "fe_mass_feature_plots"
        fe_plot_dir.mkdir(parents=True, exist_ok=True)

        n_plotted = 0
        for mf_id, mf in myLCMSobj.mass_features.items():
            if mf_id not in fe_mf_ids:
                continue
            try:
                fig = mf.plot(to_plot=["EIC", "MS1", "MS2"], return_fig=True)
                if fig is not None:
                    plot_path = fe_plot_dir / f"mf_{mf_id}_mz{mf.mz:.4f}_rt{mf.retention_time:.2f}.png"
                    fig.savefig(str(plot_path), dpi=150, bbox_inches="tight")
                    plt.close(fig)
                    n_plotted += 1
            except Exception as exc:
                print(f"  WARNING: Could not plot mf_id={mf_id}: {exc}")

        print(f"  Plots saved: {n_plotted}  →  {fe_plot_dir}")
    else:
        print("  [Skipped – PLOT_FE_MASS_FEATURES=False]")

    # =============================================================================
    # STEP 12 – Build MolecularNetwork
    # =============================================================================
    print("\n" + "=" * 65)
    print("STEP 12 – Build MolecularNetwork (open + neutral_loss)")
    print("=" * 65)

    from corems.molecular_networking import MolecularNetwork, SimilarityMatrix  # noqa: E402

    query_spectra, query_ids, query_precursor_mzs = (
        MolecularNetwork.prepare_query_spectra_from_lcms_object(
            myLCMSobj, mf_ids=fe_mf_ids
        )
    )
    print(f"  Query spectra extracted (Fe-containing): {len(query_spectra)}")

    LIBRARY_SIMILARITY_THRESHOLD = 0.3

    networks = {}

    def run_network(search_type: str, label: str):
        print(f"\n  --- {label} ---")

        network = MolecularNetwork(
            fe_lib=fe_lib,                          # None → query-vs-query only
            search_type=search_type,
            additional_similarities=["cosine"],
            similarity_thresholds={
                "entropy_similarity": 0.6,
                "cosine": 0.6,
            },
            use_parallel=False,
            n_jobs=1,
        )

        if fe_lib is not None:
            # Full tiered query: query-vs-query + query-vs-library (+ optional lib-vs-lib)
            network.query_vs_library(
                query_spectra=query_spectra,
                query_ids=query_ids,
                query_precursor_mzs=query_precursor_mzs,
                hydrate_library_similarities=True,
                library_similarity_threshold=LIBRARY_SIMILARITY_THRESHOLD,
            )
        else:
            # No library – query-vs-query only
            network.run_query_vs_query_only(
                query_spectra=query_spectra,
                query_ids=query_ids,
                query_precursor_mzs=query_precursor_mzs,
            )

        print(f"  {network}")

        # ── Stats ─────────────────────────────────────────────────────────────
        query_id_set = set(query_ids)
        for metric in ["entropy_similarity", "cosine"]:
            mat = network.similarity_matrices[metric]
            all_pairs = mat.get_pairs_above_threshold(0.0)
            threshold = network._threshold_for(metric)
            qq_edges = [(a, b, s) for a, b, s in all_pairs
                        if a in query_id_set and b in query_id_set and s >= threshold]
            ql_edges = [(a, b, s) for a, b, s in all_pairs
                        if (a in query_id_set) != (b in query_id_set) and s >= threshold]
            print(f"    [{metric}] Q-Q edges: {len(qq_edges)}  Q-L edges: {len(ql_edges)}")

        # ── Save edge lists + matrices ─────────────────────────────────────────
        for metric in ["entropy_similarity", "cosine"]:
            csv_path = str(OUT_DIR / f"{search_type}_edges_{metric}.csv")
            network.save_edge_list(csv_path, metric=metric)
            mat_path = str(OUT_DIR / f"{search_type}_matrix_{metric}.csv")
            network.save_similarity_matrix(mat_path, metric=metric, threshold=0.0)

        # ── Save NPZ + round-trip check ────────────────────────────────────────
        npz_path = str(OUT_DIR / f"{search_type}_entropy_similarity_matrix.npz")
        network.similarity_matrices["entropy_similarity"].save(npz_path)
        reloaded = SimilarityMatrix.load(npz_path)
        assert reloaded.n_spectra == network.similarity_matrices["entropy_similarity"].n_spectra
        print(f"  ✓ NPZ save/load round-trip OK ({reloaded.n_spectra} spectra)")

        # ── Cluster + HTML network ─────────────────────────────────────────────
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
                f"  ✓ [{metric}] {cluster_summary['n_clusters']} clusters "
                f"across {cluster_summary['n_nodes']} nodes"
            )

            cluster_paths = network.save_network_clusters(
                str(OUT_DIR),
                metric=metric,
                run_id=search_type,
            )
            print(f"  ✓ [{metric}] cluster manifest: {cluster_paths['manifest']}")

            png_path = OUT_DIR / f"{search_type}_network_{metric}.png"
            network.plot_network(
                metric=metric,
                path=str(png_path),
                return_fig=True,
                max_edges=500,
                library_label_field=("compound_name", "name", "spectra_id"),
                bypass_clustering=False,
            )
            print(f"  ✓ [{metric}] static network PNG: {png_path}")

            try:
                import ipysigma  # noqa: F401
            except ImportError:
                print(
                    f"  Skipping interactive HTML for [{metric}] "
                    '(pip install "corems[networking]")'
                )
            else:
                html_path = OUT_DIR / f"{search_type}_network_{metric}.html"
                network.plot_interactive_network(
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
                print(f"  ✓ [{metric}] interactive HTML: {html_path}")

        return network

    if len(query_spectra) == 0:
        print("  Skipping MolecularNetwork – no query spectra available.")
    else:
        networks["open"] = run_network(search_type="open", label="Open Search")
        networks["neutral_loss"] = run_network(search_type="neutral_loss", label="Neutral Loss Search")

else:
    # RUN_MOLECULAR_NETWORKING = False
    fe_lib = None
    query_spectra = []
    networks = {}
    print("\n  [Molecular networking skipped – RUN_MOLECULAR_NETWORKING=False]")

# =============================================================================
# STEP 13 – Summary
# =============================================================================
print("\n" + "=" * 65)
print("DONE")
print("=" * 65)
print(f"  Mass features total       : {len(myLCMSobj.mass_features)}")
print(f"  Query spectra for network : {len(query_spectra)}")
print(f"  Networks built            : {list(networks.keys())}")
print(f"  Outputs written to        : {OUT_DIR}")
