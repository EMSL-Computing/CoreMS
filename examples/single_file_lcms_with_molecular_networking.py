"""
Single-File LC-MS with Molecular Networking
============================================

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
    STEP 9  – MS2 spectral library search (MSPInterface → fe_search, Fe features only) [if RUN_MOLECULAR_NETWORKING]
    STEP 10 – Extract MS2 spectra from Fe-containing mass features → query lists  [if RUN_MOLECULAR_NETWORKING]
    STEP 10b– Save EIC/MS1/MS2 plots for Fe-containing mass features [if RUN_MOLECULAR_NETWORKING and PLOT_FE_MASS_FEATURES]
    STEP 11 – Build MolecularNetwork (open + neutral_loss)            [if RUN_MOLECULAR_NETWORKING]
    STEP 12 – Summary

Run from the repo root:
    python examples/single_file_lcms_with_molecular_networking.py

# Some hits we are expecting to see in the LCMS data (based on publication):
rt 16.1 min, m/zs at 670.151, 721.0677, 723.0630
molecular formulas C30H27N3O15 
"""

import sys
import numpy as np
from pathlib import Path

# ── Repo root on path ─────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

# =============================================================================
# STEP 1 – Config / paths
# =============================================================================
print("=" * 65)
print("STEP 1 – Config / paths")
print("=" * 65)

# ── Input file ────────────────────────────────────────────────────────────────
RAW_FILE = Path(
    "/Volumes/LaCie/boiteau_data/Prosser soil enrichments/RMB_CWD_180608_prosserM9enrich_hrms2_5.raw"
)

# ── MSP spectral library (same as molecular_networking_demo.py) ───────────────
# Falls back to fe_lib=None (query-vs-query only) if file is missing.
MSP_FILE = REPO_ROOT / "tmp_data" / "20250407_database.msp"

# ── Feature flags ─────────────────────────────────────────────────────────────
# Set to False to skip the molecular networking steps (Steps 9–11) while
# debugging the upstream LCMS processing pipeline.
RUN_MOLECULAR_NETWORKING = True

# Set to False to skip saving EIC/MS1/MS2 plots for Fe-containing mass features.
PLOT_FE_MASS_FEATURES = True

# ── Output directory ──────────────────────────────────────────────────────────
OUT_DIR = REPO_ROOT / "temp.corems" / "single_file_lcms_networking"
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

# Pull MS1 spectra into the internal dataframe (lazy – no MassSpectrum objects yet)
myLCMSobj = parser.get_lcms_obj(spectra="ms1")

print(f"  Polarity : {myLCMSobj.polarity}")
print(f"  Scans    : {len(myLCMSobj.scan_df)}")
ms1_count = (myLCMSobj.scan_df.ms_level == 1).sum()
ms2_count = (myLCMSobj.scan_df.ms_level == 2).sum()
print(f"  MS1 scans: {ms1_count}   MS2 scans: {ms2_count}")

# =============================================================================
# STEP 3 – Set parameters inline
# =============================================================================
# NOTE: For a production workflow, we recommend setting parameters via YAML config files
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

print("  Running find_mass_features …")
myLCMSobj.find_mass_features()
print(f"  Found {len(myLCMSobj.mass_features)} mass features")

print("  Integrating mass features …")
myLCMSobj.integrate_mass_features(drop_if_fail=True)
print(f"  After integration: {len(myLCMSobj.mass_features)} mass features")

print("  Adding peak shape metrics …")
myLCMSobj.add_peak_metrics()

# =============================================================================
# STEP 5 – Add associated MS1 spectra
# =============================================================================
print("\n" + "=" * 65)
print("STEP 5 – Add associated MS1 spectra")
print("=" * 65)

if ms1_is_centroid:
    # Centroided Thermo RAW: use parser for centroid spectra
    myLCMSobj.add_associated_ms1(
        auto_process=True, use_parser=True, spectrum_mode="centroid"
    )
else:
    # Profile mode: reconstruct from raw data
    myLCMSobj.add_associated_ms1(
        auto_process=True, use_parser=False, spectrum_mode="profile"
    )

print(f"  Mass features with MS1 spectra: {sum(1 for mf in myLCMSobj.mass_features.values() if mf.mass_spectrum is not None)}")

# =============================================================================
# STEP 6 – Remove unprocessed data (free memory)
# =============================================================================
print("\n" + "=" * 65)
print("STEP 6 – Remove unprocessed data")
print("=" * 65)

myLCMSobj.remove_unprocessed_data()
print("  Done.")

# =============================================================================
# STEP 7 – Molecular formula search
# =============================================================================
print("\n" + "=" * 65)
print("STEP 7 – Molecular formula search (SearchMolecularFormulasLC)")
print("=" * 65)

from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulasLC  # noqa: E402

mol_search = SearchMolecularFormulasLC(myLCMSobj)
mol_search.run_mass_feature_search()
print("  Molecular formula search complete.")

# Quick summary
n_assigned = sum(
    1 for mf in myLCMSobj.mass_features.values()
    if mf.mass_spectrum is not None and mf.ms1_peak is not None
    and len(mf.ms1_peak.molecular_formulas) > 0
)
print(f"  Mass features with at least one formula: {n_assigned}")

# =============================================================================
# STEP 7b – Identify Fe-containing mass features (putative siderophores)
# =============================================================================
print("\n" + "=" * 65)
print("STEP 7b – Identify Fe-containing mass features")
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
# STEP 7c – Comprehensive export (HDF5 + report CSV)
# =============================================================================
print("\n" + "=" * 65)
print("STEP 7c – Comprehensive export (HDF5 + report CSV)")
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
# STEP 8 – Add MS2 spectra (auto-detect centroid / profile per scan)
# =============================================================================
print("\n" + "=" * 65)
print("STEP 8 – Add MS2 spectra (add_associated_ms2_dda)")
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
    # STEP 9 – MS2 spectral library search (fe_search)
    # =============================================================================
    print("\n" + "=" * 65)
    print("STEP 9 – MS2 spectral library search")
    print("=" * 65)

    fe_lib = None           # will stay None if MSP file is missing
    metabolite_metadata = {}

    if MSP_FILE.exists():
        from corems.molecular_id.search.database_interfaces import MSPInterface  # noqa: E402

        print(f"  Parsing MSP file: {MSP_FILE.name} …")
        my_msp = MSPInterface(file_path=str(MSP_FILE))

        # Build FlashEntropy library for the correct polarity
        fe_lib, metabolite_metadata = my_msp.get_metabolomics_spectra_library(
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

        # Collect MS2 scan numbers belonging to Fe-containing mass features only
        fe_ms2_scans = []
        for mf_id in fe_mf_ids:
            mf = myLCMSobj.mass_features[mf_id]
            for scan_num in mf.ms2_scan_numbers:
                if scan_num in myLCMSobj._ms:
                    fe_ms2_scans.append(scan_num)
        fe_ms2_scans = list(set(fe_ms2_scans))
        print(f"  MS2 scans from Fe-containing features available for search: {len(fe_ms2_scans)}")

        if len(fe_ms2_scans) > 0:
            myLCMSobj.fe_search(
                scan_list=fe_ms2_scans,
                fe_lib=fe_lib,
                peak_sep_da=0.002,
            )
            print("  fe_search complete.")
        else:
            print("  WARNING: No Fe-feature MS2 scans loaded – skipping fe_search.")

        # ── Re-export with MS2 library annotations now populated ─────────────
        exporter_final = LCMSMetabolomicsExport(str(_out_stem), myLCMSobj)
        exporter_final.to_hdf(overwrite=True)
        exporter_final.report_to_csv(molecular_metadata=metabolite_metadata)
        print(f"  Final HDF5 + annotated report CSV written to: {_out_stem}.corems/")
    else:
        print(f"  MSP file not found ({MSP_FILE}). Skipping library search.")
        print("  Molecular networking will run query-vs-query only (fe_lib=None).")

    # =============================================================================
    # STEP 10 – Extract MS2 spectra from mass features → query lists
    # =============================================================================
    print("\n" + "=" * 65)
    print("STEP 10 – Extract MS2 spectra from mass features")
    print("=" * 65)

    query_spectra = []
    query_ids = []
    query_precursor_mzs = []

    n_total_ms2 = 0
    n_skipped_no_fe = 0

    for mf_id, mf in myLCMSobj.mass_features.items():
        ms2 = mf.best_ms2
        if ms2 is None:
            continue
        if not hasattr(ms2, "mz_exp") or len(ms2.mz_exp) == 0:
            continue
        n_total_ms2 += 1
        if mf_id not in fe_mf_ids:
            n_skipped_no_fe += 1
            continue
        query_spectra.append(ms2)                   # MassSpectrum has .mz_exp and .abundance
        query_ids.append(str(mf_id))
        query_precursor_mzs.append(float(mf.mz))    # MS1 m/z of the mass feature

    print(f"  Mass features with MS2: {n_total_ms2}")
    print(f"  Skipped (no Fe in formula): {n_skipped_no_fe}")
    print(f"  Query spectra extracted (Fe-containing): {len(query_spectra)}")
    if len(query_spectra) == 0:
        print("  WARNING: No MS2 spectra found in mass features.")
        print("  Check that add_associated_ms2_dda ran successfully and the file has DDA scans.")
    else:
        # Print a few examples
        for i in range(min(3, len(query_spectra))):
            print(
                f"    [{i}] mf_id={query_ids[i]}  precursor_mz={query_precursor_mzs[i]:.4f}"
                f"  n_peaks={len(query_spectra[i].mz_exp)}"
            )

    # =============================================================================
    # STEP 10b – Save EIC/MS1/MS2 plots for Fe-containing mass features
    # =============================================================================
    print("\n" + "=" * 65)
    print("STEP 10b – Plot Fe-containing mass features")
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
    # STEP 11 – Build MolecularNetwork
    # =============================================================================
    print("\n" + "=" * 65)
    print("STEP 11 – Build MolecularNetwork (open + neutral_loss)")
    print("=" * 65)

    from corems.molecular_networking import MolecularNetwork, SimilarityMatrix  # noqa: E402

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
            print(f"  ✓ [{metric}] HTML network: {html_path}")

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
# STEP 12 – Summary
# =============================================================================
print("\n" + "=" * 65)
print("DONE")
print("=" * 65)
print(f"  Mass features total       : {len(myLCMSobj.mass_features)}")
print(f"  Query spectra for network : {len(query_spectra)}")
print(f"  Networks built            : {list(networks.keys())}")
print(f"  Outputs written to        : {OUT_DIR}")
