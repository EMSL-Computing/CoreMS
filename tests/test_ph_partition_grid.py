"""Persistent-homology peak picking keeps one m/z grid across scan windows."""

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import KDTree

from corems.encapsulation.factory.parameters import LCMSParameters
from corems.mass_spectra.calc.lc_calc import PHCalculations, find_closest


def _unaligned_profile(n_scans=40, n_bins=900, bin_width=1.5e-4, mz0=200.0):
    """Profile points that are not on one m/z axis, with three chromatographic peaks."""
    jitter = 0.25 * bin_width
    mz_rows, scan_rows, int_rows = [], [], []
    centers = (
        (n_scans * 0.35, n_bins * 0.25, 8.0, 2.0e6),
        (n_scans * 0.55, n_bins * 0.50, 12.0, 1.2e6),
        (n_scans * 0.72, n_bins * 0.75, 6.0, 7.0e5),
    )
    for scan_i in range(n_scans):
        shift = jitter * np.sin(scan_i / 3.0)
        mz = mz0 + np.arange(n_bins) * bin_width + shift
        intensity = np.full(n_bins, 50.0)
        for scan_c, bin_c, width, height in centers:
            rt_env = np.exp(-0.5 * ((scan_i - scan_c) / width) ** 2)
            start = max(0, int(bin_c) - 3)
            stop = min(n_bins, int(bin_c) + 4)
            bump = np.exp(-0.5 * ((np.arange(start, stop) - bin_c) / 1.5) ** 2)
            intensity[start:stop] += height * bump * rt_env
        mz_rows.append(mz)
        scan_rows.append(np.full(n_bins, scan_i + 1))
        int_rows.append(intensity)
    data = pd.DataFrame(
        {
            "mz": np.concatenate(mz_rows),
            "scan": np.concatenate(scan_rows).astype(np.int64),
            "intensity": np.concatenate(int_rows),
        }
    )
    scan_df = pd.DataFrame(
        {
            "scan": np.arange(1, n_scans + 1),
            "scan_time": np.arange(n_scans) * 0.008,
        }
    )
    return data, scan_df


class _PhHost(PHCalculations):
    def __init__(self, data, scan_df):
        self.parameters = LCMSParameters()
        self.parameters.lc_ms.ph_inten_min_rel = 1e-6
        self.parameters.lc_ms.ph_persis_min_rel = 1e-6
        self.parameters.lc_ms.ph_smooth_it = 1
        self.parameters.lc_ms.ph_smooth_radius_mz = 0
        self.parameters.lc_ms.ph_smooth_radius_scan = 1
        self.parameters.lc_ms.verbose_processing = False
        self.parameters.lc_ms.mass_feature_cluster_mz_tolerance_rel = 5e-6
        self.parameters.lc_ms.mass_feature_cluster_rt_tolerance = 0.3
        self._ms_unprocessed = {1: data}
        self.scan_df = scan_df
        self.mass_features = {}


def _peak_rows(host):
    rows = [
        (int(feature.apex_scan), float(feature.mz), float(feature.intensity))
        for feature in host.mass_features.values()
    ]
    return sorted(rows, key=lambda row: (-row[2], row[0], row[1]))


def test_persistent_homology_uses_one_mz_grid_across_scan_windows():
    """Apex scans match a single m/z grid when the run is split into windows.

    Gridding each window chooses different bin centers, so the same ion can
    be reported on a different scan.
    """
    data, scan_df = _unaligned_profile()
    picked = _PhHost(data.copy(), scan_df)
    picked.find_mass_features_ph()

    reference = _PhHost(data.copy(), scan_df)
    reference._ms_unprocessed[1] = reference.grid_data(data)
    reference.check_if_grid = lambda frame: True
    reference.find_mass_features_ph()

    assert _peak_rows(picked) == _peak_rows(reference)


def _centers_from_distance_pairs(mz_desc, radius):
    """Kept m/z values from the higher-intensity radius graph.

    An m/z is dropped when a kept, higher-intensity m/z lies within ``radius``.
    """
    mz = np.asarray(mz_desc, dtype=np.float64)
    if len(mz) == 0:
        return mz
    tree = KDTree(mz.reshape(-1, 1))
    distances = sparse.triu(
        tree.sparse_distance_matrix(tree, radius, output_type="coo_matrix"), k=1
    )
    pairs = (
        np.stack((distances.row, distances.col), axis=1)
        if distances.nnz
        else np.empty((0, 2), dtype=int)
    )
    to_drop = []
    while len(pairs) > 0:
        root_parents = np.setdiff1d(np.unique(pairs[:, 0]), np.unique(pairs[:, 1]))
        children = np.unique(pairs[np.isin(pairs[:, 0], root_parents), 1])
        to_drop = np.append(to_drop, children)
        pairs = pairs[~np.isin(pairs[:, 1], to_drop)]
        pairs = pairs[~np.isin(pairs[:, 0], to_drop)]
    keep = np.ones(len(mz), dtype=bool)
    if len(to_drop):
        keep[np.asarray(to_drop, dtype=int)] = False
    return np.sort(mz[keep])


def test_mz_grid_matches_higher_intensity_radius_rule():
    """Grid bin centers are the kept m/z values from one radius on the whole run."""
    rng = np.random.default_rng(0)
    frames = []
    for scan in range(1, 13):
        mz = 300.0 + np.sort(rng.normal(0.0, 0.02, size=40))
        frames.append(
            pd.DataFrame({"mz": mz, "scan": scan, "intensity": rng.random(40) * 1e5})
        )
    data = pd.concat(frames, ignore_index=True)
    host = _PhHost(
        data,
        pd.DataFrame({"scan": np.arange(1, 13), "scan_time": np.arange(12) * 0.01}),
    )
    data_w = data.copy().reset_index(drop=True)
    data_w["mz_diff"] = np.abs(data_w["mz"].diff())
    radius = data_w.groupby("scan")["mz_diff"].min().min() * 0.99999
    ordered = (
        data_w[["mz", "intensity"]]
        .sort_values(by=["intensity"], ascending=False)[["mz"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    centers = _centers_from_distance_pairs(ordered["mz"].to_numpy(), radius)
    gridded = host._grid_data(data)
    snapped = centers[find_closest(centers, np.sort(data["mz"].to_numpy()))]
    assert set(np.round(gridded["mz"], 10)).issubset(set(np.round(centers, 10)))
    assert len(gridded) <= len(data)
    # Every input point snaps to one of those centers, and the grid uses them.
    assert np.isin(np.round(gridded["mz"], 10), np.round(centers, 10)).all()
    assert np.isin(np.round(snapped, 10), np.round(gridded["mz"], 10)).all()
