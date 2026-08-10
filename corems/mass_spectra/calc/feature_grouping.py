"""
Consensus feature grouping: isotopes (Stage 1), adducts/ISF later.

Approach A': bulk KDTree + sparse pair matrix (RT ∩ unit isotope Δm, like
``find_c13_mass_features``) → Pearson **apex height** gate → roll-up labeling
from mono through successive unit steps (so ¹³C₂ only via ¹³C₁). Do not invent
parents.

Quant gate is fixed for Stage 1 (no runtime method switch):

- Correlation: Pearson only (pairwise-complete on samples with both heights > 0)
- Abundance: mass-feature apex ``intensity`` only (not integrated area)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import KDTree
from scipy.stats import pearsonr

from corems.encapsulation.constant import Atoms

GROUP_COLUMNS = (
    "feature_group_id",
    "ion_role",
    "isotope_state",
    "parent_cluster_id",
)

# Fixed quant-gate policy (Stage 1). Not user-selectable switches.
CORR_METHOD = "pearson"
HEIGHT_COL = "intensity"  # apex peak height; not integrated area


@dataclass(frozen=True)
class FeatureGroupParams:
    """Parameters for consensus feature grouping (Stage 1 isotopes).

    ``rt_tol`` and ``mz_tol_ppm`` are taken from collection
    ``alignment_rt_tol`` / ``alignment_mz_tol_ppm`` when built via
    :meth:`from_lcms_collection_settings`.

    Absolute charges from ``min_charge`` through ``max_charge`` (inclusive)
    are tried when matching isotope Δm (spacing = (m_heavy − m_mono) / |z|).

    Correlation is always Pearson on apex ``intensity`` heights
    (see module-level ``CORR_METHOD`` / ``HEIGHT_COL``). There is no setting
    for alternate metrics or for using integrated area.
    """

    rt_tol: float = 0.4
    mz_tol_ppm: float = 5.0
    min_charge: int = 1
    max_charge: int = 1
    isotope_atoms: Tuple[str, ...] = ("C",)
    min_isotope_abundance: float = 0.01
    max_isotope_offset: int = 4
    corr_threshold: float = 0.80
    min_shared_sample_fraction: float = 0.15
    mono_height_fraction: float = 0.3
    partition_size: int = 5000
    cores: int = 1

    def charge_values(self) -> Tuple[int, ...]:
        """Inclusive absolute charges to search, low to high."""
        lo = abs(int(self.min_charge))
        hi = abs(int(self.max_charge))
        if lo > hi:
            lo, hi = hi, lo
        return tuple(range(lo, hi + 1))

    def min_shared_count(self, n_samples: int) -> int:
        """Minimum shared non-zero sample count for the correlation gate."""
        if n_samples < 1:
            return 1
        return max(1, int(np.ceil(self.min_shared_sample_fraction * n_samples)))

    @classmethod
    def from_lcms_collection_settings(cls, settings) -> "FeatureGroupParams":
        """Build params from LCMSCollectionSettings (or duck-typed object).

        RT and m/z tolerances reuse ``alignment_rt_tol`` and
        ``alignment_mz_tol_ppm`` (not separate feature-group settings).
        Charge range uses ``feature_group_min_charge`` /
        ``feature_group_max_charge``.
        """
        return cls(
            rt_tol=float(settings.alignment_rt_tol),
            mz_tol_ppm=float(settings.alignment_mz_tol_ppm),
            min_charge=int(settings.feature_group_min_charge),
            max_charge=int(settings.feature_group_max_charge),
            isotope_atoms=tuple(settings.isotope_atoms),
            min_isotope_abundance=float(
                settings.feature_group_min_isotope_abundance
            ),
            max_isotope_offset=int(settings.feature_group_max_isotope_offset),
            corr_threshold=float(settings.feature_group_corr_threshold),
            min_shared_sample_fraction=float(
                settings.feature_group_min_shared_sample_fraction
            ),
            mono_height_fraction=float(settings.feature_group_mono_height_fraction),
            partition_size=int(settings.feature_group_partition_size),
            cores=int(getattr(settings, "cores", 1)),
        )


def validate_feature_group_params(params: FeatureGroupParams) -> None:
    """Raise ValueError if parameters are invalid."""
    if params.rt_tol <= 0:
        raise ValueError("alignment_rt_tol (feature grouping RT window) must be > 0")
    if params.mz_tol_ppm <= 0:
        raise ValueError("alignment_mz_tol_ppm (feature grouping m/z tol) must be > 0")
    if params.min_charge == 0 or params.max_charge == 0:
        raise ValueError(
            "feature_group_min_charge and feature_group_max_charge must be non-zero"
        )
    if abs(params.min_charge) < 1 or abs(params.max_charge) < 1:
        raise ValueError(
            "feature_group_min/max_charge absolute values must be >= 1"
        )
    if not params.isotope_atoms:
        raise ValueError("isotope_atoms must be non-empty")
    if not (0.0 <= params.min_isotope_abundance <= 1.0):
        raise ValueError(
            "feature_group_min_isotope_abundance must be in [0, 1]"
        )
    if params.max_isotope_offset < 1:
        raise ValueError("feature_group_max_isotope_offset must be >= 1")
    if not (0.0 < params.mono_height_fraction <= 1.0):
        raise ValueError("feature_group_mono_height_fraction must be in (0, 1]")
    if not (0.0 < params.min_shared_sample_fraction <= 1.0):
        raise ValueError(
            "feature_group_min_shared_sample_fraction must be in (0, 1]"
        )
    if params.corr_threshold < -1.0 or params.corr_threshold > 1.0:
        raise ValueError("feature_group_corr_threshold must be in [-1, 1]")
    if params.partition_size < 1:
        raise ValueError("feature_group_partition_size must be >= 1")
    for atom in params.isotope_atoms:
        rare_isotope_entries(
            atom, min_abundance=params.min_isotope_abundance
        )  # raises if unknown / none above floor


def _validate_mono_element(mono_symbol: str) -> None:
    """Raise ValueError if mono_symbol is not a usable Atoms mono element."""
    if mono_symbol not in Atoms.isotopes:
        raise ValueError(
            f"Unknown mono element '{mono_symbol}' for isotope_atoms; "
            "must be a key in Atoms.isotopes"
        )
    if mono_symbol not in Atoms.atomic_masses:
        raise ValueError(
            f"Mono element '{mono_symbol}' missing from Atoms.atomic_masses"
        )


def rare_isotope_entries(
    mono_symbol: str,
    min_abundance: float = 0.01,
) -> Tuple[Tuple[str, float, float], ...]:
    """
    Rare isotopes for a mono element with signed mass deltas from Atoms.

    Includes **every** rare form listed for the element in ``Atoms.isotopes``
    whose natural abundance is at least ``min_abundance`` (from
    ``Atoms.isotopic_abundance``). Multi-isotope elements (e.g. Se) therefore
    contribute multiple Δm targets, not only the first listed rare form.

    Parameters
    ----------
    mono_symbol : str
        Most-abundant isotope symbol (e.g. ``"C"``, ``"Fe"``, ``"Se"``).
    min_abundance : float
        Minimum natural abundance fraction (0–1). Default 0.01. Isotopes missing
        from ``Atoms.isotopic_abundance`` are skipped.

    Returns
    -------
    tuple of (rare_label, signed_delta, abundance)
        ``signed_delta = m(rare) - m(mono)``. Positive for heavier rare forms
        (¹³C), negative when the listed rare isotope is lighter (⁵⁴Fe).
        Ordered by decreasing abundance (then by |signed_delta|).
    """
    _validate_mono_element(mono_symbol)
    heavies = Atoms.isotopes[mono_symbol][1]
    out = []
    for rare in heavies:
        if rare is None:
            continue
        if rare not in Atoms.atomic_masses:
            raise ValueError(
                f"Rare isotope '{rare}' for '{mono_symbol}' missing from "
                "Atoms.atomic_masses"
            )
        abun = Atoms.isotopic_abundance.get(rare)
        if abun is None:
            continue
        if float(abun) < float(min_abundance):
            continue
        signed = Atoms.atomic_masses[rare] - Atoms.atomic_masses[mono_symbol]
        if signed == 0:
            continue
        out.append((rare, float(signed), float(abun)))
    if not out:
        raise ValueError(
            f"Element '{mono_symbol}' has no rare isotope in Atoms with "
            f"natural abundance >= {min_abundance}"
        )
    out.sort(key=lambda t: (-t[2], abs(t[1])))
    return tuple(out)


def _heavy_isotope_label(
    mono_symbol: str, min_abundance: float = 0.01
) -> str:
    """Primary rare isotope label (highest abundance above the floor)."""
    return rare_isotope_entries(mono_symbol, min_abundance=min_abundance)[0][0]


def isotope_mass_delta(
    mono_symbol: str, charge: int = 1, min_abundance: float = 0.01
) -> float:
    """
    Signed mass difference (primary rare − mono) / |charge| from Atoms.

    Primary rare = highest-abundance rare form meeting ``min_abundance``.
    May be negative (e.g. ⁵⁴Fe − ⁵⁶Fe). Never hard-coded.
    """
    _rare, signed, _ab = rare_isotope_entries(
        mono_symbol, min_abundance=min_abundance
    )[0]
    return signed / abs(int(charge))


def isotope_state_label(rare_or_mono_symbol: str, n: int) -> str:
    """Build isotope_state string, e.g. ``13C1`` or ``54Fe1``."""
    if n <= 0:
        return "M+0"
    return f"{rare_or_mono_symbol}{n}"


def empty_group_labels(cluster_ids: Sequence) -> pd.DataFrame:
    """Return unlabeled group columns for the given cluster ids."""
    idx = pd.Index(cluster_ids, name="cluster")
    return pd.DataFrame(
        {
            "feature_group_id": pd.Series(pd.NA, index=idx, dtype="Int64"),
            "ion_role": pd.Series(None, index=idx, dtype=object),
            "isotope_state": pd.Series(None, index=idx, dtype=object),
            "parent_cluster_id": pd.Series(pd.NA, index=idx, dtype="Int64"),
        }
    )


def _mz_tol_abs(mz_a: float, mz_b: float, mz_tol_ppm: float) -> float:
    return max(mz_a, mz_b) * mz_tol_ppm * 1e-6


def find_isotope_edges(
    cluster_ids: np.ndarray,
    mz: np.ndarray,
    rt: np.ndarray,
    params: FeatureGroupParams,
) -> pd.DataFrame:
    """
    Find unit-step mono→isotope edges via bulk KDTree + sparse pair matrices.

    Pattern matches ``LCMSBase.find_c13_mass_features``:

    1. Sort features by ascending m/z so ``triu`` means light → heavy.
    2. Sparse RT pairs within ``params.rt_tol``.
    3. Sparse m/z pairs within the largest unit isotope spacing (+ ppm tol).
    4. Keep the intersection; retain only pairs whose Δm matches a **unit**
       Atoms spacing (one rare substitution, not 2×, 3×, …).

    Higher-order isotopologues (¹³C₂, …) are **not** linked as mono→M+n here.
    They appear later by roll-up along successive unit edges
    (mono→¹³C₁→¹³C₂), so you never get M+n without M+(n−1).

    Parent is the chemical mono side of the unit step (lighter for ¹³C, heavier
    for ⁵⁴Fe). Tries each charge in ``min_charge``…``max_charge`` and each rare
    form above the abundance floor.

    Returns
    -------
    DataFrame
        Unit edges only. Columns: parent_idx, child_idx, parent_cluster,
        child_cluster, n (=1), atom, rare_label, charge, abs_dm.
    """
    edge_columns = [
        "parent_idx",
        "child_idx",
        "parent_cluster",
        "child_cluster",
        "n",
        "atom",
        "rare_label",
        "charge",
        "abs_dm",
    ]
    n_feat = len(cluster_ids)
    if n_feat < 2:
        return pd.DataFrame(columns=edge_columns)

    # ------------------------------------------------------------------
    # 1) Unit isotope spacings only (signed). Higher n comes from roll-up.
    #    signed_unit = (m_rare - m_mono) / |z|
    # ------------------------------------------------------------------
    spacings = []
    for atom in params.isotope_atoms:
        for rare_label, signed_delta, _abun in rare_isotope_entries(
            atom, min_abundance=params.min_isotope_abundance
        ):
            for z in params.charge_values():
                signed_unit = signed_delta / abs(int(z))
                if signed_unit == 0:
                    continue
                spacings.append((atom, rare_label, int(z), float(signed_unit)))

    if not spacings:
        return pd.DataFrame(columns=edge_columns)

    max_unit = max(abs(s[3]) for s in spacings)

    # ------------------------------------------------------------------
    # 2) Sort by m/z so sparse triu is always light → heavy (same as C13).
    # ------------------------------------------------------------------
    order = np.argsort(mz, kind="mergesort")
    mz_s = np.asarray(mz, dtype=float)[order]
    rt_s = np.asarray(rt, dtype=float)[order]
    # Map sorted row/col back to original array indices
    # order[k] = original index of sorted position k

    # Absolute ppm window used to bound the m/z KD query (generous).
    mz_tol_pad = float(np.max(mz_s)) * params.mz_tol_ppm * 1e-6
    max_mz_gap = max_unit + mz_tol_pad

    # ------------------------------------------------------------------
    # 3) Sparse pair matrices (KDTree.sparse_distance_matrix + triu).
    # ------------------------------------------------------------------
    # RT: pairs with |ΔRT| <= rt_tol
    tree_rt = KDTree(rt_s.reshape(-1, 1))
    sdm_rt = tree_rt.sparse_distance_matrix(
        tree_rt, params.rt_tol, output_type="coo_matrix"
    )
    sdm_rt = sparse.triu(sdm_rt, k=1)
    sdm_rt.data = np.ones_like(sdm_rt.data)

    # m/z: pairs with 0 < Δm/z <= max unit spacing (+ pad)
    tree_mz = KDTree(mz_s.reshape(-1, 1))
    sdm_mz = tree_mz.sparse_distance_matrix(
        tree_mz, max_mz_gap, output_type="coo_matrix"
    )
    sdm_mz = sparse.triu(sdm_mz, k=1)

    # Intersection: keep m/z distances only where RT also coelutes
    cand = sdm_mz.multiply(sdm_rt).tocoo()
    if cand.nnz == 0:
        return pd.DataFrame(columns=edge_columns)

    # ------------------------------------------------------------------
    # 4) Keep pairs whose Δm matches a unit isotope spacing (Atoms).
    # ------------------------------------------------------------------
    rows = []
    # Track (parent_orig, child_orig, atom, rare) to avoid duplicate edges
    seen = set()

    for r, c, dm in zip(cand.row, cand.col, cand.data):
        # r, c are indices into sorted arrays; r is lighter, c is heavier
        i_light = int(order[r])
        i_heavy = int(order[c])
        dm = float(dm)
        tol = _mz_tol_abs(mz_s[r], mz_s[c], params.mz_tol_ppm)

        # Collect all unit spacings that fit; pick smallest residual so near-isobar
        # shifts (e.g. |Δ⁷⁸Se| ≈ |Δ⁸²Se|) disambiguate by exact Atoms mass.
        best = None  # (residual, atom, rare_label, z, signed_unit, parent, child)
        for atom, rare_label, z, signed_unit in spacings:
            unit = abs(signed_unit)
            residual = abs(dm - unit)
            if residual > tol:
                continue
            if signed_unit > 0:
                parent_idx, child_idx = i_light, i_heavy
            else:
                parent_idx, child_idx = i_heavy, i_light
            cand_t = (residual, atom, rare_label, z, signed_unit, parent_idx, child_idx)
            if best is None or residual < best[0]:
                best = cand_t

        if best is None:
            continue

        residual, atom, rare_label, z, signed_unit, parent_idx, child_idx = best
        key = (parent_idx, child_idx, atom, rare_label)
        if key in seen:
            continue
        seen.add(key)

        rows.append(
            {
                "parent_idx": parent_idx,
                "child_idx": child_idx,
                "parent_cluster": cluster_ids[parent_idx],
                "child_cluster": cluster_ids[child_idx],
                "n": 1,  # unit step; depth assigned at roll-up
                "atom": atom,
                "rare_label": rare_label,
                "charge": z,
                "abs_dm": dm,
            }
        )

    if not rows:
        return pd.DataFrame(columns=edge_columns)
    return pd.DataFrame(rows)


def filter_edges_by_height_correlation(
    edges: pd.DataFrame,
    heights: np.ndarray,
    params: FeatureGroupParams,
) -> pd.DataFrame:
    """
    Keep edges with enough shared non-zero samples and Pearson r >= threshold.

    Correlation method is fixed to Pearson (``CORR_METHOD``). Heights must be
    apex peak intensities (``HEIGHT_COL``), not integrated areas.

    Pearson is computed only on samples where **both** heights are > 0
    (pairwise-complete; zeros / missing do not enter the correlation).

    Parameters
    ----------
    heights : ndarray, shape (n_features, n_samples)
        Apex intensity matrix (cluster × sample).
    """
    if edges.empty:
        return edges

    n_samples = int(heights.shape[1]) if heights.ndim == 2 else 1
    min_shared = params.min_shared_count(n_samples)

    keep = []
    for row in edges.itertuples(index=False):
        hi = heights[row.parent_idx]
        hj = heights[row.child_idx]
        # Correlate only where both features were detected (drop zeros / missing)
        shared = (hi > 0) & (hj > 0)
        n_shared = int(np.count_nonzero(shared))
        if n_shared < min_shared:
            keep.append(False)
            continue
        hi_s = hi[shared]
        hj_s = hj[shared]
        if np.std(hi_s) == 0 or np.std(hj_s) == 0:
            keep.append(False)
            continue
        r, _ = pearsonr(hi_s, hj_s)
        if np.isnan(r) or r < params.corr_threshold:
            keep.append(False)
            continue
        keep.append(True)

    return edges.loc[np.asarray(keep)].reset_index(drop=True)


def assign_isotope_labels(
    cluster_ids: np.ndarray,
    mz: np.ndarray,
    rt: np.ndarray,
    heights: np.ndarray,
    edges: pd.DataFrame,
    params: FeatureGroupParams,
) -> pd.DataFrame:
    """
    Assign feature groups by roll-up along unit isotope edges.

    Same idea as ``find_c13_mass_features``:

    - Roots = features that appear as edge parents but never as children
      (chemical monoisotopes of a family).
    - Wave 1: unit children of roots → rare¹ (e.g. ¹³C₁).
    - Wave 2: unit children of wave 1 → rare² (e.g. ¹³C₂), etc.

    Higher-order labels only appear if intermediate steps exist, so there is
    no ¹³C₂ without ¹³C₁. Depth is capped by ``max_isotope_offset``.

    Mono must pass the height prior vs the max height in its rolled-up family.
    """
    labels = empty_group_labels(cluster_ids)
    if edges.empty:
        return labels

    height_scalar = heights.max(axis=1) if heights.ndim == 2 else heights

    # Adjacency: parent_idx -> list of (child_idx, rare_label, atom, z)
    children_of = {}
    all_parents = set()
    all_children = set()
    for row in edges.itertuples(index=False):
        p, c = int(row.parent_idx), int(row.child_idx)
        rare = getattr(row, "rare_label", None) or row.atom
        z = int(getattr(row, "charge", params.charge_values()[0]))
        children_of.setdefault(p, []).append((c, rare, row.atom, z))
        all_parents.add(p)
        all_children.add(c)

    # Roots: monoisotopic candidates (never a unit-step child)
    roots = sorted(all_parents - all_children)
    if not roots:
        return labels

    assigned = set()  # feature indices already labeled into some family
    next_group_id = 0

    for mono_idx in roots:
        if mono_idx in assigned:
            continue

        # Roll-up waves from this root (depth = isotope count for that rare chain)
        # child_idx -> (depth, rare_label, mono_atom, z)
        family_children = {}
        # frontier: list of (node_idx, depth, rare_label used to reach node)
        frontier = [(mono_idx, 0, None)]
        family_nodes = {mono_idx}

        while frontier:
            next_frontier = []
            for node, depth, _rare_in in frontier:
                if depth >= params.max_isotope_offset:
                    continue
                for c, rare_label, atom, z in children_of.get(node, []):
                    if c in family_nodes or c in assigned:
                        continue
                    new_depth = depth + 1
                    if new_depth > params.max_isotope_offset:
                        continue
                    # Prefer closer RT if c already staged with worse RT (rare)
                    if c in family_children:
                        continue
                    family_children[c] = (new_depth, rare_label, atom, z)
                    family_nodes.add(c)
                    next_frontier.append((c, new_depth, rare_label))
            frontier = next_frontier

        if not family_children:
            continue

        # Height prior on chemical mono vs family max
        family_list = [mono_idx] + list(family_children.keys())
        h_max = float(np.max(height_scalar[family_list]))
        if h_max <= 0:
            continue
        if float(height_scalar[mono_idx]) < params.mono_height_fraction * h_max:
            continue

        mono_cluster = int(cluster_ids[mono_idx])
        gid = next_group_id
        next_group_id += 1

        labels.loc[cluster_ids[mono_idx], "feature_group_id"] = gid
        labels.loc[cluster_ids[mono_idx], "ion_role"] = "mono"
        labels.loc[cluster_ids[mono_idx], "isotope_state"] = "M+0"
        labels.loc[cluster_ids[mono_idx], "parent_cluster_id"] = mono_cluster
        assigned.add(mono_idx)

        for c, (depth, rare_label, _atom, _z) in family_children.items():
            cid = cluster_ids[c]
            labels.loc[cid, "feature_group_id"] = gid
            labels.loc[cid, "ion_role"] = "isotope"
            labels.loc[cid, "isotope_state"] = isotope_state_label(rare_label, depth)
            labels.loc[cid, "parent_cluster_id"] = mono_cluster
            assigned.add(c)

    return labels


def group_features_arrays(
    cluster_ids: np.ndarray,
    mz: np.ndarray,
    rt: np.ndarray,
    heights: np.ndarray,
    params: FeatureGroupParams,
) -> pd.DataFrame:
    """
    Run full Stage 1 isotope grouping on array inputs.

    Parameters
    ----------
    cluster_ids : array-like, shape (N,)
    mz, rt : array-like, shape (N,)
    heights : array-like, shape (N, S)
    params : FeatureGroupParams

    Returns
    -------
    DataFrame
        Index = cluster_ids; columns GROUP_COLUMNS.
    """
    validate_feature_group_params(params)

    cluster_ids = np.asarray(cluster_ids)
    mz = np.asarray(mz, dtype=float)
    rt = np.asarray(rt, dtype=float)
    heights = np.asarray(heights, dtype=float)
    if heights.ndim == 1:
        heights = heights.reshape(-1, 1)

    if len(cluster_ids) != len(mz) or len(mz) != len(rt):
        raise ValueError("cluster_ids, mz, and rt must have the same length")
    if heights.shape[0] != len(cluster_ids):
        raise ValueError("heights must have shape (N, S) matching cluster_ids")

    if len(cluster_ids) < 2:
        return empty_group_labels(cluster_ids)

    # Single-process path (partition/multicore can be added without API change)
    edges = find_isotope_edges(cluster_ids, mz, rt, params)
    edges = filter_edges_by_height_correlation(edges, heights, params)
    labels = assign_isotope_labels(cluster_ids, mz, rt, heights, edges, params)
    return labels


def build_height_matrix_from_features(
    mass_features_df: pd.DataFrame,
    cluster_ids: Sequence,
    sample_ids: Sequence,
    induced_df: Optional[pd.DataFrame] = None,
    intensity_col: str = HEIGHT_COL,
    cluster_col: str = "cluster",
    sample_col: str = "sample_id",
) -> np.ndarray:
    """
    Build (N, S) height matrix for ordered cluster_ids × sample_ids.

    Library grouping always passes apex ``intensity`` (``HEIGHT_COL``). There is
    no integrated-area path in ``group_consensus_features``. The ``intensity_col``
    argument is only for offline diagnostics outside the public API.

    Missing cluster/sample pairs are 0.0. If multiple rows exist for the same
    (cluster, sample), the maximum intensity is used.
    """
    frames = [mass_features_df]
    if induced_df is not None and len(induced_df) > 0:
        frames.append(induced_df)
    mf = pd.concat(frames, axis=0, ignore_index=False)
    if cluster_col not in mf.columns:
        raise ValueError(f"mass features missing '{cluster_col}' column")
    if intensity_col not in mf.columns:
        raise ValueError(f"mass features missing '{intensity_col}' column")
    if sample_col not in mf.columns:
        raise ValueError(f"mass features missing '{sample_col}' column")

    mf = mf.dropna(subset=[cluster_col])
    mf = mf.copy()
    mf[cluster_col] = mf[cluster_col].astype(int)

    pivot = (
        mf.groupby([cluster_col, sample_col], as_index=False)[intensity_col]
        .max()
        .pivot(index=cluster_col, columns=sample_col, values=intensity_col)
    )
    # Reindex to full cluster × sample grid
    pivot = pivot.reindex(index=list(cluster_ids), columns=list(sample_ids))
    pivot = pivot.fillna(0.0)
    return pivot.to_numpy(dtype=float)
