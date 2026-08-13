"""
Consensus feature grouping: natural-abundance isotopes + adducts.

**Isotopes:** natural-abundance isotopologues only (e.g. ¹²C/¹³C; rare forms in
``Atoms`` above a natural-abundance floor). Not for tracer/enriched labeling.

**Adducts:** alternate ion forms linked by **pairwise** mass offsets from
``ion_type_dict`` (``corems.mass_spectra.output.export``) among
``feature_group_ion_types`` (e.g. ``[M+H]+`` and ``[M+NH4]+``). No designated
“base” form — any pair with matching Δm can link. Same-analyte forms share one
``feature_group_id`` with isotopes of each form.

Approach: RT ∩ Δm edges (isotope unit steps and/or adduct shifts) → Pearson
**apex height** gate → isotope roll-up → merge across adduct edges.

Quant gate is fixed (no runtime method switch):

- Correlation: Pearson only (pairwise-complete on samples with both heights > 0)
- Abundance: mass-feature apex ``intensity`` only (not integrated area)
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import re
import time
from typing import Dict, MutableMapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import KDTree
from scipy.stats import pearsonr

from corems.encapsulation.constant import Atoms

GROUP_COLUMNS = (
    "feature_group_id",
    "ion_role",
    "ion_type",
    "isotope_state",
    "mono_cluster_id",
)

# Fixed quant-gate policy. Not user-selectable switches.
CORR_METHOD = "pearson"
HEIGHT_COL = "intensity"  # apex peak height; not integrated area

# Default ion forms for adduct linking (keys in ion_type_dict).
# Ordered most → least common (literature frequency). Order is preserved
# after polarity filtering and is used to break exact Δm / neutral-mass ties
# (prefer lower index / more common forms; e.g. [M+H-H2O]+ over [M+H+H2O]+
# when both fit the same spacing to [M+H]+). Mixed pos/neg is OK: opposite
# polarity keys are dropped at group_consensus_features time.
DEFAULT_ION_TYPES: Tuple[str, ...] = (
    "[M+H]+",
    "[M+2H]2+",
    "[M+H-H2O]+",
    "[M-H]-",
    "[M+Na]+",
    "[M+H-NH3]+",
    "[M+NH4]+",
    "[M-H-H2O]-",
    "[M-H+2Na]+",
    "[M-H+H2O]-",
    "[M+NH4-H2O]+",
    "[M+H+H2O]+",
    "[M+H+Na]2+",
    "[M+H+K]2+",
    "[M-2H]2-",
    "[M+2Na]2+",
    "[M+2H-NH3]2+",
    "[M+K]+",
    "[M+H-2H2O]+",
    "[M+3H]3+",
    "[M+2H-H2O]2+",
    "[M]+",
)

PolarityLike = Union[str, int, None]


def normalize_ms_polarity(polarity: PolarityLike) -> Optional[str]:
    """Normalize sample polarity to ``'positive'``, ``'negative'``, or ``None``.

    Accepts common LCMS encodings: ``'positive'`` / ``'negative'``,
    ``1`` / ``-1``, ``'+'`` / ``'-'``, and short forms ``'pos'`` / ``'neg'``.
    Empty or unrecognized values return ``None``.
    """
    if polarity is None:
        return None
    if isinstance(polarity, (int, np.integer)):
        if int(polarity) > 0:
            return "positive"
        if int(polarity) < 0:
            return "negative"
        return None
    s = str(polarity).strip().lower()
    if not s:
        return None
    if s in ("positive", "pos", "+", "1"):
        return "positive"
    if s in ("negative", "neg", "-", "-1"):
        return "negative"
    return None


_ION_TYPE_CHARGE_RE = re.compile(r"(\d*)([+-])\s*$")


def ion_type_polarity(ion_type: str) -> Optional[str]:
    """Infer polarity of an ``ion_type_dict`` key from its trailing charge sign.

    Keys ending in ``+`` / ``2+`` / ``3+`` are positive (e.g. ``[M+H]+``,
    ``[M+2H]2+``); ending in ``-`` / ``2-`` are negative (e.g. ``[M+HCOO]-``,
    ``[M-2H]2-``). Returns ``None`` if no trailing sign is present
    (e.g. bare ``protonated``).
    """
    if ion_type is None:
        return None
    s = str(ion_type).strip()
    if not s:
        return None
    m = _ION_TYPE_CHARGE_RE.search(s)
    if m is None:
        return None
    return "positive" if m.group(2) == "+" else "negative"


def ion_type_charge(ion_type: str) -> int:
    """Absolute charge state encoded in an ion-type key.

    Examples: ``[M+H]+`` → 1, ``[M+2H]2+`` → 2, ``[M-2H]2-`` → 2,
    ``[M+3H]3+`` → 3. Keys without a trailing charge marker default to 1.
    """
    if ion_type is None:
        return 1
    s = str(ion_type).strip()
    m = _ION_TYPE_CHARGE_RE.search(s)
    if m is None:
        return 1
    digits = m.group(1)
    return int(digits) if digits else 1


def max_charge_from_ion_types(ion_types: Sequence[str]) -> int:
    """Largest absolute charge among ion-type keys (at least 1)."""
    if not ion_types:
        return 1
    return max(1, max(ion_type_charge(t) for t in ion_types))


# Allowed mono↔multi (or related multi-charge) ion-type pairs for Pass 2.
# Undirected: frozenset of two ion_type_dict keys. Same-|z| pairs are Pass 1.
SERIES_PAIRS: frozenset = frozenset(
    {
        frozenset(("[M+H]+", "[M+2H]2+")),
        frozenset(("[M+H]+", "[M+3H]3+")),
        frozenset(("[M+Na]+", "[M+2Na]2+")),
        frozenset(("[M+Na]+", "[M+H+Na]2+")),
        frozenset(("[M+K]+", "[M+H+K]2+")),
        frozenset(("[M-H]-", "[M-2H]2-")),
    }
)


def is_allowed_series_pair(type_a: str, type_b: str) -> bool:
    """True if (type_a, type_b) is an allowed mono↔multi series pair."""
    if type_a == type_b:
        return False
    return frozenset((type_a, type_b)) in SERIES_PAIRS


def is_allowed_adduct_type_pair(type_a: str, type_b: str) -> bool:
    """Pass 1 (same |z|) or Pass 2 (series map) type-pair filter."""
    if type_a == type_b:
        return False
    za, zb = ion_type_charge(type_a), ion_type_charge(type_b)
    if za == zb:
        return True
    return is_allowed_series_pair(type_a, type_b)


def filter_ion_types_for_polarity(
    ion_types: Sequence[str],
    polarity: PolarityLike,
) -> Tuple[str, ...]:
    """Keep only ion types whose charge sign matches sample polarity.

    Prevents negative adducts such as ``[M+HCOO]-`` / ``[M+CH3COO]-`` from
    being considered on positive-mode data (and the reverse).

    Parameters
    ----------
    ion_types :
        Candidate ``ion_type_dict`` keys (order preserved).
    polarity :
        Sample/collection polarity. If unknown (``None`` / unrecognized),
        ``ion_types`` is returned unchanged so unit tests and offline array
        callers are not forced to set polarity.

    Returns
    -------
    tuple of str
        Filtered ion types. When polarity is known, types without a trailing
        ``+``/``-`` charge marker are dropped (cannot be assigned safely).
    """
    pol = normalize_ms_polarity(polarity)
    if pol is None:
        return tuple(ion_types)
    kept: list[str] = []
    for it in ion_types:
        sign = ion_type_polarity(it)
        if sign == pol:
            kept.append(it)
    return tuple(kept)


def params_with_polarity_filtered_ion_types(
    params: "FeatureGroupParams",
    polarity: PolarityLike,
) -> "FeatureGroupParams":
    """Return params with ``ion_types`` filtered to ``polarity`` (no-op if unknown)."""
    filtered = filter_ion_types_for_polarity(params.ion_types, polarity)
    if filtered == params.ion_types:
        return params
    return replace(params, ion_types=filtered)


@dataclass(frozen=True)
class FeatureGroupParams:
    """Parameters for consensus feature grouping (isotopes + adducts).

    **Isotopes:** natural-abundance rare forms in ``Atoms`` at or above
    ``min_isotope_abundance``. Not for tracer / labeled experiments.

    **Adducts:** ``ion_types`` is an ordered sequence of ``ion_type_dict`` keys
    (most → least common by default). Edges are sought for every pair via
    neutral-mass consistency (per-type ``|z|`` from the key). There is no
    designated base form and no intensity prior. When two type assignments
    fit equally well, the earlier (more common) types in this sequence win.

    ``rt_tol`` / ``mz_tol_ppm`` come from collection alignment settings via
    :meth:`from_lcms_collection_settings`. Correlation is always Pearson on
    apex ``intensity``.
    """

    rt_tol: float = 0.4
    mz_tol_ppm: float = 5.0
    min_charge: int = 1
    max_charge: int = 1
    # Mono elements for natural-abundance isotope edge search in feature grouping
    # (maps from LCMSCollectionSettings.feature_group_isotope_atoms).
    isotope_atoms: Tuple[str, ...] = ("C",)
    min_isotope_abundance: float = 0.01
    max_isotope_offset: int = 4
    corr_threshold: float = 0.80
    min_shared_sample_fraction: float = 0.15
    # ion_type_dict keys to link; order = most → least common (tie-break)
    ion_types: Tuple[str, ...] = DEFAULT_ION_TYPES
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
    def from_lcms_collection_settings(
        cls,
        settings,
        polarity: PolarityLike = None,
    ) -> "FeatureGroupParams":
        """Build params from LCMSCollectionSettings (or duck-typed object).

        RT and m/z tolerances reuse ``alignment_rt_tol`` and
        ``alignment_mz_tol_ppm`` (not separate feature-group settings).
        Charge range uses ``feature_group_min_charge`` /
        ``feature_group_max_charge``. Mono elements for natural-abundance
        isotope spacing use ``feature_group_isotope_atoms``.

        Parameters
        ----------
        settings :
            ``LCMSCollectionSettings`` (or duck-typed equivalent).
        polarity :
            Optional sample/collection polarity. When provided, ``ion_types``
            are filtered with :func:`filter_ion_types_for_polarity` so that
            e.g. formate/acetate negative adducts are not used on positive
            data. ``group_consensus_features()`` supplies collection polarity.
        """
        ion_types = getattr(
            settings, "feature_group_ion_types", DEFAULT_ION_TYPES
        )
        ion_types = filter_ion_types_for_polarity(tuple(ion_types), polarity)
        # Multi-charge ion types (e.g. [M+2H]2+) expand the isotope charge
        # search so ¹³C spacing is tried at the encoded |z| as well.
        z_from_types = max_charge_from_ion_types(ion_types)
        min_charge = abs(int(settings.feature_group_min_charge))
        max_charge = max(abs(int(settings.feature_group_max_charge)), z_from_types)
        return cls(
            rt_tol=float(settings.alignment_rt_tol),
            mz_tol_ppm=float(settings.alignment_mz_tol_ppm),
            min_charge=min_charge,
            max_charge=max_charge,
            isotope_atoms=tuple(settings.feature_group_isotope_atoms),
            min_isotope_abundance=float(
                settings.feature_group_min_isotope_abundance
            ),
            max_isotope_offset=int(settings.feature_group_max_isotope_offset),
            corr_threshold=float(settings.feature_group_corr_threshold),
            min_shared_sample_fraction=float(
                settings.feature_group_min_shared_sample_fraction
            ),
            ion_types=tuple(ion_types),
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
        raise ValueError(
            "feature_group_isotope_atoms must be non-empty "
            "(mono elements for natural-abundance feature-group isotope edges)"
        )
    if not (0.0 <= params.min_isotope_abundance <= 1.0):
        raise ValueError(
            "feature_group_min_isotope_abundance must be in [0, 1]"
        )
    if params.max_isotope_offset < 1:
        raise ValueError("feature_group_max_isotope_offset must be >= 1")
    if not (0.0 < params.min_shared_sample_fraction <= 1.0):
        raise ValueError(
            "feature_group_min_shared_sample_fraction must be in (0, 1]"
        )
    if params.corr_threshold < -1.0 or params.corr_threshold > 1.0:
        raise ValueError("feature_group_corr_threshold must be in [-1, 1]")
    if params.partition_size < 1:
        raise ValueError("feature_group_partition_size must be >= 1")
    # Validate ion types resolve in ion_type_dict (empty = isotopes only)
    seen_it = set()
    for it in params.ion_types:
        if it in seen_it:
            raise ValueError(
                f"Duplicate ion_type {it!r} in feature_group_ion_types"
            )
        seen_it.add(it)
        _ion_type_mass_offset(it)
    for atom in params.isotope_atoms:
        rare_isotope_entries(
            atom, min_abundance=params.min_isotope_abundance
        )  # raises if unknown / none above floor


def _validate_mono_element(mono_symbol: str) -> None:
    """Raise ValueError if mono_symbol is not a usable Atoms mono element."""
    if mono_symbol not in Atoms.isotopes:
        raise ValueError(
            f"Unknown mono element '{mono_symbol}' in feature_group_isotope_atoms; "
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
    Natural-abundance rare isotopes for a mono element with signed mass deltas.

    Includes **every** rare form listed for the element in ``Atoms.isotopes``
    whose **natural** (terrestrial) abundance is at least ``min_abundance`` (from
    ``Atoms.isotopic_abundance``). Multi-isotope elements (e.g. Se) therefore
    contribute multiple Δm targets, not only the first listed rare form.

    These are natural-abundance isotopologues only—not enriched/tracer labels.

    Parameters
    ----------
    mono_symbol : str
        Most-abundant (natural) isotope symbol (e.g. ``"C"``, ``"Fe"``, ``"Se"``).
    min_abundance : float
        Minimum natural abundance fraction (0–1). Default 0.01. Isotopes missing
        from ``Atoms.isotopic_abundance`` are skipped.

    Returns
    -------
    tuple of (rare_label, signed_delta, abundance)
        ``signed_delta = m(rare) - m(mono)``. Positive for heavier rare forms
        (¹³C), negative when the listed rare isotope is lighter (⁵⁴Fe).
        Ordered by decreasing natural abundance (then by |signed_delta|).
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
    """Primary natural-abundance rare isotope label (highest abundance above the floor)."""
    return rare_isotope_entries(mono_symbol, min_abundance=min_abundance)[0][0]


def isotope_mass_delta(
    mono_symbol: str, charge: int = 1, min_abundance: float = 0.01
) -> float:
    """
    Signed mass difference (primary natural-abundance rare − mono) / |charge|.

    Primary rare = highest natural-abundance rare form meeting ``min_abundance``
    in ``Atoms``. May be negative (e.g. ⁵⁴Fe − ⁵⁶Fe). Never hard-coded.
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
            "ion_type": pd.Series(None, index=idx, dtype=object),
            "isotope_state": pd.Series(None, index=idx, dtype=object),
            "mono_cluster_id": pd.Series(pd.NA, index=idx, dtype="Int64"),
        }
    )


def _mz_tol_abs(mz_a: float, mz_b: float, mz_tol_ppm: float) -> float:
    return max(mz_a, mz_b) * mz_tol_ppm * 1e-6


def _get_ion_type_dict() -> Dict:
    """Lazy import to avoid heavy export module at package import time."""
    from corems.mass_spectra.output.export import ion_type_dict

    return ion_type_dict


def _atom_count_mass(atom_counts: Dict[str, int]) -> float:
    """Exact mass for an atom-count dict using ``Atoms.atomic_masses``."""
    total = 0.0
    for symbol, n in atom_counts.items():
        if n == 0:
            continue
        if symbol not in Atoms.atomic_masses:
            raise ValueError(
                f"Unknown atom {symbol!r} in ion_type_dict entry; "
                "must be a key in Atoms.atomic_masses"
            )
        total += float(Atoms.atomic_masses[symbol]) * int(n)
    return total


def _ion_type_mass_offset(ion_type: str) -> float:
    """
    Neutral-formula mass offset for an ion_type_dict key.

    offset = mass(atoms to add) − mass(atoms to subtract).
    Same convention as ``LCMSMetabolomicsExport.get_ion_formula``.
    """
    ion_type_dict = _get_ion_type_dict()
    if ion_type not in ion_type_dict:
        raise ValueError(
            f"Unknown ion_type {ion_type!r} for feature grouping; "
            f"must be a key in ion_type_dict "
            f"(e.g. '[M+H]+', '[M+NH4]+'). "
            f"Known: {sorted(ion_type_dict.keys())}"
        )
    add_dict, sub_dict = ion_type_dict[ion_type]
    return _atom_count_mass(add_dict) - _atom_count_mass(sub_dict)


def ion_type_mass_delta(
    ion_type_a: str,
    ion_type_b: str,
    charge: int = 1,
) -> float:
    """
    m/z spacing for **same** absolute charge: ion_type_b − ion_type_a.

    Example: ``[M+NH4]+`` − ``[M+H]+`` = m(N) + 3·m(H) at |z|=1.
    Sign indicates which form is heavier at that charge; no preferred “base”.

    For pairs with **different** charges (e.g. ``[M+H]+`` vs ``[M+2H]2+``),
    m/z spacing depends on neutral mass ``M`` and is **not** a constant;
    use :func:`neutral_mass_from_mz` / adduct edge search instead.
    """
    signed = _ion_type_mass_offset(ion_type_b) - _ion_type_mass_offset(
        ion_type_a
    )
    return float(signed) / abs(int(charge))


# Backward-compatible alias
adduct_mass_delta = ion_type_mass_delta


def neutral_mass_from_mz(mz: float, ion_type: str) -> float:
    """Neutral mass implied by observed m/z and an ion-type assignment.

    ``M = |z| * m/z − offset(ion_type)`` with ``offset`` from
    ``ion_type_dict`` atom add/subtract and ``|z|`` from the key suffix
    (e.g. ``2+`` → 2).
    """
    z = ion_type_charge(ion_type)
    offset = _ion_type_mass_offset(ion_type)
    return float(z) * float(mz) - float(offset)


def find_isotope_edges(
    cluster_ids: np.ndarray,
    mz: np.ndarray,
    rt: np.ndarray,
    params: FeatureGroupParams,
) -> pd.DataFrame:
    """
    Find unit-step mono→natural-abundance isotope edges via KDTree pair matrices.

    Stage 1 only: natural-abundance rare isotopes from ``Atoms`` (not tracer
    enrichment). Pattern matches ``LCMSBase.find_c13_mass_features``:

    1. Sort features by ascending m/z so ``triu`` means light → heavy.
    2. Sparse RT pairs within ``params.rt_tol``.
    3. Sparse m/z pairs within the largest unit isotope spacing (+ ppm tol).
    4. Keep the intersection; retain only pairs whose Δm matches a **unit**
       natural-abundance Atoms spacing (one rare substitution, not 2×, 3×, …).

    Higher-order natural isotopologues (¹³C₂, …) are **not** linked as mono→M+n
    here. They appear later by roll-up along successive unit edges
    (mono→¹³C₁→¹³C₂), so you never get M+n without M+(n−1).

    Parent is the chemical mono (most-abundant natural form) side of the unit
    step (lighter for ¹³C, heavier for ⁵⁴Fe). Tries each charge in
    ``min_charge``…``max_charge`` and each rare form above the natural-abundance
    floor.

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


_ADDUCT_EDGE_COLUMNS = [
    "parent_idx",
    "child_idx",
    "parent_cluster",
    "child_cluster",
    "parent_ion_type",
    "child_ion_type",
    "charge",
    "abs_dm",
    "residual",
    "rank_sum",
]


def find_adduct_edges(
    cluster_ids: np.ndarray,
    mz: np.ndarray,
    rt: np.ndarray,
    params: FeatureGroupParams,
) -> pd.DataFrame:
    """
    Find ion-type edges using pairwise neutral-mass consistency.

    ``M = |z| · m/z − offset(ion_type)`` must agree within ppm for the two
    assigned types on a coeluting peak pair.

    **Pass 1:** only type pairs with the same absolute charge (same-|z|).
    **Pass 2:** type pairs in ``SERIES_PAIRS`` (mono↔multi of the same series,
    e.g. ``[M+H]+``↔``[M+2H]2+``, ``[M+Na]+``↔``[M+2Na]2+`` / ``[M+H+Na]2+``).

    Arbitrary cross-charge pairs (e.g. light peak as ``[M+2H-NH3]2+`` of a
    heavy ``[M+H]+`` stranger) are **not** considered.

    When multiple type assignments fit, prefer lower residual then lower
    rank-sum (earlier = more common in ``params.ion_types``).

    Returns
    -------
    DataFrame
        Columns include parent/child indices and ion types, ``residual``,
        ``rank_sum``, ``abs_dm``. Parent = lower m/z of the pair.
    """
    edge_columns = list(_ADDUCT_EDGE_COLUMNS)
    n_feat = len(cluster_ids)
    ion_types = tuple(params.ion_types)
    if n_feat < 2 or len(ion_types) < 2:
        return pd.DataFrame(columns=edge_columns)

    type_rank = {t: i for i, t in enumerate(ion_types)}
    typed: list[Tuple[str, int, float]] = []
    for t in ion_types:
        typed.append((t, ion_type_charge(t), float(_ion_type_mass_offset(t))))

    mz = np.asarray(mz, dtype=float)
    rt = np.asarray(rt, dtype=float)
    order = np.argsort(mz, kind="mergesort")
    rt_s = rt[order]

    tree_rt = KDTree(rt_s.reshape(-1, 1))
    sdm_rt = tree_rt.sparse_distance_matrix(
        tree_rt, params.rt_tol, output_type="coo_matrix"
    )
    sdm_rt = sparse.triu(sdm_rt, k=1)
    if sdm_rt.nnz == 0:
        return pd.DataFrame(columns=edge_columns)

    rows = []
    seen = set()
    ppm = float(params.mz_tol_ppm)

    for r, c in zip(sdm_rt.row, sdm_rt.col):
        i_lo = int(order[r])
        i_hi = int(order[c])
        mz_lo = float(mz[i_lo])
        mz_hi = float(mz[i_hi])
        dm = abs(mz_hi - mz_lo)

        best = None
        for t_lo, z_lo, off_lo in typed:
            M_lo = z_lo * mz_lo - off_lo
            if M_lo <= 0:
                continue
            rank_lo = type_rank[t_lo]
            for t_hi, z_hi, off_hi in typed:
                if not is_allowed_adduct_type_pair(t_lo, t_hi):
                    continue
                M_hi = z_hi * mz_hi - off_hi
                if M_hi <= 0:
                    continue
                residual = abs(M_lo - M_hi)
                tol = max(M_lo, M_hi) * ppm * 1e-6
                if residual > tol:
                    continue
                rank_hi = type_rank[t_hi]
                cand_key = (
                    residual,
                    rank_lo + rank_hi,
                    rank_lo,
                    rank_hi,
                    t_lo,
                    t_hi,
                    z_lo,
                    z_hi,
                )
                if best is None or cand_key < best:
                    best = cand_key

        if best is None:
            continue
        residual, rank_sum, _rlo, _rhi, t_lo, t_hi, z_lo, z_hi = best
        key = (i_lo, i_hi, t_lo, t_hi)
        if key in seen:
            continue
        seen.add(key)
        rows.append(
            {
                "parent_idx": i_lo,
                "child_idx": i_hi,
                "parent_cluster": cluster_ids[i_lo],
                "child_cluster": cluster_ids[i_hi],
                "parent_ion_type": t_lo,
                "child_ion_type": t_hi,
                "charge": int(max(z_lo, z_hi)),
                "abs_dm": dm,
                "residual": float(residual),
                "rank_sum": int(rank_sum),
            }
        )

    if not rows:
        return pd.DataFrame(columns=edge_columns)
    return pd.DataFrame(rows)


def sort_adduct_edges(edges: pd.DataFrame) -> pd.DataFrame:
    """Order edges for merge: residual ↑, rank_sum ↑, then cluster ids."""
    if edges is None or edges.empty:
        return edges if edges is not None else pd.DataFrame(columns=_ADDUCT_EDGE_COLUMNS)
    out = edges.copy()
    if "residual" not in out.columns:
        out["residual"] = 0.0
    if "rank_sum" not in out.columns:
        out["rank_sum"] = 0
    return out.sort_values(
        by=["residual", "rank_sum", "parent_cluster", "child_cluster"],
        kind="mergesort",
    ).reset_index(drop=True)


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
    Assign feature groups by roll-up along unit natural-abundance isotope edges.

    Same idea as ``find_c13_mass_features`` for natural-abundance envelopes:

    - Roots = features that appear as edge parents but never as children
      (chemical monoisotopes of a family).
    - Wave 1: unit children of roots → rare¹ (e.g. natural ¹³C₁).
    - Wave 2: unit children of wave 1 → rare² (e.g. ¹³C₂), etc.

    Higher-order labels only appear if intermediate steps exist, so there is
    no ¹³C₂ without ¹³C₁. Depth is capped by ``max_isotope_offset``.

    Chemical mono (``ion_role="mono"`` / ``M+0``) is the roll-up root from
    geometry (Atoms side of each unit step), not necessarily the tallest peak
    in the envelope. No mono-vs-family height prior is applied; geometry +
    correlation are the gates. Not designed for labeled/enriched isotope series.

    ``heights`` is accepted for API compatibility with the correlation stage but
    is not used during roll-up labeling.
    """
    labels = empty_group_labels(cluster_ids)
    if edges.empty:
        return labels

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

        mono_cluster = int(cluster_ids[mono_idx])
        gid = next_group_id
        next_group_id += 1

        # ion_type filled later from pairwise adduct edges (or single configured type)
        single_ion_type = (
            params.ion_types[0] if len(params.ion_types) == 1 else None
        )
        labels.loc[cluster_ids[mono_idx], "feature_group_id"] = gid
        labels.loc[cluster_ids[mono_idx], "ion_role"] = "mono"
        labels.loc[cluster_ids[mono_idx], "ion_type"] = single_ion_type
        labels.loc[cluster_ids[mono_idx], "isotope_state"] = "M+0"
        labels.loc[cluster_ids[mono_idx], "mono_cluster_id"] = mono_cluster
        assigned.add(mono_idx)

        for c, (depth, rare_label, _atom, _z) in family_children.items():
            cid = cluster_ids[c]
            labels.loc[cid, "feature_group_id"] = gid
            labels.loc[cid, "ion_role"] = "isotope"
            labels.loc[cid, "ion_type"] = single_ion_type
            labels.loc[cid, "isotope_state"] = isotope_state_label(rare_label, depth)
            labels.loc[cid, "mono_cluster_id"] = mono_cluster
            assigned.add(c)

    return labels


def _remap_group_id(labels: pd.DataFrame, old_gid, new_gid) -> None:
    """In-place remap feature_group_id old → new."""
    if pd.isna(old_gid) or pd.isna(new_gid) or old_gid == new_gid:
        return
    mask = labels["feature_group_id"] == old_gid
    labels.loc[mask, "feature_group_id"] = new_gid


def _is_null_ion_type(val) -> bool:
    return val is None or (isinstance(val, float) and np.isnan(val)) or pd.isna(val)


def _form_members_for_endpoint(
    labels: pd.DataFrame,
    endpoint: int,
    target_ion_type: str,
) -> list:
    """
    Isotope form subtree for an adduct-edge endpoint.

    Members share the endpoint's chemical ``mono_cluster_id`` (or the endpoint
    itself) and have ``ion_type`` null or equal to ``target_ion_type``.
    Does **not** include other ion forms already in the same feature_group_id.
    """
    if endpoint not in labels.index:
        return [endpoint]

    mono_id = labels.loc[endpoint, "mono_cluster_id"]
    if pd.isna(mono_id):
        mono_id = endpoint
    else:
        mono_id = int(mono_id)

    members = []
    for cid in labels.index:
        cid_i = int(cid)
        mid = labels.loc[cid, "mono_cluster_id"]
        if pd.isna(mid):
            same_mono = cid_i == endpoint or cid_i == mono_id
        else:
            same_mono = int(mid) == mono_id or cid_i == mono_id
        if not same_mono:
            continue
        it = labels.loc[cid, "ion_type"]
        if _is_null_ion_type(it) or it == target_ion_type:
            members.append(cid_i)
    if endpoint not in members:
        members.append(int(endpoint))
    return members


def _form_mono_id(labels: pd.DataFrame, form_members: Sequence[int], endpoint: int) -> int:
    monos = [
        c
        for c in form_members
        if c in labels.index and labels.loc[c, "ion_role"] == "mono"
    ]
    if monos:
        return int(monos[0])
    if endpoint in form_members:
        return int(endpoint)
    return int(form_members[0])


def merge_adduct_edges_into_labels(
    labels: pd.DataFrame,
    cluster_ids: np.ndarray,
    adduct_edges: pd.DataFrame,
    params: FeatureGroupParams,
    mz_by_cluster: Optional[Dict[int, float]] = None,
) -> pd.DataFrame:
    """
    Merge forms across ion-type edges under correctness constraints.

    - Paint **form subtrees only** (mono + isotopes of that form), never
      re-type other ion forms already in a multi-form group.
    - Require shared neutral mass between form monos (if ``mz_by_cluster``).
    - Reject if an endpoint already has a conflicting non-null ``ion_type``.
    - Edges should be pre-sorted (see :func:`sort_adduct_edges`).
    """
    if adduct_edges is None or adduct_edges.empty:
        return labels

    labels = labels.copy()
    next_gid = 0
    if labels["feature_group_id"].notna().any():
        next_gid = int(labels["feature_group_id"].max()) + 1

    ppm = float(params.mz_tol_ppm)
    if mz_by_cluster is None:
        mz_by_cluster = {}

    def _M(cid: int, ion_type: str) -> Optional[float]:
        if cid not in mz_by_cluster:
            return None
        try:
            return neutral_mass_from_mz(mz_by_cluster[cid], ion_type)
        except Exception:
            return None

    for row in adduct_edges.itertuples(index=False):
        parent_c = int(row.parent_cluster)
        child_c = int(row.child_cluster)
        type_light = row.parent_ion_type
        type_heavy = row.child_ion_type

        if parent_c not in labels.index or child_c not in labels.index:
            continue

        # One interpretation per cluster: reject conflicting types
        it_p = labels.loc[parent_c, "ion_type"]
        it_c = labels.loc[child_c, "ion_type"]
        if not _is_null_ion_type(it_p) and it_p != type_light:
            continue
        if not _is_null_ion_type(it_c) and it_c != type_heavy:
            continue

        light_form = _form_members_for_endpoint(labels, parent_c, type_light)
        heavy_form = _form_members_for_endpoint(labels, child_c, type_heavy)

        # Seed unlabeled endpoints as mono of their form
        for cid in (parent_c, child_c):
            if _is_null_ion_type(labels.loc[cid, "ion_role"]):
                labels.loc[cid, "ion_role"] = "mono"
                labels.loc[cid, "isotope_state"] = "M+0"

        light_mono = _form_mono_id(labels, light_form, parent_c)
        heavy_mono = _form_mono_id(labels, heavy_form, child_c)

        # Shared-M check between form monos
        M_l = _M(light_mono, type_light)
        M_h = _M(heavy_mono, type_heavy)
        if M_l is not None and M_h is not None:
            if M_l <= 0 or M_h <= 0:
                continue
            tol = max(M_l, M_h) * ppm * 1e-6
            if abs(M_l - M_h) > tol:
                continue

        g_p = labels.loc[parent_c, "feature_group_id"]
        g_c = labels.loc[child_c, "feature_group_id"]

        if pd.isna(g_p) and pd.isna(g_c):
            keep_gid = next_gid
            next_gid += 1
        elif pd.isna(g_p):
            keep_gid = int(g_c)
        elif pd.isna(g_c):
            keep_gid = int(g_p)
        else:
            keep_gid = int(g_p)
            if int(g_c) != keep_gid:
                _remap_group_id(labels, int(g_c), keep_gid)

        # Ensure form monos are mono
        for mid in (light_mono, heavy_mono):
            if mid in labels.index and labels.loc[mid, "ion_role"] != "isotope":
                labels.loc[mid, "ion_role"] = "mono"
                if _is_null_ion_type(labels.loc[mid, "isotope_state"]):
                    labels.loc[mid, "isotope_state"] = "M+0"

        for cid in light_form:
            if cid not in labels.index:
                continue
            it = labels.loc[cid, "ion_type"]
            if not _is_null_ion_type(it) and it != type_light:
                continue
            labels.loc[cid, "feature_group_id"] = keep_gid
            labels.loc[cid, "ion_type"] = type_light
            labels.loc[cid, "mono_cluster_id"] = light_mono

        for cid in heavy_form:
            if cid not in labels.index:
                continue
            it = labels.loc[cid, "ion_type"]
            if not _is_null_ion_type(it) and it != type_heavy:
                continue
            labels.loc[cid, "feature_group_id"] = keep_gid
            labels.loc[cid, "ion_type"] = type_heavy
            labels.loc[cid, "mono_cluster_id"] = heavy_mono

    return labels


def validate_feature_group_labels(
    labels: pd.DataFrame,
    mz_by_cluster: Dict[int, float],
    params: FeatureGroupParams,
) -> pd.DataFrame:
    """
    Enforce group invariants; unlabel violator clusters.

    - At most one mono per ion_type within a group.
    - Form monos share neutral mass within ppm.
    - Unique (ion_type, isotope_state) pairs within a group.
    """
    if labels is None or labels.empty:
        return labels
    labels = labels.copy()
    ppm = float(params.mz_tol_ppm)

    grouped = labels.dropna(subset=["feature_group_id"])
    if grouped.empty:
        return labels

    for gid, sub in grouped.groupby("feature_group_id"):
        to_clear = set()

        # Unique (ion_type, isotope_state)
        if sub["ion_type"].notna().any():
            key = sub.apply(
                lambda r: (r["ion_type"], r["isotope_state"]), axis=1
            )
            for k, cnt in key.value_counts().items():
                if cnt > 1 and k[0] is not None and not (isinstance(k[0], float) and np.isnan(k[0])):
                    dups = sub.index[key == k].tolist()
                    # keep lowest cluster id
                    for cid in sorted(dups)[1:]:
                        to_clear.add(int(cid))

        # One mono per ion_type
        for ion_type, form_sub in sub.groupby(sub["ion_type"], dropna=True):
            monos = form_sub.index[form_sub["ion_role"] == "mono"].tolist()
            if len(monos) > 1:
                # keep mono with highest intensity proxy: lowest m/z as mono of form
                # deterministic: keep smallest cluster id
                keep = min(int(c) for c in monos)
                for cid in monos:
                    if int(cid) != keep:
                        to_clear.add(int(cid))
                monos = [keep]

            if len(monos) != 1:
                continue
            mono_c = int(monos[0])
            # isotopes must point at this mono
            for cid, row in form_sub.iterrows():
                if row["ion_role"] == "isotope":
                    mid = row["mono_cluster_id"]
                    if pd.isna(mid) or int(mid) != mono_c:
                        to_clear.add(int(cid))

        # Shared M among form monos still in group
        form_Ms = []
        for ion_type, form_sub in sub.groupby(sub["ion_type"], dropna=True):
            monos = [
                int(c)
                for c in form_sub.index[form_sub["ion_role"] == "mono"].tolist()
                if int(c) not in to_clear
            ]
            if not monos:
                continue
            mono_c = monos[0]
            if mono_c not in mz_by_cluster:
                continue
            try:
                M = neutral_mass_from_mz(mz_by_cluster[mono_c], str(ion_type))
            except Exception:
                continue
            form_Ms.append((str(ion_type), mono_c, M))

        if len(form_Ms) >= 2:
            ref_M = form_Ms[0][2]
            for ion_type, mono_c, M in form_Ms[1:]:
                tol = max(abs(ref_M), abs(M), 1.0) * ppm * 1e-6
                if abs(M - ref_M) > tol:
                    form_mask = sub["ion_type"] == ion_type
                    for cid in sub.index[form_mask]:
                        to_clear.add(int(cid))

        for cid in to_clear:
            if cid not in labels.index:
                continue
            labels.loc[cid, "feature_group_id"] = pd.NA
            labels.loc[cid, "ion_role"] = None
            labels.loc[cid, "ion_type"] = None
            labels.loc[cid, "isotope_state"] = None
            labels.loc[cid, "mono_cluster_id"] = pd.NA

    return labels


def group_features_arrays(
    cluster_ids: np.ndarray,
    mz: np.ndarray,
    rt: np.ndarray,
    heights: np.ndarray,
    params: FeatureGroupParams,
    timings_out: Optional[MutableMapping[str, float]] = None,
) -> pd.DataFrame:
    """
    Run natural-abundance isotope + adduct feature grouping on array inputs.

    1. Unit isotope edges (``Atoms``) + Pearson gate + roll-up.
    2. Base↔adduct edges (``ion_type_dict`` mass offsets) + Pearson gate.
    3. Merge groups across adduct edges; set ``ion_type`` / roles.

    Does not model tracer/labeled enrichment. No mono-vs-family height prior.

    Parameters
    ----------
    cluster_ids : array-like, shape (N,)
    mz, rt : array-like, shape (N,)
    heights : array-like, shape (N, S)
        Apex intensity matrix (not area).
    params : FeatureGroupParams
    timings_out : mutable mapping, optional
        If provided, filled with stage wall times in seconds
        (``isotope_edges``, ``isotope_corr``, ``isotope_labels``,
        ``adduct_edges``, ``adduct_corr``, ``adduct_merge``, ``total``,
        plus edge/type counts).

    Returns
    -------
    DataFrame
        Index = cluster_ids; columns GROUP_COLUMNS.
    """
    t_all = time.perf_counter()
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

    if timings_out is not None:
        timings_out.clear()
        timings_out["n_clusters"] = float(len(cluster_ids))
        timings_out["n_ion_types"] = float(len(params.ion_types))
        timings_out["n_samples"] = float(heights.shape[1])

    if len(cluster_ids) < 2:
        if timings_out is not None:
            timings_out["total"] = time.perf_counter() - t_all
        return empty_group_labels(cluster_ids)

    # Isotope stage
    t0 = time.perf_counter()
    iso_edges = find_isotope_edges(cluster_ids, mz, rt, params)
    if timings_out is not None:
        timings_out["isotope_edges"] = time.perf_counter() - t0
        timings_out["n_isotope_edges_geom"] = float(len(iso_edges))

    t0 = time.perf_counter()
    iso_edges = filter_edges_by_height_correlation(iso_edges, heights, params)
    if timings_out is not None:
        timings_out["isotope_corr"] = time.perf_counter() - t0
        timings_out["n_isotope_edges"] = float(len(iso_edges))

    t0 = time.perf_counter()
    labels = assign_isotope_labels(
        cluster_ids, mz, rt, heights, iso_edges, params
    )
    if timings_out is not None:
        timings_out["isotope_labels"] = time.perf_counter() - t0

    mz_by_cluster = {
        int(cid): float(m) for cid, m in zip(cluster_ids, mz)
    }
    # Adduct stage: Pass1 same-|z| + Pass2 series multi (inside find_adduct_edges)
    if len(params.ion_types) >= 2:
        t0 = time.perf_counter()
        add_edges = find_adduct_edges(cluster_ids, mz, rt, params)
        if timings_out is not None:
            timings_out["adduct_edges"] = time.perf_counter() - t0
            timings_out["n_adduct_edges_geom"] = float(len(add_edges))

        t0 = time.perf_counter()
        add_edges = filter_edges_by_height_correlation(add_edges, heights, params)
        if timings_out is not None:
            timings_out["adduct_corr"] = time.perf_counter() - t0
            timings_out["n_adduct_edges"] = float(len(add_edges))

        t0 = time.perf_counter()
        add_edges = sort_adduct_edges(add_edges)
        if timings_out is not None:
            timings_out["adduct_sort"] = time.perf_counter() - t0

        t0 = time.perf_counter()
        labels = merge_adduct_edges_into_labels(
            labels,
            cluster_ids,
            add_edges,
            params,
            mz_by_cluster=mz_by_cluster,
        )
        if timings_out is not None:
            timings_out["adduct_merge"] = time.perf_counter() - t0
    elif timings_out is not None:
        timings_out["adduct_edges"] = 0.0
        timings_out["adduct_corr"] = 0.0
        timings_out["adduct_sort"] = 0.0
        timings_out["adduct_merge"] = 0.0
        timings_out["n_adduct_edges"] = 0.0

    t0 = time.perf_counter()
    labels = validate_feature_group_labels(labels, mz_by_cluster, params)
    if timings_out is not None:
        timings_out["validate"] = time.perf_counter() - t0
        timings_out["total"] = time.perf_counter() - t_all

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
