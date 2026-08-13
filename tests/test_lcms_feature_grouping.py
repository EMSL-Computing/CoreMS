"""Unit tests for consensus feature grouping (natural-abundance isotopes Stage 1)."""

import numpy as np
import pandas as pd
import pytest

from corems.encapsulation.constant import Atoms
from corems.mass_spectra.calc.feature_grouping import (
    FeatureGroupParams,
    adduct_mass_delta,
    empty_group_labels,
    filter_edges_by_height_correlation,
    filter_ion_types_for_polarity,
    find_adduct_edges,
    find_isotope_edges,
    group_features_arrays,
    ion_type_polarity,
    isotope_mass_delta,
    isotope_state_label,
    normalize_ms_polarity,
    params_with_polarity_filtered_ion_types,
    validate_feature_group_params,
)


def _delta_c13(charge=1):
    return (Atoms.atomic_masses["13C"] - Atoms.atomic_masses["C"]) / abs(charge)


def _delta_nh4_vs_h(charge=1):
    """[M+NH4]+ − [M+H]+ spacing from ion_type_dict / Atoms."""
    return adduct_mass_delta("[M+H]+", "[M+NH4]+", charge=charge)


def test_isotope_mass_delta_uses_atoms_not_hardcoded():
    d = isotope_mass_delta("C", charge=1)
    expected = Atoms.atomic_masses["13C"] - Atoms.atomic_masses["C"]
    assert d == pytest.approx(expected)
    assert d != pytest.approx(1.003355) or d == pytest.approx(
        Atoms.atomic_masses["13C"] - Atoms.atomic_masses["C"]
    )
    d2 = isotope_mass_delta("C", charge=2)
    assert d2 == pytest.approx(expected / 2)


def test_validate_params_rejects_bad_values():
    with pytest.raises(ValueError, match="alignment_rt_tol"):
        validate_feature_group_params(FeatureGroupParams(rt_tol=0))
    with pytest.raises(ValueError, match="alignment_mz_tol_ppm"):
        validate_feature_group_params(FeatureGroupParams(mz_tol_ppm=0))
    with pytest.raises(ValueError, match="charge"):
        validate_feature_group_params(FeatureGroupParams(min_charge=0, max_charge=1))
    with pytest.raises(ValueError, match="charge"):
        validate_feature_group_params(FeatureGroupParams(min_charge=1, max_charge=0))
    with pytest.raises(ValueError, match="feature_group_isotope_atoms"):
        validate_feature_group_params(FeatureGroupParams(isotope_atoms=()))
    with pytest.raises(ValueError, match="feature_group_isotope_atoms|Unknown mono"):
        validate_feature_group_params(FeatureGroupParams(isotope_atoms=("NotAnElement",)))


def test_from_lcms_collection_settings_uses_alignment_tols():
    from corems.encapsulation.factory.processingSetting import LCMSCollectionSettings

    s = LCMSCollectionSettings()
    s.alignment_rt_tol = 0.15
    s.alignment_mz_tol_ppm = 7
    params = FeatureGroupParams.from_lcms_collection_settings(s)
    assert params.rt_tol == pytest.approx(0.15)
    assert params.mz_tol_ppm == pytest.approx(7.0)
    assert "feature_group_rt_tol" not in getattr(s, "__dataclass_fields__", {})
    assert "feature_group_mz_tol_ppm" not in getattr(s, "__dataclass_fields__", {})


def test_default_feature_group_settings_locked_in():
    """Pearson + intensity policy with adduct-ready numeric defaults."""
    from corems.encapsulation.factory.processingSetting import LCMSCollectionSettings
    from corems.mass_spectra.calc.feature_grouping import CORR_METHOD, HEIGHT_COL

    assert CORR_METHOD == "pearson"
    assert HEIGHT_COL == "intensity"

    s = LCMSCollectionSettings()
    assert s.feature_group_corr_threshold == pytest.approx(0.80)
    assert s.feature_group_min_shared_sample_fraction == pytest.approx(0.15)
    assert not hasattr(s, "feature_group_mono_height_fraction")
    assert s.feature_group_max_isotope_offset == 4
    assert s.feature_group_min_isotope_abundance == pytest.approx(0.01)
    assert s.feature_group_min_charge == 1
    assert s.feature_group_max_charge == 1
    assert s.feature_group_ion_types == ("[M+H]+", "[M+NH4]+")

    params = FeatureGroupParams.from_lcms_collection_settings(s)
    assert params.corr_threshold == pytest.approx(0.80)
    assert params.min_shared_sample_fraction == pytest.approx(0.15)
    assert not hasattr(params, "mono_height_fraction")

    # FeatureGroupParams dataclass defaults match settings defaults
    bare = FeatureGroupParams()
    assert bare.corr_threshold == pytest.approx(0.80)
    assert bare.min_shared_sample_fraction == pytest.approx(0.15)


def test_filter_ion_types_for_polarity_drops_wrong_sign():
    """Negative adducts must not be candidates on positive data (and vice versa)."""
    mixed = (
        "[M+H]+",
        "[M+NH4]+",
        "[M+Na]+",
        "[M-H]-",
        "[M+Cl]-",
        "[M+HCOO]-",
        "[M+CH3COO]-",
    )
    assert filter_ion_types_for_polarity(mixed, "positive") == (
        "[M+H]+",
        "[M+NH4]+",
        "[M+Na]+",
    )
    assert filter_ion_types_for_polarity(mixed, "negative") == (
        "[M-H]-",
        "[M+Cl]-",
        "[M+HCOO]-",
        "[M+CH3COO]-",
    )
    assert filter_ion_types_for_polarity(mixed, 1) == (
        "[M+H]+",
        "[M+NH4]+",
        "[M+Na]+",
    )
    assert filter_ion_types_for_polarity(mixed, -1) == (
        "[M-H]-",
        "[M+Cl]-",
        "[M+HCOO]-",
        "[M+CH3COO]-",
    )
    # Unknown polarity: leave list unchanged (array/unit-test path)
    assert filter_ion_types_for_polarity(mixed, None) == mixed
    assert filter_ion_types_for_polarity(mixed, "") == mixed

    assert ion_type_polarity("[M+HCOO]-") == "negative"
    assert ion_type_polarity("[M+H]+") == "positive"
    assert ion_type_polarity("protonated") is None
    assert normalize_ms_polarity("pos") == "positive"
    assert normalize_ms_polarity("neg") == "negative"


def test_from_settings_filters_ion_types_by_polarity():
    from corems.encapsulation.factory.processingSetting import LCMSCollectionSettings

    s = LCMSCollectionSettings()
    s.feature_group_ion_types = (
        "[M+H]+",
        "[M+NH4]+",
        "[M+HCOO]-",
        "[M+CH3COO]-",
    )
    pos = FeatureGroupParams.from_lcms_collection_settings(s, polarity="positive")
    assert pos.ion_types == ("[M+H]+", "[M+NH4]+")
    neg = FeatureGroupParams.from_lcms_collection_settings(s, polarity="negative")
    assert neg.ion_types == ("[M+HCOO]-", "[M+CH3COO]-")
    # Without polarity, configured list is kept (caller responsibility)
    raw = FeatureGroupParams.from_lcms_collection_settings(s)
    assert raw.ion_types == s.feature_group_ion_types

    params = FeatureGroupParams(ion_types=s.feature_group_ion_types)
    filtered = params_with_polarity_filtered_ion_types(params, "positive")
    assert filtered.ion_types == ("[M+H]+", "[M+NH4]+")
    assert params.ion_types == s.feature_group_ion_types  # original frozen params intact


def test_mono_plus_c13_high_corr_groups():
    """Mono, ¹³C, NH₄ adduct, and ¹³C–NH₄ all share one feature_group_id."""
    dm_c = _delta_c13()
    dm_nh4 = _delta_nh4_vs_h()
    # 10: [M+H]+ mono, 11: [M+H]+ 13C1, 13: [M+NH4]+, 14: [M+NH4]+ 13C1, 12: noise
    cluster_ids = np.array([10, 11, 12, 13, 14])
    mz = np.array(
        [
            200.0,
            200.0 + dm_c,
            350.0,
            200.0 + dm_nh4,
            200.0 + dm_nh4 + dm_c,
        ]
    )
    rt = np.array([5.0, 5.02, 5.0, 5.01, 5.03])
    # Shared correlation pattern for true family; noise uncorrelated
    heights = np.array(
        [
            [10.0, 20.0, 30.0, 40.0],
            [5.0, 10.0, 15.0, 20.0],
            [1.0, 50.0, 1.0, 50.0],
            [8.0, 16.0, 24.0, 32.0],
            [4.0, 8.0, 12.0, 16.0],
        ]
    )
    params = FeatureGroupParams(
        rt_tol=0.1,
        mz_tol_ppm=20.0,
        corr_threshold=0.9,
        min_shared_sample_fraction=0.75,
        ion_types=("[M+H]+", "[M+NH4]+"),
    )
    labels = group_features_arrays(cluster_ids, mz, rt, heights, params)

    # [M+H]+ form (lighter offset)
    assert labels.loc[10, "ion_role"] == "mono"
    assert labels.loc[10, "ion_type"] == "[M+H]+"
    assert labels.loc[10, "isotope_state"] == "M+0"
    assert labels.loc[10, "mono_cluster_id"] == 10
    assert labels.loc[11, "ion_role"] == "isotope"
    assert labels.loc[11, "ion_type"] == "[M+H]+"
    assert labels.loc[11, "isotope_state"] == "13C1"
    assert labels.loc[11, "mono_cluster_id"] == 10

    # [M+NH4]+ form (heavier offset): still chemical mono/isotope, not demoted
    assert labels.loc[13, "ion_role"] == "mono"
    assert labels.loc[13, "ion_type"] == "[M+NH4]+"
    assert labels.loc[13, "isotope_state"] == "M+0"
    assert labels.loc[13, "mono_cluster_id"] == 13
    assert labels.loc[14, "ion_role"] == "isotope"
    assert labels.loc[14, "ion_type"] == "[M+NH4]+"
    assert labels.loc[14, "isotope_state"] == "13C1"
    assert labels.loc[14, "mono_cluster_id"] == 13

    gid = labels.loc[10, "feature_group_id"]
    assert pd.notna(gid)
    assert labels.loc[11, "feature_group_id"] == gid
    assert labels.loc[13, "feature_group_id"] == gid
    assert labels.loc[14, "feature_group_id"] == gid
    assert pd.isna(labels.loc[12, "feature_group_id"])


def test_rt_outside_tol_no_group():
    dm = _delta_c13()
    cluster_ids = np.array([1, 2])
    mz = np.array([200.0, 200.0 + dm])
    rt = np.array([5.0, 5.5])  # 0.5 min apart
    heights = np.array([[10.0, 20.0, 30.0], [5.0, 10.0, 15.0]])
    params = FeatureGroupParams(rt_tol=0.1, mz_tol_ppm=20.0, min_shared_sample_fraction=0.5)
    labels = group_features_arrays(cluster_ids, mz, rt, heights, params)
    assert labels["feature_group_id"].isna().all()


def test_low_correlation_no_group():
    dm = _delta_c13()
    cluster_ids = np.array([1, 2])
    mz = np.array([200.0, 200.0 + dm])
    rt = np.array([5.0, 5.01])
    heights = np.array(
        [
            [10.0, 20.0, 30.0, 40.0],
            [40.0, 5.0, 40.0, 5.0],  # anti/low corr on shared samples
        ]
    )
    params = FeatureGroupParams(
        rt_tol=0.1, mz_tol_ppm=20.0, corr_threshold=0.9, min_shared_sample_fraction=0.75
    )
    labels = group_features_arrays(cluster_ids, mz, rt, heights, params)
    assert labels["feature_group_id"].isna().all()


def test_correlation_uses_only_shared_nonzero_samples():
    """Zeros (missing) must not enter Pearson; shared positive pattern can pass."""
    dm = _delta_c13()
    cluster_ids = np.array([1, 2])
    mz = np.array([200.0, 200.0 + dm])
    rt = np.array([5.0, 5.01])
    # Shared samples track well; zeros would drag correlation if included
    heights = np.array(
        [
            [100.0, 80.0, 0.0, 60.0, 0.0],
            [50.0, 40.0, 0.0, 30.0, 0.0],
        ]
    )
    params = FeatureGroupParams(
        rt_tol=0.1,
        mz_tol_ppm=20.0,
        corr_threshold=0.7,
        min_shared_sample_fraction=0.4,  # need 2 of 5
    )
    labels = group_features_arrays(cluster_ids, mz, rt, heights, params)
    assert labels.loc[1, "ion_role"] == "mono"
    assert labels.loc[2, "ion_role"] == "isotope"


def test_chemical_mono_labeled_even_if_shorter_than_13c1():
    """Large-envelope case: M+0 can be shorter than ¹³C₁ and still be mono."""
    dm = _delta_c13()
    cluster_ids = np.array([1, 2])
    mz = np.array([200.0, 200.0 + dm])
    rt = np.array([5.0, 5.01])
    # mono (lighter) smaller; 13C taller — still label via geometry + corr
    heights = np.array(
        [
            [40.0, 32.0, 24.0, 16.0],
            [100.0, 80.0, 60.0, 40.0],
        ]
    )
    params = FeatureGroupParams(
        rt_tol=0.1,
        mz_tol_ppm=20.0,
        corr_threshold=0.9,
        min_shared_sample_fraction=0.75,
    )
    labels = group_features_arrays(cluster_ids, mz, rt, heights, params)
    assert labels.loc[1, "ion_role"] == "mono"
    assert labels.loc[1, "isotope_state"] == "M+0"
    assert labels.loc[2, "ion_role"] == "isotope"
    assert labels.loc[2, "isotope_state"] == "13C1"


def test_fe54_lighter_than_mono_groups():
    """⁵⁴Fe is lighter than ⁵⁶Fe; mono is the higher-m/z most-abundant form."""
    dm = abs(Atoms.atomic_masses["54Fe"] - Atoms.atomic_masses["Fe"])
    # cluster 1 = 54Fe (lighter, smaller), cluster 2 = Fe mono (heavier, taller)
    cluster_ids = np.array([1, 2])
    mz = np.array([400.0, 400.0 + dm])
    rt = np.array([5.0, 5.01])
    heights = np.array(
        [
            [20.0, 16.0, 12.0, 8.0],
            [100.0, 80.0, 60.0, 40.0],
        ]
    )
    params = FeatureGroupParams(
        isotope_atoms=("Fe",),
        rt_tol=0.1,
        mz_tol_ppm=20.0,
        corr_threshold=0.9,
        min_shared_sample_fraction=0.75,
    )
    labels = group_features_arrays(cluster_ids, mz, rt, heights, params)
    assert labels.loc[2, "ion_role"] == "mono"
    assert labels.loc[2, "isotope_state"] == "M+0"
    assert labels.loc[1, "ion_role"] == "isotope"
    assert labels.loc[1, "isotope_state"] == "54Fe1"
    assert labels.loc[1, "mono_cluster_id"] == 2


def test_c13_chain_13c1_13c2():
    """¹³C₂ only via roll-up through ¹³C₁ (unit edges), not a direct mono→M+2 link."""
    dm = _delta_c13()
    cluster_ids = np.array([1, 2, 3])
    mz = np.array([200.0, 200.0 + dm, 200.0 + 2 * dm])
    rt = np.array([5.0, 5.01, 5.02])
    heights = np.array(
        [
            [100.0, 80.0, 60.0, 40.0],
            [50.0, 40.0, 30.0, 20.0],
            [20.0, 16.0, 12.0, 8.0],
        ]
    )
    params = FeatureGroupParams(
        rt_tol=0.1,
        mz_tol_ppm=20.0,
        corr_threshold=0.9,
        min_shared_sample_fraction=0.75,
        max_isotope_offset=4,
    )
    labels = group_features_arrays(cluster_ids, mz, rt, heights, params)
    assert labels.loc[1, "ion_role"] == "mono"
    assert labels.loc[2, "isotope_state"] == "13C1"
    assert labels.loc[3, "isotope_state"] == "13C2"
    assert labels.loc[2, "mono_cluster_id"] == 1
    assert labels.loc[3, "mono_cluster_id"] == 1
    assert labels.loc[1, "feature_group_id"] == labels.loc[3, "feature_group_id"]

    # Unit edges only: mono–M+1 and M+1–M+2, no mono–M+2
    edges = find_isotope_edges(
        cluster_ids, mz, rt, params
    )
    assert (edges["n"] == 1).all()
    pairs = set(zip(edges["parent_cluster"], edges["child_cluster"]))
    assert (1, 2) in pairs
    assert (2, 3) in pairs
    assert (1, 3) not in pairs


def test_no_13c2_without_13c1():
    """Direct 2×Δm without intermediate ¹³C₁ must not label as ¹³C₂."""
    dm = _delta_c13()
    cluster_ids = np.array([1, 2])
    mz = np.array([200.0, 200.0 + 2 * dm])
    rt = np.array([5.0, 5.01])
    heights = np.array(
        [
            [100.0, 80.0, 60.0, 40.0],
            [20.0, 16.0, 12.0, 8.0],
        ]
    )
    params = FeatureGroupParams(
        rt_tol=0.1,
        mz_tol_ppm=20.0,
        corr_threshold=0.9,
        min_shared_sample_fraction=0.75,
        max_isotope_offset=4,
    )
    labels = group_features_arrays(cluster_ids, mz, rt, heights, params)
    assert labels["feature_group_id"].isna().all()


def test_singleton_unlabeled():
    labels = group_features_arrays(
        np.array([7]),
        np.array([100.0]),
        np.array([1.0]),
        np.array([[1.0, 2.0]]),
        FeatureGroupParams(),
    )
    assert pd.isna(labels.loc[7, "feature_group_id"])
    assert labels.loc[7, "ion_role"] is None or pd.isna(labels.loc[7, "ion_role"])


def test_charge_two_spacing():
    dm = _delta_c13(charge=2)
    cluster_ids = np.array([1, 2])
    mz = np.array([400.0, 400.0 + dm])
    rt = np.array([5.0, 5.01])
    heights = np.array([[10.0, 20.0, 30.0], [5.0, 10.0, 15.0]])
    params = FeatureGroupParams(
        min_charge=2,
        max_charge=2,
        rt_tol=0.1,
        mz_tol_ppm=20.0,
        min_shared_sample_fraction=0.5,
        corr_threshold=0.9,
    )
    labels = group_features_arrays(cluster_ids, mz, rt, heights, params)
    assert labels.loc[1, "ion_role"] == "mono"
    assert labels.loc[2, "ion_role"] == "isotope"


def test_min_max_charge_range_matches_z2_when_enabled():
    """With min=1 max=3, a pure |z|=2 spacing still groups."""
    dm = _delta_c13(charge=2)
    cluster_ids = np.array([1, 2])
    mz = np.array([400.0, 400.0 + dm])
    rt = np.array([5.0, 5.01])
    heights = np.array([[10.0, 20.0, 30.0], [5.0, 10.0, 15.0]])
    params = FeatureGroupParams(
        min_charge=1,
        max_charge=3,
        rt_tol=0.1,
        mz_tol_ppm=20.0,
        min_shared_sample_fraction=0.5,
        corr_threshold=0.9,
    )
    labels = group_features_arrays(cluster_ids, mz, rt, heights, params)
    assert labels.loc[1, "ion_role"] == "mono"
    assert labels.loc[2, "ion_role"] == "isotope"


def test_from_settings_min_max_charge():
    from corems.encapsulation.factory.processingSetting import LCMSCollectionSettings

    s = LCMSCollectionSettings()
    s.feature_group_min_charge = 1
    s.feature_group_max_charge = 3
    params = FeatureGroupParams.from_lcms_collection_settings(s)
    assert params.min_charge == 1
    assert params.max_charge == 3
    assert params.charge_values() == (1, 2, 3)


def test_rerun_clears_via_empty_then_assign():
    dm = _delta_c13()
    cluster_ids = np.array([1, 2])
    mz = np.array([200.0, 200.0 + dm])
    rt = np.array([5.0, 5.01])
    heights = np.array([[10.0, 20.0, 30.0], [5.0, 10.0, 15.0]])
    params = FeatureGroupParams(
        rt_tol=0.1, mz_tol_ppm=20.0, min_shared_sample_fraction=0.5, corr_threshold=0.9
    )
    labels1 = group_features_arrays(cluster_ids, mz, rt, heights, params)
    assert labels1["feature_group_id"].notna().any()
    # second run with impossible RT tol yields empty
    params2 = FeatureGroupParams(
        rt_tol=1e-9, mz_tol_ppm=20.0, min_shared_sample_fraction=0.5, corr_threshold=0.9
    )
    labels2 = group_features_arrays(cluster_ids, mz, rt, heights, params2)
    assert labels2["feature_group_id"].isna().all()


def test_isotope_state_label():
    assert isotope_state_label("13C", 0) == "M+0"
    assert isotope_state_label("13C", 1) == "13C1"
    assert isotope_state_label("13C", 2) == "13C2"
    assert isotope_state_label("54Fe", 1) == "54Fe1"


def test_find_edges_uses_atoms_spacing():
    dm = _delta_c13()
    cluster_ids = np.array([1, 2])
    mz = np.array([100.0, 100.0 + dm])
    rt = np.array([1.0, 1.0])
    params = FeatureGroupParams(rt_tol=0.1, mz_tol_ppm=10.0)
    edges = find_isotope_edges(cluster_ids, mz, rt, params)
    assert len(edges) == 1
    assert edges.iloc[0]["n"] == 1
    assert edges.iloc[0]["atom"] == "C"


def test_rare_isotope_entries_filters_by_abundance():
    from corems.mass_spectra.calc.feature_grouping import rare_isotope_entries

    # 74Se abundance is ~0.0089; excluded at default 0.01, included at 0.0
    se_default = rare_isotope_entries("Se", min_abundance=0.01)
    labels = {r[0] for r in se_default}
    assert "74Se" not in labels
    assert "78Se" in labels
    assert "76Se" in labels
    assert "77Se" in labels
    assert "82Se" in labels

    se_all = rare_isotope_entries("Se", min_abundance=0.0)
    labels_all = {r[0] for r in se_all}
    assert "74Se" in labels_all

    # Ordered by decreasing abundance
    abs_order = [r[2] for r in se_default]
    assert abs_order == sorted(abs_order, reverse=True)


def test_se_multiple_isotopes_create_multiple_edges():
    """Se mono + several rare forms above abundance floor → multiple edges."""
    from corems.mass_spectra.calc.feature_grouping import rare_isotope_entries

    mono_mz = 500.0
    entries = rare_isotope_entries("Se", min_abundance=0.01)
    # mono + each rare as a separate cluster
    cluster_ids = [0] + list(range(1, len(entries) + 1))
    mz = [mono_mz] + [
        mono_mz + (Atoms.atomic_masses[lab] - Atoms.atomic_masses["Se"])
        for lab, _, _ in entries
    ]
    rt = [5.0] * len(cluster_ids)
    params = FeatureGroupParams(
        isotope_atoms=("Se",),
        min_isotope_abundance=0.01,
        rt_tol=0.1,
        mz_tol_ppm=20.0,
        max_isotope_offset=1,
    )
    edges = find_isotope_edges(
        np.array(cluster_ids), np.array(mz), np.array(rt), params
    )
    expected = {lab for lab, _, _ in entries}
    # Edges with chemical mono (cluster 0) as parent cover each rare form
    from_mono = edges[edges["parent_cluster"] == 0]
    rares_from_mono = set(from_mono["rare_label"].tolist())
    assert expected.issubset(rares_from_mono)
    assert len(from_mono) >= len(expected)
