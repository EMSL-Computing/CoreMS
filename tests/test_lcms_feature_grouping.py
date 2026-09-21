"""Unit tests for consensus feature grouping (natural-abundance isotopes Stage 1)."""

import numpy as np
import pandas as pd
import pytest

from corems.encapsulation.constant import Atoms
from corems.mass_spectra.calc.feature_grouping import (
    DEFAULT_ION_TYPES,
    FeatureGroupParams,
    adduct_mass_delta,
    empty_group_labels,
    filter_edges_by_height_correlation,
    filter_ion_types_for_polarity,
    find_adduct_edges,
    find_isotope_edges,
    group_features_arrays,
    ion_type_charge,
    is_allowed_adduct_type_pair,
    isotope_mass_delta,
    isotope_state_label,
    neutral_mass_from_mz,
    normalize_ms_polarity,
    params_with_polarity_filtered_ion_types,
    validate_feature_group_params,
)

# Default ordered set (most → least common) — must resolve in ion_type_dict.
COMMON_FEATURE_GROUP_ION_TYPES = DEFAULT_ION_TYPES


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
    with pytest.raises(ValueError, match="feature_group_isotope_atoms"):
        validate_feature_group_params(FeatureGroupParams(isotope_atoms=()))
    with pytest.raises(ValueError, match="feature_group_isotope_atoms|Unknown mono"):
        validate_feature_group_params(FeatureGroupParams(isotope_atoms=("NotAnElement",)))


def test_defaults_are_singly_charged_only():
    from corems.encapsulation.factory.processingSetting import LCMSCollectionSettings

    assert not hasattr(LCMSCollectionSettings(), "feature_group_min_charge")
    assert not hasattr(LCMSCollectionSettings(), "feature_group_max_charge")
    assert not hasattr(FeatureGroupParams(), "min_charge")
    assert not hasattr(FeatureGroupParams(), "max_charge")
    for it in DEFAULT_ION_TYPES:
        assert ion_type_charge(it) == 1, it
    params = FeatureGroupParams.from_lcms_collection_settings(LCMSCollectionSettings())
    assert params.ion_types == DEFAULT_ION_TYPES
    assert not hasattr(params, "charge_values")


def test_validate_rejects_multicharge_ion_type():
    with pytest.raises(ValueError, match=r"singly-charged|charge"):
        validate_feature_group_params(
            FeatureGroupParams(ion_types=("[M+H]+", "[M+2H]2+"))
        )


def test_c13_chain_13c1_13c2_13c3():
    dm = _delta_c13()
    cluster_ids = np.array([1, 2, 3, 4])
    mz = np.array([200.0, 200.0 + dm, 200.0 + 2 * dm, 200.0 + 3 * dm])
    rt = np.array([5.0, 5.01, 5.02, 5.03])
    pat = np.array([100.0, 80.0, 60.0, 40.0])
    heights = np.vstack([pat, pat * 0.5, pat * 0.2, pat * 0.08])
    params = FeatureGroupParams(
        rt_tol=0.1,
        mz_tol_ppm=20.0,
        corr_threshold=0.9,
        min_shared_sample_fraction=0.75,
        max_isotope_offset=4,
        ion_types=("[M+H]+",),
    )
    labels = group_features_arrays(cluster_ids, mz, rt, heights, params)
    assert labels.loc[1, "ion_role"] == "mono"
    assert labels.loc[2, "isotope_state"] == "13C1"
    assert labels.loc[3, "isotope_state"] == "13C2"
    assert labels.loc[4, "isotope_state"] == "13C3"
    gid = labels.loc[1, "feature_group_id"]
    assert pd.notna(gid)
    assert (labels.loc[[2, 3, 4], "feature_group_id"] == gid).all()
    assert (labels.loc[[2, 3, 4], "mono_cluster_id"] == 1).all()

    edges = find_isotope_edges(cluster_ids, mz, rt, params)
    assert (edges["n"] == 1).all()
    assert (edges["charge"] == 1).all()
    pairs = set(zip(edges["parent_cluster"], edges["child_cluster"]))
    assert (1, 2) in pairs and (2, 3) in pairs and (3, 4) in pairs
    assert (1, 3) not in pairs and (1, 4) not in pairs


def test_z2_c13_spacing_does_not_group():
    """Δm = ¹³C/2 must not be treated as a unit isotope edge."""
    dm = _delta_c13(charge=2)
    cluster_ids = np.array([1, 2])
    mz = np.array([400.0, 400.0 + dm])
    rt = np.array([5.0, 5.01])
    heights = np.array([[10.0, 20.0, 30.0], [5.0, 10.0, 15.0]])
    params = FeatureGroupParams(
        rt_tol=0.1,
        mz_tol_ppm=20.0,
        min_shared_sample_fraction=0.5,
        corr_threshold=0.9,
        ion_types=("[M+H]+",),
    )
    labels = group_features_arrays(cluster_ids, mz, rt, heights, params)
    assert labels["feature_group_id"].isna().all()
    edges = find_isotope_edges(cluster_ids, mz, rt, params)
    assert edges.empty


def test_mh_mna_group_without_isotopes():
    """[M+H]+ and [M+Na]+ (z=1 only) share one feature_group_id."""
    from corems.mass_spectra.calc.feature_grouping import _ion_type_mass_offset

    M = 400.0
    mz_mh = M + _ion_type_mass_offset("[M+H]+")
    mz_na = M + _ion_type_mass_offset("[M+Na]+")
    cluster_ids = np.array([0, 1])
    mz = np.array([mz_mh, mz_na])
    rt = np.array([10.0, 10.01])
    pat = np.array([10.0, 20.0, 30.0, 40.0, 50.0])
    heights = np.vstack([pat, pat * 0.6])
    params = FeatureGroupParams(
        rt_tol=0.1,
        mz_tol_ppm=5.0,
        corr_threshold=0.9,
        min_shared_sample_fraction=0.6,
        ion_types=("[M+H]+", "[M+Na]+"),
    )
    labels = group_features_arrays(cluster_ids, mz, rt, heights, params)
    assert labels.loc[0, "ion_type"] == "[M+H]+"
    assert labels.loc[1, "ion_type"] == "[M+Na]+"
    assert labels.loc[0, "feature_group_id"] == labels.loc[1, "feature_group_id"]
    assert pd.notna(labels.loc[0, "feature_group_id"])
    # No isotope observed → not claimed as chemical mono
    assert pd.isna(labels.loc[0, "mono_cluster_id"])
    assert pd.isna(labels.loc[1, "mono_cluster_id"])
    assert labels.loc[0, "ion_role"] is None or pd.isna(labels.loc[0, "ion_role"])
    assert labels.loc[1, "ion_role"] is None or pd.isna(labels.loc[1, "ion_role"])


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
    assert "[M+2H]2+" not in s.feature_group_ion_types
    assert s.feature_group_ion_types == DEFAULT_ION_TYPES
    assert "[M+H]+" in s.feature_group_ion_types
    # Allow-list can be narrowed in settings
    s.feature_group_ion_types = ("[M+H]+", "[M+Na]+")
    assert FeatureGroupParams.from_lcms_collection_settings(s).ion_types == (
        "[M+H]+",
        "[M+Na]+",
    )
    s.feature_group_ion_types = DEFAULT_ION_TYPES

    params = FeatureGroupParams.from_lcms_collection_settings(s)
    assert params.corr_threshold == pytest.approx(0.80)
    assert params.min_shared_sample_fraction == pytest.approx(0.15)
    assert not hasattr(params, "mono_height_fraction")
    assert params.ion_types == DEFAULT_ION_TYPES

    # FeatureGroupParams dataclass defaults match settings defaults
    bare = FeatureGroupParams()
    assert bare.corr_threshold == pytest.approx(0.80)
    assert bare.min_shared_sample_fraction == pytest.approx(0.15)
    assert bare.ion_types == DEFAULT_ION_TYPES


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

    from corems.encapsulation.constant import ION_TYPE_DICT

    assert ION_TYPE_DICT["[M+HCOO]-"]["polarity"] == "negative"
    assert ION_TYPE_DICT["[M+H]+"]["polarity"] == "positive"
    assert ION_TYPE_DICT["[M+2H]2+"]["polarity"] == "positive"
    assert ION_TYPE_DICT["[M-2H]2-"]["polarity"] == "negative"
    assert ION_TYPE_DICT["protonated"]["polarity"] == "positive"
    assert normalize_ms_polarity("pos") == "positive"
    assert normalize_ms_polarity("neg") == "negative"


def test_common_ion_types_in_dict_and_charge_parse():
    """Allow-list defaults are |z|=1 keys in ION_TYPE_DICT with polarity."""
    from corems.encapsulation.constant import ION_TYPE_DICT
    from corems.mass_spectra.output.export import ion_type_dict

    assert ion_type_dict is ION_TYPE_DICT
    assert len(DEFAULT_ION_TYPES) == 14
    for it in DEFAULT_ION_TYPES:
        assert it in ION_TYPE_DICT, f"missing ION_TYPE_DICT key: {it}"
        entry = ION_TYPE_DICT[it]
        assert entry["polarity"] in ("positive", "negative")
        assert "add" in entry and "sub" in entry
        assert "feature_group_order" not in entry
        assert ion_type_charge(it) == 1
        validate_feature_group_params(
            FeatureGroupParams(ion_types=(it, "[M+H]+") if it != "[M+H]+" else (it, "[M+Na]+"))
        )
    # Narrow allow-list is valid
    validate_feature_group_params(
        FeatureGroupParams(ion_types=("[M+H]+", "[M+Na]+"))
    )

    # Parser still understands multi-charge suffixes (formula/export keys)
    assert ion_type_charge("[M+H]+") == 1
    assert ion_type_charge("[M+2H]2+") == 2
    assert ion_type_charge("[M+3H]3+") == 3
    assert ion_type_charge("[M-2H]2-") == 2
    assert ion_type_charge("[M]+") == 1


def test_water_loss_preferred_over_water_adduct_on_delta_tie():
    """Exact 18.01 Da spacing: both interps kept in possible_ion_types.

    Edge parent/child_ion_type is a geometry handle only; labels do not
    pick a preferred type among residual ties.
    """
    M = 400.0
    from corems.mass_spectra.calc.feature_grouping import _ion_type_mass_offset

    mz_loss = M + _ion_type_mass_offset("[M+H-H2O]+")
    mz_mh = M + _ion_type_mass_offset("[M+H]+")
    cluster_ids = np.array([0, 1])
    mz = np.array([mz_loss, mz_mh])
    rt = np.array([5.0, 5.01])
    params = FeatureGroupParams(
        rt_tol=0.1,
        mz_tol_ppm=5.0,
        ion_types=DEFAULT_ION_TYPES,
    )
    edges = find_adduct_edges(cluster_ids, mz, rt, params)
    assert len(edges) >= 1
    # Lower m/z is parent
    row = edges.iloc[0]
    assert row["parent_ion_type"] == "[M+H-H2O]+"
    assert row["child_ion_type"] == "[M+H]+"
    # Water-adduct interpretation also fits residual → kept as alternate
    parent_poss = str(row["parent_possible_ion_types"])
    child_poss = str(row["child_possible_ion_types"])
    assert "[M+H-H2O]+" in parent_poss
    assert "[M+H]+" in parent_poss  # alternate light type for adduct interp
    assert "[M+H]+" in child_poss
    assert "[M+H+H2O]+" in child_poss
    assert int(row["n_interpretations"]) >= 2


def test_validate_clears_orphan_mono_without_feature_group():
    """ion_role=mono with no feature_group_id must not survive validation."""
    from corems.mass_spectra.calc.feature_grouping import (
        empty_group_labels,
        validate_feature_group_labels,
    )

    labels = empty_group_labels([518, 3352])
    labels.loc[518, "ion_role"] = "mono"
    labels.loc[518, "isotope_state"] = "M+0"
    labels.loc[3352, "ion_role"] = "mono"
    labels.loc[3352, "isotope_state"] = "M+0"
    out = validate_feature_group_labels(
        labels,
        {518: 416.2, 3352: 957.8},
        FeatureGroupParams(),
    )
    assert out.loc[518, "ion_role"] is None or pd.isna(out.loc[518, "ion_role"])
    assert out.loc[3352, "ion_role"] is None or pd.isna(out.loc[3352, "ion_role"])
    assert pd.isna(out.loc[518, "feature_group_id"])


def test_nh3_vs_nh4_keeps_ambiguous_possible_ion_types():
    """Δm = m(NH3): both NH3-loss/MH and MH/NH4 fit; keep both as possibles."""
    from corems.mass_spectra.calc.feature_grouping import _ion_type_mass_offset

    M = 583.5896
    mz_nh3 = M + _ion_type_mass_offset("[M+H-NH3]+")
    mz_mh = M + _ion_type_mass_offset("[M+H]+")
    cluster_ids = np.array([144, 179])
    mz = np.array([mz_nh3, mz_mh])
    rt = np.array([48.71, 48.72])
    pat = np.array([10.0, 20.0, 30.0, 40.0])
    heights = np.vstack([pat, pat * 0.5])
    params = FeatureGroupParams(
        rt_tol=0.1,
        mz_tol_ppm=5.0,
        corr_threshold=0.9,
        min_shared_sample_fraction=0.75,
        ion_types=DEFAULT_ION_TYPES,
    )
    edges = find_adduct_edges(cluster_ids, mz, rt, params)
    assert len(edges) == 1
    row = edges.iloc[0]
    # Preferred among residual ties (alphabetical): NH3-loss + MH
    assert row["parent_ion_type"] == "[M+H-NH3]+"
    assert row["child_ion_type"] == "[M+H]+"
    assert "[M+H]+" in str(row["parent_possible_ion_types"])
    assert "[M+H-NH3]+" in str(row["parent_possible_ion_types"])
    assert "[M+NH4]+" in str(row["child_possible_ion_types"])
    assert "[M+H]+" in str(row["child_possible_ion_types"])
    assert int(row["n_interpretations"]) >= 2

    labels = group_features_arrays(cluster_ids, mz, rt, heights, params)
    assert pd.isna(labels.loc[144, "ion_type"])
    assert pd.isna(labels.loc[179, "ion_type"])
    assert "[M+H]+" in str(labels.loc[144, "possible_ion_types"])
    assert "[M+H-NH3]+" in str(labels.loc[144, "possible_ion_types"])
    assert "[M+NH4]+" in str(labels.loc[179, "possible_ion_types"])
    assert "[M+H]+" in str(labels.loc[179, "possible_ion_types"])
    assert labels.loc[144, "feature_group_id"] == labels.loc[179, "feature_group_id"]


def test_unrelated_high_mz_not_grouped():
    """H/Na family must not glue to an uncorrelated high-m/z stranger."""
    M = 390.2766
    h = Atoms.atomic_masses["H"]
    na = Atoms.atomic_masses["Na"]
    dm_c = _delta_c13()
    mz_mh = M + h
    mz_na = M + na
    mz_hi = 797.5808 + h
    cluster_ids = np.array([0, 1, 2, 3])
    mz = np.array([mz_mh, mz_mh + dm_c, mz_na, mz_hi])
    rt = np.array([29.31, 29.31, 29.39, 29.33])
    pat = np.array([10.0, 20.0, 30.0, 40.0])
    heights = np.vstack([pat, pat * 0.5, pat * 0.7, pat * 0.3])
    params = FeatureGroupParams(
        rt_tol=0.5,
        mz_tol_ppm=5.0,
        corr_threshold=0.9,
        min_shared_sample_fraction=0.75,
        ion_types=DEFAULT_ION_TYPES,
    )
    labels = group_features_arrays(cluster_ids, mz, rt, heights, params)
    gid = labels.loc[0, "feature_group_id"]
    assert pd.notna(gid)
    assert labels.loc[1, "feature_group_id"] == gid
    assert labels.loc[2, "feature_group_id"] == gid
    assert "[M+H]+" in str(labels.loc[0, "possible_ion_types"])
    assert "[M+Na]+" in str(labels.loc[2, "possible_ion_types"])
    assert pd.isna(labels.loc[3, "feature_group_id"]) or labels.loc[
        3, "feature_group_id"
    ] != gid


def test_form_paint_does_not_retype_other_forms():
    """Second adduct edge must not overwrite first form's ion_type on isotopes."""
    dm_c = _delta_c13()
    dm_nh4 = _delta_nh4_vs_h()
    dm_na = Atoms.atomic_masses["Na"] - Atoms.atomic_masses["H"]
    # 0 MH, 1 MH-13C, 2 NH4, 3 Na — all same M base 200
    base = 200.0
    cluster_ids = np.array([0, 1, 2, 3])
    mz = np.array(
        [base, base + dm_c, base + dm_nh4, base + dm_na]
    )
    rt = np.array([5.0, 5.01, 5.02, 5.01])
    pat = np.array([10.0, 20.0, 30.0, 40.0])
    heights = np.vstack([pat, pat * 0.5, pat * 0.8, pat * 0.6])
    params = FeatureGroupParams(
        rt_tol=0.1,
        mz_tol_ppm=20.0,
        corr_threshold=0.9,
        min_shared_sample_fraction=0.75,
        ion_types=("[M+H]+", "[M+NH4]+", "[M+Na]+"),
    )
    labels = group_features_arrays(cluster_ids, mz, rt, heights, params)
    assert labels.loc[0, "ion_type"] == "[M+H]+"
    assert labels.loc[1, "ion_type"] == "[M+H]+"
    assert labels.loc[1, "ion_role"] == "isotope"
    # NH4 and Na may each group with H; isotopes of H stay [M+H]+
    assert labels.loc[1, "ion_type"] != "[M+Na]+"
    assert labels.loc[1, "ion_type"] != "[M+NH4]+"


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


def test_annotation_eligibility_helpers():
    from corems.mass_spectra.calc.feature_grouping import (
        allowed_ion_types_for_row,
        constrain_annotation_active,
        empty_group_labels,
        ion_type_allowed,
        is_adduct_endpoint_eligible,
        should_skip_isotope_for_annotation,
    )

    assert is_adduct_endpoint_eligible(None)
    assert is_adduct_endpoint_eligible("mono")
    assert not is_adduct_endpoint_eligible("isotope")

    assert not constrain_annotation_active(None, True)
    assert not constrain_annotation_active(empty_group_labels([]), True)
    labs = empty_group_labels([1])
    assert constrain_annotation_active(labs, True)
    assert not constrain_annotation_active(labs, False)

    assert should_skip_isotope_for_annotation("isotope", True)
    assert not should_skip_isotope_for_annotation("mono", True)
    assert not should_skip_isotope_for_annotation("isotope", False)

    assert allowed_ion_types_for_row(None, None) is None
    assert allowed_ion_types_for_row("[M+H]+", None) == {"[M+H]+"}
    assert allowed_ion_types_for_row("[M+H]+", "[M+H]+;[M+Na]+") == {
        "[M+H]+",
        "[M+Na]+",
    }
    assert ion_type_allowed("[M+Na]+", {"[M+H]+", "[M+Na]+"})
    assert ion_type_allowed("[m+h]+", {"[M+H]+"})
    assert not ion_type_allowed("[M+K]+", {"[M+H]+"})
    assert ion_type_allowed("[M+K]+", None)
    assert not ion_type_allowed(None, {"[M+H]+"})
    # FlashEntropy stores parallel hit lists on SpectrumSearchResults
    assert ion_type_allowed(["[M+H]+", "[M+K]+"], {"[M+H]+"})
    assert not ion_type_allowed(["[M+K]+"], {"[M+H]+"})
    assert ion_type_allowed(np.array(["[M+Na]+"]), {"[M+Na]+"})

    from corems.mass_spectra.calc.feature_grouping import (
        subset_hits_by_allowed_ion_types,
    )
    from types import SimpleNamespace

    hits = SimpleNamespace(
        ref_ion_type=["[M+H]+", "[M+K]+", "[M+Na]+"],
        entropy_similarity=np.array([0.9, 0.8, 0.7]),
        ref_mol_id=["a", "b", "c"],
        precursor_mz=200.0,
    )
    kept = subset_hits_by_allowed_ion_types(hits, {"[M+H]+", "[M+Na]+"})
    assert kept is hits
    assert kept.ref_ion_type == ["[M+H]+", "[M+Na]+"]
    assert list(kept.entropy_similarity) == [0.9, 0.7]
    assert kept.ref_mol_id == ["a", "c"]
    assert kept.precursor_mz == 200.0
    assert subset_hits_by_allowed_ion_types(hits, {"[M+Li]+"}) is None


def test_isotope_not_used_as_adduct_endpoint():
    """After isotope labeling, adduct edges must not use ion_role=isotope endpoints."""
    from corems.mass_spectra.calc.feature_grouping import _ion_type_mass_offset

    M = 400.0
    dm_c = _delta_c13(charge=1)
    # 0: MH mono, 1: MH 13C1, 2: MNa mono — isotope must not pair as adduct alone
    cluster_ids = np.array([0, 1, 2])
    mz = np.array(
        [
            M + _ion_type_mass_offset("[M+H]+"),
            M + _ion_type_mass_offset("[M+H]+") + dm_c,
            M + _ion_type_mass_offset("[M+Na]+"),
        ]
    )
    rt = np.array([5.0, 5.01, 5.02])
    pat = np.array([10.0, 20.0, 30.0, 40.0])
    heights = np.vstack([pat, pat * 0.5, pat * 0.7])
    params = FeatureGroupParams(
        rt_tol=0.1,
        mz_tol_ppm=5.0,
        corr_threshold=0.9,
        min_shared_sample_fraction=0.6,
        ion_types=("[M+H]+", "[M+Na]+"),
    )
    labels = group_features_arrays(cluster_ids, mz, rt, heights, params)
    assert labels.loc[0, "ion_role"] == "mono"
    assert labels.loc[1, "ion_role"] == "isotope"
    assert labels.loc[2, "ion_type"] == "[M+Na]+"
    assert labels.loc[0, "feature_group_id"] == labels.loc[2, "feature_group_id"]
    assert labels.loc[1, "feature_group_id"] == labels.loc[0, "feature_group_id"]
    # Isotope should not be typed as a separate adduct form
    assert labels.loc[1, "ion_type"] == "[M+H]+"
    # [M+H]+ form has 13C → mono_cluster_id is set; Na form has no isotope
    assert labels.loc[0, "mono_cluster_id"] == 0
    assert labels.loc[1, "mono_cluster_id"] == 0
    assert pd.isna(labels.loc[2, "mono_cluster_id"])
    assert labels.loc[2, "ion_role"] is None or pd.isna(labels.loc[2, "ion_role"])


def test_feature_group_constrain_annotation_setting_default():
    from corems.encapsulation.factory.processingSetting import LCMSCollectionSettings

    s = LCMSCollectionSettings()
    assert s.feature_group_constrain_annotation is True
