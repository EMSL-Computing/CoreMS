"""Unit tests for SpectralSimilaritySearchSettings nested on MSParameters."""

import dataclasses
import textwrap

import pytest

from corems.encapsulation.factory.parameters import (
    LCMSParameters,
    MSParameters,
    apply_spectral_similarity_search_to_collection,
    reset_lcms_parameters,
    reset_ms_parameters,
    settings_from_lcms,
    settings_from_lcms_collection,
)
from corems.encapsulation.factory.processingSetting import (
    LEGACY_PARAMETER,
    LiquidChromatographSetting,
    SpectralSimilaritySearchSettings,
    legacy_field,
    settings_to_export_dict,
)
from corems.encapsulation.input.parameter_from_json import (
    load_and_set_toml_parameters_lcms,
)
from corems.encapsulation.output import parameter_to_dict
from corems.encapsulation.output.parameter_to_json import dump_lcms_settings_toml
from corems.molecular_id.search.database_interfaces import _resolve_fe_kwargs
from corems.molecular_networking.network_builder import MolecularNetwork

_LEGACY_ANNOTATION_KEYS = (
    "ms2_min_fe_score",
    "search_as_lipids",
    "include_fragment_types",
)


class _FakeLCMS:
    def __init__(self):
        self.parameters = LCMSParameters(use_defaults=True)


def test_spectral_similarity_search_settings_defaults():
    s = SpectralSimilaritySearchSettings()
    assert s.max_ms2_tolerance_in_da == 0.01
    assert s.min_ms2_difference_in_da == 0.02
    assert s.ms2_min_fe_score == 0.2
    assert s.search_type == "open"
    assert s.resolved_peak_sep_da == 0.02


def test_as_fe_kwargs_and_derived_min():
    s = SpectralSimilaritySearchSettings(max_ms2_tolerance_in_da=0.05)
    assert s.min_ms2_difference_in_da == 0.1
    kw = s.as_fe_kwargs()
    assert kw["max_ms2_tolerance_in_da"] == 0.05
    assert kw["min_ms2_difference_in_da"] == 0.1
    assert "min_ms2_difference_in_da" not in {
        f.name for f in dataclasses.fields(SpectralSimilaritySearchSettings)
    }


def test_ms_parameters_has_spectral_similarity_search():
    reset_ms_parameters()
    p = MSParameters(use_defaults=True)
    assert isinstance(p.spectral_similarity_search, SpectralSimilaritySearchSettings)
    p2 = p.copy()
    p2.spectral_similarity_search.ms2_min_fe_score = 0.9
    assert p.spectral_similarity_search.ms2_min_fe_score == 0.2


def test_lcms_nested_path_and_no_top_level_spectral_similarity_search():
    reset_lcms_parameters()
    p = LCMSParameters(use_defaults=True)
    assert not hasattr(LCMSParameters, "spectral_similarity_search")
    assert "spectral_similarity_search" not in p.__dict__
    assert isinstance(p.mass_spectrum["ms2"].spectral_similarity_search, SpectralSimilaritySearchSettings)


def test_settings_from_lcms_helpers():
    reset_lcms_parameters()
    obj = _FakeLCMS()
    obj.parameters.mass_spectrum["ms2"].spectral_similarity_search.ms2_min_fe_score = 0.44
    assert settings_from_lcms(obj).ms2_min_fe_score == 0.44

    class _Coll(list):
        pass

    coll = _Coll([obj, obj])
    assert settings_from_lcms_collection(coll).ms2_min_fe_score == 0.44


def test_apply_spectral_similarity_search_to_collection_broadcasts_copies():
    """Mutating settings_from_lcms_collection only hits sample 0; apply fixes all."""
    reset_lcms_parameters()
    a, b = _FakeLCMS(), _FakeLCMS()
    # independent parameter trees
    assert a.parameters is not b.parameters
    assert (
        a.parameters.mass_spectrum["ms2"].spectral_similarity_search
        is not b.parameters.mass_spectrum["ms2"].spectral_similarity_search
    )

    ss = SpectralSimilaritySearchSettings()
    ss.ms2_min_fe_score = 0.77
    ss.max_ms2_tolerance_in_da = 0.05
    coll = [a, b]
    apply_spectral_similarity_search_to_collection(coll, ss, profile="ms2")

    assert a.parameters.mass_spectrum["ms2"].spectral_similarity_search.ms2_min_fe_score == 0.77
    assert b.parameters.mass_spectrum["ms2"].spectral_similarity_search.ms2_min_fe_score == 0.77
    assert a.parameters.mass_spectrum["ms2"].spectral_similarity_search.max_ms2_tolerance_in_da == 0.05
    assert b.parameters.mass_spectrum["ms2"].spectral_similarity_search.max_ms2_tolerance_in_da == 0.05
    # each sample has its own copy
    assert (
        a.parameters.mass_spectrum["ms2"].spectral_similarity_search
        is not b.parameters.mass_spectrum["ms2"].spectral_similarity_search
    )
    a.parameters.mass_spectrum["ms2"].spectral_similarity_search.ms2_min_fe_score = 0.1
    assert b.parameters.mass_spectrum["ms2"].spectral_similarity_search.ms2_min_fe_score == 0.77


def test_resolve_fe_kwargs_settings_base_fe_overrides():
    s = SpectralSimilaritySearchSettings(max_ms2_tolerance_in_da=0.01)
    merged = _resolve_fe_kwargs(settings=s, fe_kwargs={"max_indexed_mz": 1000})
    assert merged["max_ms2_tolerance_in_da"] == 0.01
    assert merged["max_indexed_mz"] == 1000


def test_export_nested_spectral_similarity_search_omits_legacy_lc_ms_annotation():
    reset_lcms_parameters()
    obj = _FakeLCMS()
    obj.parameters.mass_spectrum["ms2"].spectral_similarity_search.ms2_min_fe_score = 0.33
    d = parameter_to_dict.get_dict_data_lcms(obj)
    assert d["mass_spectrum"]["ms2"]["spectral_similarity_search"]["ms2_min_fe_score"] == 0.33
    for key in ("ms2_min_fe_score", "search_as_lipids", "include_fragment_types"):
        assert key not in d["LiquidChromatograph"]


def test_legacy_toml_annotation_imports_to_ms2_profile(tmp_path):
    reset_lcms_parameters()
    toml_path = tmp_path / "legacy.toml"
    toml_path.write_text(
        textwrap.dedent(
            """\
            [LiquidChromatograph]
            eic_tolerance_ppm = 7.5
            ms2_min_fe_score = 0.42
            search_as_lipids = true
            include_fragment_types = true

            [mass_spectrum]
            """
        ),
        encoding="utf-8",
    )
    obj = _FakeLCMS()
    load_and_set_toml_parameters_lcms(obj, parameters_path=str(toml_path))
    assert obj.parameters.mass_spectrum["ms2"].spectral_similarity_search.ms2_min_fe_score == 0.42
    assert obj.parameters.mass_spectrum["ms2"].spectral_similarity_search.search_as_lipids is True
    assert obj.parameters.lc_ms.eic_tolerance_ppm == 7.5


def test_legacy_annotation_fields_load_then_save_omitted_from_toml(tmp_path):
    """Old TOML puts annotation under LiquidChromatograph; re-export must not rewrite them there.

    Values are preserved on spectral_similarity_search (and still loadable from
    the nested section). Legacy fields use ``legacy_field`` so
    ``settings_to_export_dict`` / dump_lcms_settings_toml skip them.
    """
    reset_lcms_parameters()
    legacy_path = tmp_path / "legacy_in.toml"
    legacy_path.write_text(
        textwrap.dedent(
            """\
            [LiquidChromatograph]
            eic_tolerance_ppm = 6.0
            ms2_min_fe_score = 0.42
            search_as_lipids = true
            include_fragment_types = true

            [mass_spectrum]
            """
        ),
        encoding="utf-8",
    )

    obj = _FakeLCMS()
    load_and_set_toml_parameters_lcms(obj, parameters_path=str(legacy_path))

    # Loaded into the nested settings object
    s = obj.parameters.mass_spectrum["ms2"].spectral_similarity_search
    assert s.ms2_min_fe_score == 0.42
    assert s.search_as_lipids is True
    assert s.include_fragment_types is True
    # Still present in-memory on legacy lc_ms aliases
    assert obj.parameters.lc_ms.ms2_min_fe_score == 0.42

    out_path = tmp_path / "reexported.toml"
    dump_lcms_settings_toml(file_path=out_path, lcms_obj=obj)
    saved = out_path.read_text(encoding="utf-8")

    # LiquidChromatograph section must not re-emit legacy annotation keys
    # (toml may order keys freely; parse via reload)
    obj2 = _FakeLCMS()
    load_and_set_toml_parameters_lcms(obj2, parameters_path=str(out_path))
    # Round-trip values via nested path
    s2 = obj2.parameters.mass_spectrum["ms2"].spectral_similarity_search
    assert s2.ms2_min_fe_score == 0.42
    assert s2.search_as_lipids is True
    assert s2.include_fragment_types is True
    assert obj2.parameters.lc_ms.eic_tolerance_ppm == 6.0

    import toml as _toml

    data = _toml.load(out_path)
    lc_section = data.get("LiquidChromatograph") or {}
    for key in _LEGACY_ANNOTATION_KEYS:
        assert key not in lc_section, f"legacy key {key!r} should not be saved under LiquidChromatograph"
    # Canonical location still has the values
    nested = data["mass_spectrum"]["ms2"]["spectral_similarity_search"]
    assert nested["ms2_min_fe_score"] == 0.42
    assert nested["search_as_lipids"] is True
    assert nested["include_fragment_types"] is True


def test_nested_spectral_similarity_search_toml_loads(tmp_path):
    reset_lcms_parameters()
    toml_path = tmp_path / "nested.toml"
    toml_path.write_text(
        textwrap.dedent(
            """\
            [LiquidChromatograph]
            eic_tolerance_ppm = 5.0

            [mass_spectrum.ms2.spectral_similarity_search]
            ms2_min_fe_score = 0.55
            search_type = "identity"
            max_ms2_tolerance_in_da = 0.02
            min_ms2_difference_in_da = 0.99

            [mass_spectrum]
            """
        ),
        encoding="utf-8",
    )
    # Fix toml structure - mass_spectrum.ms2 needs to exist
    toml_path.write_text(
        textwrap.dedent(
            """\
            [LiquidChromatograph]
            eic_tolerance_ppm = 5.0

            [mass_spectrum.ms2.spectral_similarity_search]
            ms2_min_fe_score = 0.55
            search_type = "identity"
            max_ms2_tolerance_in_da = 0.02
            min_ms2_difference_in_da = 0.99
            """
        ),
        encoding="utf-8",
    )
    obj = _FakeLCMS()
    load_and_set_toml_parameters_lcms(obj, parameters_path=str(toml_path))
    s = obj.parameters.mass_spectrum["ms2"].spectral_similarity_search
    assert s.ms2_min_fe_score == 0.55
    assert s.search_type == "identity"
    assert s.max_ms2_tolerance_in_da == 0.02
    assert s.min_ms2_difference_in_da == 0.04  # derived; stale min ignored


def test_legacy_field_metadata():
    for name in ("ms2_min_fe_score", "search_as_lipids", "include_fragment_types"):
        f = next(
            x for x in dataclasses.fields(LiquidChromatographSetting) if x.name == name
        )
        assert f.metadata.get(LEGACY_PARAMETER) is True
    exported = settings_to_export_dict(LiquidChromatographSetting())
    assert "ms2_min_fe_score" not in exported


def test_molecular_network_settings():
    s = SpectralSimilaritySearchSettings(search_type="identity")
    net = MolecularNetwork(fe_lib=None, settings=s)
    assert net.search_type == "identity"
    net2 = MolecularNetwork(fe_lib=None, settings=s, search_type="open")
    assert net2.search_type == "open"


def test_annotation_sync_helpers():
    reset_lcms_parameters()
    p = LCMSParameters(use_defaults=True)
    p.lc_ms.ms2_min_fe_score = 0.61
    p.sync_ms2_annotation_from_lc_ms()
    assert p.mass_spectrum["ms2"].spectral_similarity_search.ms2_min_fe_score == 0.61
    p.mass_spectrum["ms2"].spectral_similarity_search.include_fragment_types = True
    p.sync_ms2_annotation_to_lc_ms()
    assert p.lc_ms.include_fragment_types is True
