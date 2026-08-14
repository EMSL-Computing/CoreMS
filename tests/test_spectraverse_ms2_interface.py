"""Tests for SpectraverseMS2Interface and spectral library parquet cache."""

from pathlib import Path

import numpy as np
import pytest
from ms_entropy import FlashEntropySearch

from corems.mass_spectra.output.export import LCMSMetabolomicsExport
from corems.molecular_formula.calc.MolecularFormulaCalc import MolecularFormulaCalc
from corems.molecular_id.search.database_interfaces import (
    MSPInterface,
    SpectraverseMS2Interface,
)
from corems.molecular_id.search import spectral_library_cache


@pytest.fixture
def spectraverse_mgf_path():
    return Path.cwd() / "tests/tests_data/lcms/test_spectraverse.mgf"


@pytest.fixture
def msp_file_location():
    return Path.cwd() / "tests/tests_data/lcms/test_db.msp"


def test_precursor_mz_from_formula_central_and_export_wrapper():
    mz = MolecularFormulaCalc.precursor_mz_from_formula("C24H48O2", "[M-H]-")
    assert abs(mz - 367.3582) < 0.01
    assert LCMSMetabolomicsExport.precursor_mz_from_formula(
        "C24H48O2", "[M-H]-"
    ) == pytest.approx(mz)
    assert LCMSMetabolomicsExport.get_ion_formula(
        "C24H48O2", "[M-H]-"
    ) == MolecularFormulaCalc.get_ion_formula("C24H48O2", "[M-H]-")
    # Spectraverse-style alias
    assert MolecularFormulaCalc.get_ion_formula(
        "C2H4O2", "[M+HCOOH-H]-"
    ) == MolecularFormulaCalc.get_ion_formula("C2H4O2", "[M+HCOO]-")


def test_spectraverse_loads_ms2_only(spectraverse_mgf_path, tmp_path):
    cache_path = tmp_path / "sv.corems-lib.parquet"
    iface = SpectraverseMS2Interface(
        spectraverse_mgf_path, cache=True, cache_path=cache_path
    )
    assert len(iface._data_frame) == 3
    assert set(iface._data_frame["polarity"]) == {"positive", "negative"}
    assert "MS1" not in set(iface._data_frame["ms_level"].astype(str))


def test_spectraverse_precursor_is_calculated(spectraverse_mgf_path, tmp_path):
    cache_path = tmp_path / "sv.corems-lib.parquet"
    iface = SpectraverseMS2Interface(
        spectraverse_mgf_path, cache=False, cache_path=cache_path
    )
    neg = iface._data_frame[iface._data_frame["polarity"] == "negative"].iloc[0]
    # Library experimental PRECURSOR_MZ is 367.359; calculated ion m/z differs slightly
    assert neg["precursor_mz_library"] == pytest.approx(367.359, rel=0, abs=1e-6)
    assert neg["precursor_mz"] != pytest.approx(367.359, rel=0, abs=1e-6)
    assert abs(neg["precursor_mz"] - 367.3582) < 0.01


def test_get_metabolomics_spectra_library_df_and_fe(spectraverse_mgf_path, tmp_path):
    cache_path = tmp_path / "sv.corems-lib.parquet"
    iface = SpectraverseMS2Interface(
        spectraverse_mgf_path, cache=True, cache_path=cache_path
    )
    df_lib, meta = iface.get_metabolomics_spectra_library(
        polarity="negative", format="df", normalize=True
    )
    assert len(df_lib) == 2
    assert "peaks" in df_lib.columns
    assert "precursor_mz" in df_lib.columns
    assert "ion_type" in df_lib.columns
    assert len(meta) == 1
    inchikey = "QZZGJDVWLFXDLK-UHFFFAOYSA-N"
    assert inchikey in meta
    assert meta[inchikey].formula == "C24H48O2"

    fe_lib, meta_fe = iface.get_metabolomics_spectra_library(
        polarity="positive",
        format="flashentropy",
        normalize=True,
        fe_kwargs={
            "min_ms2_difference_in_da": 0.02,
            "max_ms2_tolerance_in_da": 0.01,
        },
    )
    assert isinstance(fe_lib, FlashEntropySearch)
    assert len(meta_fe) == 1


def test_spectraverse_library_settings_not_constructor(spectraverse_mgf_path, tmp_path):
    """FE knobs belong on get_metabolomics_spectra_library, not __init__."""
    from corems.encapsulation.factory.processingSetting import (
        SpectralSimilaritySearchSettings,
    )

    ss = SpectralSimilaritySearchSettings()
    ss.max_ms2_tolerance_in_da = 0.01

    with pytest.raises(TypeError):
        SpectraverseMS2Interface(
            spectraverse_mgf_path,
            settings=ss,
            cache=False,
            cache_path=tmp_path / "sv_bad.parquet",
        )

    iface = SpectraverseMS2Interface(
        spectraverse_mgf_path, cache=False, cache_path=tmp_path / "sv_ok.parquet"
    )
    fe_lib, meta = iface.get_metabolomics_spectra_library(
        polarity="negative",
        format="flashentropy",
        normalize=True,
        settings=ss,
    )
    assert isinstance(fe_lib, FlashEntropySearch)
    assert len(meta) == 1


def test_spectraverse_msp_flashentropy_search_parity(
    spectraverse_mgf_path, msp_file_location, tmp_path
):
    """Shared Lignoceric MS2 spectrum: MSP and Spectraverse FE libs return same top hit and score."""
    fe_kwargs = {
        "normalize_intensity": True,
        "min_ms2_difference_in_da": 0.02,
        "max_ms2_tolerance_in_da": 0.01,
        "max_indexed_mz": 3000,
        "precursor_ions_removal_da": None,
        "noise_threshold": 0,
    }
    msp = MSPInterface(msp_file_location, cache=False)
    sv = SpectraverseMS2Interface(
        spectraverse_mgf_path,
        cache=False,
        cache_path=tmp_path / "sv_parity.parquet",
    )
    msp_lib, _ = msp.get_metabolomics_spectra_library(
        polarity="negative",
        format="flashentropy",
        normalize=True,
        fe_kwargs=fe_kwargs,
    )
    sv_lib, _ = sv.get_metabolomics_spectra_library(
        polarity="negative",
        format="flashentropy",
        normalize=True,
        fe_kwargs=fe_kwargs,
    )

    # Identical peaks in test_db.msp (Lignoceric Acid) and test_spectraverse.mgf (NEG_1)
    peaks = np.array(
        [[367.359192, 999.0], [368.362396, 245.0]], dtype=np.float64
    )
    precursor_mz = 367.359
    expected_inchikey = "QZZGJDVWLFXDLK-UHFFFAOYSA-N"

    def top_identity_hit(fe_lib):
        query = fe_lib.clean_spectrum_for_search(
            precursor_mz=precursor_mz,
            peaks=peaks,
            precursor_ions_removal_da=None,
            noise_threshold=0,
            min_ms2_difference_in_da=0.02,
        )
        scores = fe_lib.search(
            precursor_mz=precursor_mz,
            peaks=query,
            ms1_tolerance_in_da=0.05,
            ms2_tolerance_in_da=0.005,
            method={"identity"},
            precursor_ions_removal_da=None,
            noise_threshold=0,
            target="cpu",
        )["identity_search"]
        idx = int(np.argmax(scores))
        return fe_lib[idx]["molecular_data_id"], float(scores[idx])

    msp_id, msp_score = top_identity_hit(msp_lib)
    sv_id, sv_score = top_identity_hit(sv_lib)

    assert msp_id == expected_inchikey
    assert sv_id == expected_inchikey
    assert msp_id == sv_id
    assert msp_score == pytest.approx(sv_score, rel=0, abs=1e-6)
    assert msp_score > 0.99


def test_spectraverse_bad_path():
    with pytest.raises(FileNotFoundError):
        SpectraverseMS2Interface("/no/such/file.mgf", cache=False)


def test_spectraverse_bad_polarity(spectraverse_mgf_path, tmp_path):
    iface = SpectraverseMS2Interface(
        spectraverse_mgf_path, cache=False, cache_path=tmp_path / "x.parquet"
    )
    with pytest.raises(ValueError, match="Polarity"):
        iface.get_metabolomics_spectra_library(polarity="both", format="df")


def test_spectraverse_cache_roundtrip(spectraverse_mgf_path, tmp_path):
    cache_path = tmp_path / "sv.corems-lib.parquet"
    iface1 = SpectraverseMS2Interface(
        spectraverse_mgf_path, cache=True, cache_path=cache_path, rebuild_cache=True
    )
    assert cache_path.is_file()
    n1 = len(iface1._data_frame)
    p1 = iface1._data_frame["precursor_mz"].tolist()

    iface2 = SpectraverseMS2Interface(
        spectraverse_mgf_path, cache="read", cache_path=cache_path
    )
    assert len(iface2._data_frame) == n1
    assert iface2._data_frame["precursor_mz"].tolist() == pytest.approx(p1)
    peaks0 = iface2._data_frame.iloc[0]["peaks"]
    assert isinstance(peaks0, np.ndarray)
    assert peaks0.ndim == 2 and peaks0.shape[1] == 2


def test_msp_cache_roundtrip(msp_file_location, tmp_path):
    cache_path = tmp_path / "msp.corems-lib.parquet"
    msp1 = MSPInterface(
        msp_file_location, cache=True, cache_path=cache_path, rebuild_cache=True
    )
    assert cache_path.is_file()
    n1 = len(msp1._data_frame)

    msp2 = MSPInterface(msp_file_location, cache="read", cache_path=cache_path)
    assert len(msp2._data_frame) == n1
    assert "peaks" in msp2._data_frame.columns


def test_cache_invalid_when_source_changes(spectraverse_mgf_path, tmp_path):
    src = tmp_path / "lib.mgf"
    src.write_text(spectraverse_mgf_path.read_text())
    cache_path = tmp_path / "lib.mgf.corems-lib.parquet"
    SpectraverseMS2Interface(src, cache=True, cache_path=cache_path)
    assert spectral_library_cache.is_valid(
        cache_path,
        src,
        extra_meta={
            "source_kind": "spectraverse_mgf",
            "precursor_policy": "formula+adduct",
        },
    )
    src.write_text(src.read_text() + "\n")
    assert not spectral_library_cache.is_valid(
        cache_path,
        src,
        extra_meta={
            "source_kind": "spectraverse_mgf",
            "precursor_policy": "formula+adduct",
        },
    )
