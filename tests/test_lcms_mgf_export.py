from pathlib import Path

import numpy as np
import pytest

import pandas as pd

from corems.chroma_peak.factory.chroma_peak_classes import LCMSMassFeature
from corems.mass_spectra.factory.lc_class import LCMSBase, LCMSCollection
from corems.mass_spectrum.input.numpyArray import ms_from_array_centroid
from corems.mass_spectra.output.mgf import (
    feature_mgf_text,
    format_ion_block,
    iter_collection_mgf_records,
    sirius_charge,
    write_feature_records_to_mgf,
)
from corems.molecular_id.factory.spectrum_search_results import SpectrumSearchResults


def test_sirius_charge():
    assert sirius_charge("positive") == "1+"
    assert sirius_charge("negative") == "1-"
    with pytest.raises(ValueError, match="positive"):
        sirius_charge("unknown")


def test_format_ion_block_sirius_keys():
    text = format_ion_block(
        {
            "FEATURE_ID": 21,
            "PEPMASS": 438.32382,
            "CHARGE": "1+",
            "MSLEVEL": 2,
        },
        mz=np.array([185.041199]),
        intensity=np.array([4034.674316]),
    )
    assert text.startswith("BEGIN IONS\n")
    assert text.strip().endswith("END IONS")
    assert "FEATURE_ID=21" in text
    assert "PEPMASS=438.32382" in text
    assert "CHARGE=1+" in text
    assert "MSLEVEL=2" in text
    assert "185.041199 4034.674316" in text


class _DummyLCMS:
    def __init__(self, polarity="positive", sample_name="sample_a"):
        self.mass_features = {}
        self.polarity = polarity
        self.sample_name = sample_name

    def get_time_of_scan_id(self, scan):
        return 5.2


def _centroid_ms(mz, abundance, name="spec"):
    rp = [10000.0] * len(mz)
    s2n = [10.0] * len(mz)
    spec = ms_from_array_centroid(
        mz, abundance, rp, s2n, name, polarity=1, auto_process=False
    )
    spec.settings.noise_threshold_method = "relative_abundance"
    spec.settings.noise_threshold_min_relative_abundance = 0
    spec.process_mass_spec()
    return spec


def _feature_with_spectra(pepmass=438.32, with_ms1=True, with_ms2=True, extra_ms2=False):
    parent = _DummyLCMS()
    feat = LCMSMassFeature(
        parent,
        mz=pepmass,
        retention_time=5.2,
        intensity=1000.0,
        apex_scan=154,
        id=21,
    )
    parent.mass_features[21] = feat
    if with_ms1:
        feat.mass_spectrum = _centroid_ms(
            [100.0, pepmass, pepmass + 1.003, 800.0],
            [5.0, 1000.0, 210.0, 8.0],
            name="ms1",
        )
    if with_ms2:
        feat.ms2_scan_numbers = [160]
        feat.ms2_mass_spectra[160] = _centroid_ms(
            [185.041199, 203.052597],
            [4034.674316, 12382.624023],
            name="ms2a",
        )
        feat.ms2_mass_spectra[160].scan_number = 160
        if extra_ms2:
            feat.ms2_scan_numbers.append(161)
            feat.ms2_mass_spectra[161] = _centroid_ms(
                [100.1],
                [50.0],
                name="ms2b",
            )
            feat.ms2_mass_spectra[161].scan_number = 161
    return feat


def test_feature_mgf_text_pairs_precursor_ms1_and_ms2():
    feat = _feature_with_spectra()
    text = feature_mgf_text(
        feat, feature_id=21, polarity="positive", sample_name="sample_a"
    )
    assert text is not None
    assert text.count("BEGIN IONS") == 2
    assert text.count("MSLEVEL=1") == 1
    assert text.count("MSLEVEL=2") == 1
    assert "FEATURE_ID=21" in text
    assert "PEPMASS=438.32" in text
    assert "CHARGE=1+" in text
    ms1_block = text.split("MSLEVEL=2")[0]
    assert "438.32 1000.0" in ms1_block
    assert "439.323" not in ms1_block
    assert "100.0 5.0" not in text
    assert "185.041199 " in text
    rt_lines = [line for line in text.splitlines() if line.startswith("RTINSECONDS=")]
    assert rt_lines
    assert abs(float(rt_lines[0].split("=", 1)[1]) - 312.0) < 1e-6
    assert "TITLE=sample_a feature 21" in text
    assert "IONMODE" not in text
    assert "SCANS=154" in text
    assert "SCANS=160" in text


def test_feature_mgf_text_ms1_does_not_need_spectrum():
    feat = _feature_with_spectra(with_ms1=False)
    text = feature_mgf_text(feat, 21, "positive")
    assert text is not None
    assert "438.32 1000.0" in text.split("MSLEVEL=2")[0]


def test_feature_mgf_text_skips_incomplete():
    assert feature_mgf_text(_feature_with_spectra(with_ms2=False), 21, "positive") is None


def test_feature_mgf_text_requires_mz():
    feat = _feature_with_spectra()
    feat._mz_exp = None
    feat._mz_cal = None
    with pytest.raises(ValueError, match="m/z"):
        feature_mgf_text(feat, 21, "positive")


def test_ms2_mode_all_writes_multiple_ms2_blocks():
    feat = _feature_with_spectra(extra_ms2=True)
    text = feature_mgf_text(feat, 21, "positive", ms2_mode="all")
    assert text.count("MSLEVEL=1") == 1
    assert text.count("MSLEVEL=2") == 2


def _similarity_result(spec, score):
    return SpectrumSearchResults(
        spec,
        438.32,
        {
            "entropy_similarity": np.array([score]),
            "ref_mol_id": ["mol"],
            "ref_ms_id": ["ref"],
            "ref_precursor_mz": [438.32],
            "precursor_mz_error_ppm": [0.0],
            "ref_ion_type": ["[M+H]+"],
        },
    )


def test_best_mode_writes_highest_similarity_scan():
    """Default ``best`` is ``LCMSMassFeature.best_ms2`` after a search.

    Two MS2 scans are attached. The scan nearer the apex would be chosen
    with no search results. A lower-scoring similarity result on that scan
    and a higher-scoring result on the farther scan select the library hit.
    """
    parent = _DummyLCMS()
    parent.get_time_of_scan_id = lambda scan: {160: 5.21, 200: 7.0}[scan]
    feat = LCMSMassFeature(
        parent,
        mz=438.32,
        retention_time=5.2,
        intensity=1000.0,
        apex_scan=154,
        id=21,
    )
    feat.ms2_scan_numbers = [160, 200]
    feat.ms2_mass_spectra[160] = _centroid_ms(
        [185.041199], [4034.674316], name="apex_ms2"
    )
    feat.ms2_mass_spectra[160].scan_number = 160
    feat.ms2_mass_spectra[200] = _centroid_ms([100.1], [50.0], name="hit_ms2")
    feat.ms2_mass_spectra[200].scan_number = 200

    before_search = feature_mgf_text(feat, 21, "positive")
    assert before_search.count("MSLEVEL=2") == 1
    assert "SCANS=160" in before_search
    assert "SCANS=200" not in before_search
    assert "185.041199 " in before_search
    assert "100.1 " not in before_search

    feat.ms2_similarity_results = [
        _similarity_result(feat.ms2_mass_spectra[160], 0.2),
        _similarity_result(feat.ms2_mass_spectra[200], 0.91),
    ]
    assert feat.best_ms2 is feat.ms2_mass_spectra[200]
    text = feature_mgf_text(feat, 21, "positive")
    assert text.count("BEGIN IONS") == 2
    assert text.count("MSLEVEL=2") == 1
    assert "SCANS=200" in text
    assert "SCANS=160" not in text
    assert "100.1 50.0" in text
    assert "185.041199" not in text


def test_write_feature_records_skips_and_raises_when_empty(tmp_path):
    complete = _feature_with_spectra()
    incomplete = _feature_with_spectra(with_ms2=False)
    incomplete.id = 22
    out = tmp_path / "out.mgf"
    with pytest.warns(UserWarning, match="Skipped 1 of 2"):
        path = write_feature_records_to_mgf(
            [(21, complete, "sample_a"), (22, incomplete, "sample_a")],
            polarity="positive",
            out_file_path=out,
        )
    assert path == out
    body = path.read_text()
    assert "FEATURE_ID=21" in body
    assert "FEATURE_ID=22" not in body

    with pytest.raises(ValueError, match="No complete"):
        write_feature_records_to_mgf(
            [(22, incomplete, "sample_a")],
            polarity="positive",
            out_file_path=tmp_path / "empty.mgf",
        )


def test_write_mgf_overwrite_false(tmp_path):
    feat = _feature_with_spectra()
    out = tmp_path / "dup.mgf"
    write_feature_records_to_mgf([(21, feat, "sample_a")], "positive", out)
    with pytest.raises(FileExistsError):
        write_feature_records_to_mgf([(21, feat, "sample_a")], "positive", out)


def _lcms_with_feature():
    obj = LCMSBase(Path(__file__), sample_name="sample_a")
    obj.polarity = "positive"
    feat = LCMSMassFeature(
        obj,
        mz=438.32,
        retention_time=5.2,
        intensity=1000.0,
        apex_scan=154,
        id=21,
    )
    feat.mass_spectrum = _centroid_ms(
        [438.32, 439.323], [1000.0, 210.0], name="ms1"
    )
    feat.ms2_scan_numbers = [160]
    feat.ms2_mass_spectra[160] = _centroid_ms([185.0], [100.0], name="ms2")
    feat.ms2_mass_spectra[160].scan_number = 160
    obj.mass_features[21] = feat
    return obj, feat


def test_lcmsbase_to_mgf_writes_selected_features(tmp_path):
    obj, feat = _lcms_with_feature()
    extra = LCMSMassFeature(
        obj,
        mz=500.0,
        retention_time=6.0,
        intensity=800.0,
        apex_scan=200,
        id=99,
    )
    extra.mass_spectrum = _centroid_ms([500.0], [800.0], name="ms1b")
    extra.ms2_scan_numbers = [201]
    extra.ms2_mass_spectra[201] = _centroid_ms([120.0], [40.0], name="ms2b")
    obj.mass_features[99] = extra
    path = obj.to_mgf(tmp_path / "export", feature_ids=[21])
    assert path.suffix == ".mgf"
    text = path.read_text()
    assert "FEATURE_ID=21" in text
    assert "FEATURE_ID=99" not in text


def test_lcmsbase_to_mgf_unknown_id_raises(tmp_path):
    obj, _ = _lcms_with_feature()
    with pytest.raises(ValueError, match="not found"):
        obj.to_mgf(tmp_path / "missing.mgf", feature_ids=[12345])


class _FakeCollection:
    def __init__(self, lcms_obj, cluster=7, mf_id=21):
        self._lcms = {lcms_obj.sample_name: lcms_obj}
        self._cluster = cluster
        self._mf_id = mf_id

    def get_representative_mass_features_for_all_clusters(self):
        return pd.DataFrame(
            [
                {
                    "cluster": self._cluster,
                    "sample_name": list(self._lcms.keys())[0],
                    "mf_id": self._mf_id,
                }
            ]
        )


def test_collection_records_use_cluster_as_feature_id():
    obj, feat = _lcms_with_feature()
    records = iter_collection_mgf_records(_FakeCollection(obj, cluster=7, mf_id=21))
    assert records[0][0] == 7
    assert records[0][1] is feat


def test_collection_records_resolve_sample_id():
    obj, feat = _lcms_with_feature()

    class _IdCollection:
        def __init__(self):
            self._lcms = {obj.sample_name: obj}
            self.samples = [obj.sample_name]

        def get_representative_mass_features_for_all_clusters(self):
            return pd.DataFrame(
                [{"cluster": 7, "sample_id": 0, "mf_id": feat.id}]
            )

    records = iter_collection_mgf_records(_IdCollection())
    assert records[0][0] == 7
    assert records[0][1] is feat
    assert records[0][2] == obj.sample_name


def test_collection_unknown_cluster_raises():
    obj, _ = _lcms_with_feature()
    with pytest.raises(ValueError, match="not found"):
        iter_collection_mgf_records(_FakeCollection(obj), cluster_ids=[999])


def test_lcmscollection_to_mgf_feature_id_is_cluster(tmp_path):
    obj, feat = _lcms_with_feature()
    coll = LCMSCollection.__new__(LCMSCollection)
    coll._lcms = {obj.sample_name: obj}
    coll._manifest_dict = {}
    coll.collection_location = "dummy"
    coll.collection_parser = None

    def _reps(representative_metric=None):
        return pd.DataFrame(
            [{"cluster": 7, "sample_name": obj.sample_name, "mf_id": feat.id}]
        )

    coll.get_representative_mass_features_for_all_clusters = _reps
    path = LCMSCollection.to_mgf(coll, tmp_path / "cons.mgf")
    text = path.read_text()
    assert "FEATURE_ID=7" in text
    assert "FEATURE_ID=21" not in text


def test_collection_unloaded_representative_raises():
    obj, feat = _lcms_with_feature()
    obj.mass_features = {}
    with pytest.raises(ValueError, match="load_representatives"):
        iter_collection_mgf_records(_FakeCollection(obj, cluster=7, mf_id=feat.id))
