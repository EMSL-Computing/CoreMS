"""HDF5 export/import tests for profile mass spectra.

None of these need the molecular formula database (no ``molecular_db`` mark):
they exercise raw_ms layout, peak-index round-tripping, calibration, the
MassSpecCentroid -> MassSpecProfile conversion, compression options, legacy
unlabelled files, ms_unprocessed column labels, and that the HDF5 reader
module imports without a .NET/pythonnet runtime.
"""

import json
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

from corems.mass_spectrum.factory.MassSpectrumClasses import (
    MassSpecCentroid,
    MassSpecProfile,
)
from corems.mass_spectrum.input.coremsHDF5 import ReadCoreMSHDF_MassSpectrum
from corems.mass_spectrum.input.numpyArray import ms_from_array_profile
from corems.mass_spectrum.output.export import HighResMassSpecExport, hdf5_dataset_kwargs
from corems.mass_spectra.factory.lc_class import LCMSBase, MassSpectraBase
from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectra
from corems.mass_spectra.output.export import HighResMassSpectraExport, LCMSExport


def make_profile_ms(dataname="mwe", polarity=-1, scan_number=0, centers=((300.1, 1e6),)):
    """Synthetic profile spectrum with well-separated, well-sampled Gaussian peaks.

    Peak sigma (0.05) spans several grid points (grid spacing ~0.01) so the
    quadratic apex fit and prominence/minima checks see a real peak shape,
    not 1-2 isolated samples. The baseline carries a little deterministic
    noise rather than a dead flat 10.0 everywhere: find_minima() (shared with
    real peak picking) walks outward from an apex while abundance is
    non-increasing, and on a perfectly flat baseline that walk never stops
    until it runs into the next peak or the array edge, overstating peak
    width/area in a way real (noisy) data wouldn't.
    """
    mz = np.linspace(200.0, 800.0, 60_000)
    abundance = 10.0 + np.random.default_rng(0).normal(0, 1e-6, size=mz.shape)
    for centre, height in centers:
        abundance += height * np.exp(-0.5 * ((mz - centre) / 0.05) ** 2)
    ms = ms_from_array_profile(mz, abundance, dataname, polarity=polarity, auto_process=False)
    ms.settings.noise_threshold_method = "absolute_abundance"
    ms.settings.noise_threshold_absolute_abundance = 1000
    ms.process_mass_spec()
    ms.scan_number = scan_number
    return ms


def export_single(ms, out_path, **to_hdf_kwargs):
    export = HighResMassSpecExport(str(out_path), ms)
    export._output_type = "hdf5"
    export.save() if not to_hdf_kwargs else export.to_hdf(**to_hdf_kwargs)
    return out_path.with_suffix(".hdf5")


# ---------------------------------------------------------------------------
# Profile round-trip
# ---------------------------------------------------------------------------


def test_single_spectrum_profile_round_trip(tmp_path):
    ms = make_profile_ms(centers=((300.1, 1e6), (450.2, 5e5), (600.3, 2e5)))
    h5path = export_single(ms, tmp_path / "spec", profile_dtype="f8")

    with ReadCoreMSHDF_MassSpectrum(str(h5path)) as reader:
        ms2 = reader.get_mass_spectrum(load_raw=True, load_molecular_formula=False)

    assert isinstance(ms2, MassSpecProfile)
    np.testing.assert_allclose(ms2.mz_exp_profile, ms._mz_exp)
    np.testing.assert_allclose(ms2.abundance_profile, ms.abundance_profile)
    assert len(ms2.mspeaks) == len(ms.mspeaks)
    assert len(ms2) > 0
    # settings carried over
    assert ms2.settings.noise_threshold_method == "absolute_abundance"
    # instrument / acquisition identity and derived stats survive to_profile()
    assert ms2.polarity == ms.polarity
    assert ms2.scan_number == ms.scan_number
    assert ms2.retention_time == ms.retention_time
    assert ms2.analyzer == ms.analyzer
    assert ms2.instrument_label == ms.instrument_label
    # filename legitimately differs: the reloaded spectrum's filename is the
    # HDF5 file it was read from (set from the reader's file_location, same
    # as any MassSpecCentroid from this reader), not the original source
    # name ("mwe") the export was created from -- that name isn't stored
    # anywhere in the file. to_profile() must carry the centroid's _filename
    # through unchanged rather than losing or resetting it.
    assert str(ms2.filename) == str(h5path)
    np.testing.assert_allclose(ms2.baseline_noise, ms.baseline_noise, rtol=1e-5)
    np.testing.assert_allclose(ms2.baseline_noise_std, ms.baseline_noise_std, rtol=1e-5)
    np.testing.assert_allclose(ms2.dynamic_range, ms.dynamic_range, rtol=1e-5)


class _DummyParser:
    """Stand-in parser class; only its __name__ is recorded on export."""


def test_collection_profile_round_trip(tmp_path):
    ms1 = make_profile_ms(scan_number=1, centers=((300.1, 1e6),))
    ms2 = make_profile_ms(scan_number=2, centers=((450.2, 5e5), (600.3, 2e5)))

    raw_stub = tmp_path / "source.raw"
    raw_stub.touch()
    coll = LCMSBase(str(raw_stub), analyzer="ICR", instrument_label="21T")
    coll.spectra_parser_class = _DummyParser
    coll.polarity = "negative"
    coll.add_mass_spectrum(ms1)
    coll.add_mass_spectrum(ms2)

    out = str(tmp_path / "coll")
    HighResMassSpectraExport(out, coll).to_hdf(export_raw=True, profile_dtype="f8")
    h5path = tmp_path / "coll.corems" / "coll.hdf5"

    # get_mass_spectra_obj() (-> import_mass_spectra()) is the real collection
    # load path and, unlike calling get_mass_spectrum() directly, is what
    # actually sets .scan_number on each loaded spectrum (get_mass_spectrum()
    # does not set it from its own scan_number argument; that is the
    # caller's job, pre-existing and unrelated to to_profile()).
    with ReadCoreMSHDFMassSpectra(str(h5path)) as reader:
        spectra_obj = reader.get_mass_spectra_obj(load_raw=True)
        for scan, original in [(1, ms1), (2, ms2)]:
            reloaded = spectra_obj[scan]
            assert isinstance(reloaded, MassSpecProfile)
            np.testing.assert_allclose(reloaded.mz_exp_profile, original._mz_exp)
            np.testing.assert_allclose(reloaded.abundance_profile, original.abundance_profile)
            assert len(reloaded.mspeaks) == len(original.mspeaks)
            assert reloaded.polarity == original.polarity
            assert reloaded.scan_number == original.scan_number
            assert reloaded.retention_time == original.retention_time
            # analyzer/instrument_label are collection-level HDF5 attrs (set
            # once from coll.analyzer/instrument_label at export time), not
            # per-spectrum: each original ms1/ms2 kept its own creation-time
            # "Generic Simulated" analyzer, which was never the collection's.
            assert reloaded.analyzer == coll.analyzer
            assert reloaded.instrument_label == coll.instrument_label
            np.testing.assert_allclose(reloaded.baseline_noise, original.baseline_noise, rtol=1e-5)
            np.testing.assert_allclose(
                reloaded.baseline_noise_std, original.baseline_noise_std, rtol=1e-5
            )
            np.testing.assert_allclose(reloaded.dynamic_range, original.dynamic_range, rtol=1e-5)


def test_centroid_only_file_unchanged(tmp_path):
    """A file with no raw_ms group must still return a plain MassSpecCentroid."""
    ms = make_profile_ms()
    # to_hdf() has no export_raw switch; drive add_mass_spectrum_to_hdf5
    # directly with export_raw=False to produce a profile-less (centroid-only)
    # file, matching how a centroid spectrum (is_centroid=True) is written.
    export = HighResMassSpecExport(str(tmp_path / "centroid_spec"), ms)
    with h5py.File((tmp_path / "centroid_spec").with_suffix(".hdf5"), "a") as handle:
        export.add_mass_spectrum_to_hdf5(handle, ms, str(ms.scan_number), export_raw=False)
    h5path = (tmp_path / "centroid_spec").with_suffix(".hdf5")

    with ReadCoreMSHDF_MassSpectrum(str(h5path)) as reader:
        ms2 = reader.get_mass_spectrum(load_raw=True, load_molecular_formula=False)

    assert type(ms2) is MassSpecCentroid
    assert len(ms2) == len(ms)


def test_force_centroid_returns_centroid_even_with_profile(tmp_path):
    ms = make_profile_ms()
    h5path = export_single(ms, tmp_path / "spec")

    with ReadCoreMSHDF_MassSpectrum(str(h5path)) as reader:
        ms2 = reader.get_mass_spectrum(
            load_raw=True, load_molecular_formula=False, force_centroid=True
        )

    assert type(ms2) is MassSpecCentroid


# ---------------------------------------------------------------------------
# Peak indices
# ---------------------------------------------------------------------------


def test_peak_indices_round_trip_exact(tmp_path):
    """With stored peak_profile_indices, area and apex m/z match exactly (f8)."""
    ms = make_profile_ms(centers=((300.1, 1e6), (450.2, 5e5)))
    h5path = export_single(ms, tmp_path / "spec", profile_dtype="f8")

    with ReadCoreMSHDF_MassSpectrum(str(h5path)) as reader:
        ms2 = reader.get_mass_spectrum(load_raw=True, load_molecular_formula=False)

    original_by_mz = sorted(ms.mspeaks, key=lambda p: p.mz_exp)
    reloaded_by_mz = sorted(ms2.mspeaks, key=lambda p: p.mz_exp)
    assert len(original_by_mz) == len(reloaded_by_mz)
    for orig, reloaded in zip(original_by_mz, reloaded_by_mz):
        assert reloaded.peak_left_index == orig.peak_left_index
        assert reloaded.peak_apex_index == orig.peak_apex_index
        assert reloaded.peak_right_index == orig.peak_right_index
        np.testing.assert_allclose(reloaded.area, orig.area, rtol=1e-10)
        np.testing.assert_allclose(reloaded.mz_exp, orig.mz_exp, rtol=1e-10)


def test_peak_indices_fallback_without_stored_indices(tmp_path):
    """Without peak_profile_indices (older file), nearest-mz + find_minima
    recovers bounds close enough that area/apex still match within tolerance.
    """
    ms = make_profile_ms(centers=((300.1, 1e6), (450.2, 5e5)))
    h5path = export_single(ms, tmp_path / "spec", profile_dtype="f8")

    with h5py.File(h5path, "a") as handle:
        del handle[str(ms.scan_number)]["peak_profile_indices"]

    with ReadCoreMSHDF_MassSpectrum(str(h5path)) as reader:
        ms2 = reader.get_mass_spectrum(load_raw=True, load_molecular_formula=False)

    original_by_mz = sorted(ms.mspeaks, key=lambda p: p.mz_exp)
    reloaded_by_mz = sorted(ms2.mspeaks, key=lambda p: p.mz_exp)
    assert len(original_by_mz) == len(reloaded_by_mz)
    for orig, reloaded in zip(original_by_mz, reloaded_by_mz):
        np.testing.assert_allclose(reloaded.mz_exp, orig.mz_exp, rtol=1e-6)
        np.testing.assert_allclose(reloaded.area, orig.area, rtol=0.05)


# ---------------------------------------------------------------------------
# Calibration
# ---------------------------------------------------------------------------


def test_calibrated_profile_round_trip_no_double_calibration(tmp_path):
    ms = make_profile_ms()
    # Simulate a profile-level calibration: a small constant ppm shift.
    shift_ppm = 2.0
    ms.mz_cal_profile = ms._mz_exp * (1 + shift_ppm * 1e-6)

    h5path = export_single(ms, tmp_path / "spec", profile_dtype="f8")

    with h5py.File(h5path, "r") as handle:
        row_labels = json.loads(handle[str(ms.scan_number)]["raw_ms"].attrs["row_labels"])
    assert row_labels == ["mz", "abundance", "mz_cal"]

    with ReadCoreMSHDF_MassSpectrum(str(h5path)) as reader:
        ms2 = reader.get_mass_spectrum(load_raw=True, load_molecular_formula=False)

    assert ms2.is_calibrated is True
    np.testing.assert_allclose(ms2.mz_exp_profile, ms.mz_cal_profile)
    # Not double-applied: the calibrated profile is exactly what was stored,
    # not calibration-of-the-calibration.
    np.testing.assert_allclose(np.asarray(ms2.mz_exp_profile), ms._mz_exp * (1 + shift_ppm * 1e-6))


def test_uncalibrated_profile_has_no_mz_cal_row(tmp_path):
    ms = make_profile_ms()
    h5path = export_single(ms, tmp_path / "spec")

    with h5py.File(h5path, "r") as handle:
        row_labels = json.loads(handle[str(ms.scan_number)]["raw_ms"].attrs["row_labels"])
    assert row_labels == ["mz", "abundance"]

    with ReadCoreMSHDF_MassSpectrum(str(h5path)) as reader:
        ms2 = reader.get_mass_spectrum(load_raw=True, load_molecular_formula=False)
    assert ms2.is_calibrated is False


# ---------------------------------------------------------------------------
# Legacy unlabelled files
# ---------------------------------------------------------------------------


def test_legacy_unlabelled_file_reads_correct_values(tmp_path):
    """A file written in the pre-fix [abundance, mz] order with no
    row_labels attribute must still read back CORRECT (unswapped) m/z and
    abundance -- not reproduce the historical swap.
    """
    ms = make_profile_ms()
    h5path = export_single(ms, tmp_path / "spec", profile_dtype="f8")

    with h5py.File(h5path, "a") as handle:
        raw_ms = handle[str(ms.scan_number)]["raw_ms"]
        data = raw_ms[()]
        # Rewrite in the legacy [abundance, mz] order and drop row_labels,
        # simulating a file written before this fix existed.
        legacy_order = np.stack([data[1], data[0]])  # [abundance, mz]
        del handle[str(ms.scan_number)]["raw_ms"]
        legacy_dataset = handle[str(ms.scan_number)].create_dataset("raw_ms", data=legacy_order)
        legacy_dataset.attrs["MassSpecAttrs"] = raw_ms.attrs.get("MassSpecAttrs", "{}")

    with ReadCoreMSHDF_MassSpectrum(str(h5path)) as reader:
        ms2 = reader.get_mass_spectrum(load_raw=True, load_molecular_formula=False)

    np.testing.assert_allclose(ms2.mz_exp_profile, ms._mz_exp)
    np.testing.assert_allclose(ms2.abundance_profile, ms.abundance_profile)


# ---------------------------------------------------------------------------
# Compression / dtype
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("compression", [0, 4, 9, "lzf", None])
def test_profile_dtype_and_compression_options_round_trip(tmp_path, compression):
    ms = make_profile_ms()
    h5path = export_single(
        ms, tmp_path / f"spec_{compression}", compression=compression, profile_dtype="f4"
    )

    with ReadCoreMSHDF_MassSpectrum(str(h5path)) as reader:
        ms2 = reader.get_mass_spectrum(load_raw=True, load_molecular_formula=False)

    assert isinstance(ms2, MassSpecProfile)
    np.testing.assert_allclose(
        ms2.mz_exp_profile, ms._mz_exp.astype(np.float32), rtol=1e-6
    )


@pytest.mark.parametrize("profile_dtype", ["f4", "f8"])
def test_profile_dtype_option(tmp_path, profile_dtype):
    ms = make_profile_ms()
    h5path = export_single(ms, tmp_path / f"spec_{profile_dtype}", profile_dtype=profile_dtype)

    with h5py.File(h5path, "r") as handle:
        raw_ms = handle[str(ms.scan_number)]["raw_ms"]
        assert raw_ms.dtype == np.dtype(profile_dtype)


@pytest.mark.parametrize("compression", [0, 4, 9, "lzf", None])
def test_hdf5_dataset_kwargs_handles_vlen_strings(tmp_path, compression):
    """h5py refuses fletcher32 on variable-length string datasets; the
    helper must omit it (and shuffle) for every compression option.
    """
    data = ["abc", "de", "fghij"]
    dtype = h5py.string_dtype(encoding="utf-8")
    kwargs = hdf5_dataset_kwargs((len(data),), compression=compression, vlen_str=True)
    assert "shuffle" not in kwargs
    assert "fletcher32" not in kwargs

    with h5py.File(tmp_path / "vlen.hdf5", "w") as handle:
        # Must not raise.
        handle.create_dataset("strings", data=data, dtype=dtype, **kwargs)


def test_hdf5_dataset_kwargs_scalar_and_empty():
    assert hdf5_dataset_kwargs(()) == {}
    assert hdf5_dataset_kwargs((0,)) == {}
    assert hdf5_dataset_kwargs((0, 3)) == {}


# ---------------------------------------------------------------------------
# ms_unprocessed column labels
# ---------------------------------------------------------------------------


def test_ms_unprocessed_column_labels_round_trip(tmp_path):
    import pandas as pd

    raw_stub = tmp_path / "source.raw"
    raw_stub.touch()
    coll = LCMSBase(str(raw_stub), analyzer="ICR", instrument_label="21T")
    coll.spectra_parser_class = _DummyParser
    coll.polarity = "negative"
    ms = make_profile_ms(scan_number=1)
    coll.add_mass_spectrum(ms)
    coll._ms_unprocessed[1] = pd.DataFrame(
        {"scan": [1, 1], "mz": [300.0, 301.0], "intensity": [10.0, 20.0]}
    )

    out = str(tmp_path / "coll_unproc")
    coll.parameters.lc_ms.export_unprocessed_ms1 = True
    LCMSExport(out, coll).to_hdf()
    h5path = tmp_path / "coll_unproc.corems" / "coll_unproc.hdf5"

    with h5py.File(h5path, "r") as handle:
        dataset = handle["ms_unprocessed"]["1"]
        assert json.loads(dataset.attrs["column_labels"]) == ["scan", "mz", "intensity"]

    with ReadCoreMSHDFMassSpectra(str(h5path)) as reader:
        ms_raw = reader.get_ms_raw()
    assert list(ms_raw[1].columns) == ["scan", "mz", "intensity"]
    np.testing.assert_array_equal(ms_raw[1]["mz"].to_numpy(), [300.0, 301.0])


# ---------------------------------------------------------------------------
# MassSpectraBase.__iter__
# ---------------------------------------------------------------------------


def test_mass_spectra_base_is_iterable(tmp_path):
    raw_stub = tmp_path / "source.raw"
    raw_stub.touch()
    spectra = MassSpectraBase(str(raw_stub), analyzer="ICR", instrument_label="21T")
    ms1 = make_profile_ms(scan_number=1)
    ms2 = make_profile_ms(scan_number=2)
    spectra.add_mass_spectrum(ms1)
    spectra.add_mass_spectrum(ms2)

    collected = list(spectra)
    assert len(collected) == 2
    assert {m.scan_number for m in collected} == {1, 2}


# ---------------------------------------------------------------------------
# Lazy Thermo import
# ---------------------------------------------------------------------------


def test_corems_hdf5_imports_without_clr():
    """Importing corems_hdf5 must not require pythonnet/.NET: rawFileReader
    (which does `import clr`) must not be imported at module scope.

    Runs in a fresh subprocess rather than monkeypatching sys.modules in this
    process: this process may already have pythonnet imported (by an earlier
    test, or just by being on a machine where corems_hdf5's old module-level
    import already ran once at collection time), and under pytest-xdist this
    test's worker process is reused by later tests, so mutating sys.modules
    here could leak into them. A subprocess gets a clean sys.modules with no
    such risk either way.
    """
    script = (
        "import sys; sys.modules['clr'] = None\n"
        "import corems.mass_spectra.input.corems_hdf5 as m\n"
        "assert hasattr(m, 'ReadCoreMSHDFMassSpectra')\n"
        "assert 'corems.mass_spectra.input.rawFileReader' not in sys.modules\n"
        "print('OK')\n"
    )
    result = subprocess.run(
        [sys.executable, "-I", "-c", script],
        capture_output=True,
        text=True,
        cwd=str(Path(__file__).resolve().parent.parent),
    )
    assert result.returncode == 0, (
        f"import failed under a clr-unavailable subprocess:\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )
    assert "OK" in result.stdout
