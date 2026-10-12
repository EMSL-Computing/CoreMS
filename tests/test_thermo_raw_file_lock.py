"""Readers opened for a Thermo .raw file must not outlive the object that owns them.

These tests use stand-ins. They do not open a .raw file and do not need Windows.
"""

import gc
import pickle
import warnings
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from corems.mass_spectra.calc.lc_calc import LCCalculations
from corems.mass_spectra.factory.lc_class import LCMSCollection, MassSpectraBase
from corems.mass_spectra.output.export import HighResMassSpectraExport


class CountingParser:
    """Stand-in parser. ``constructed`` and ``closed`` survive the instance."""

    constructed = 0
    closed = 0

    def __init__(self, path):
        CountingParser.constructed += 1
        self._file_location = Path(path)
        self._sample_name = Path(path).stem
        self._analyzer = "Unknown"
        self._instrument_label = "Unknown"
        self.reads = []

    @property
    def file_location(self):
        self.reads.append("file_location")
        return self._file_location

    @property
    def sample_name(self):
        self.reads.append("sample_name")
        return self._sample_name

    @property
    def analyzer(self):
        self.reads.append("analyzer")
        return self._analyzer

    @property
    def instrument_label(self):
        self.reads.append("instrument_label")
        return self._instrument_label

    def close_file(self):
        CountingParser.closed += 1


class _Disposable:
    def __init__(self, events):
        self.events = events

    def Dispose(self):
        self.events.append("dispose")


class _NoDispose:
    pass


def _reset_counts():
    CountingParser.constructed = 0
    CountingParser.closed = 0


def test_constructor_reuses_passed_parser(tmp_path):
    _reset_counts()
    path = tmp_path / "sample.raw"
    path.write_bytes(b"")
    existing = CountingParser(path)
    existing.reads.clear()

    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        mass_spectra = MassSpectraBase(
            path,
            sample_name=path.stem,
            spectra_parser=existing,
        )
    assert CountingParser.constructed == 1
    assert mass_spectra.spectra_parser is existing
    assert CountingParser.constructed == 1
    for name in ("sample_name", "analyzer", "instrument_label", "file_location"):
        assert name in existing.reads
    mass_spectra.close()
    assert existing.closed == 0
    assert mass_spectra.spectra_parser is existing


def test_raw_file_location_change_closes_cached_parser_once(tmp_path):
    _reset_counts()
    path = tmp_path / "sample.raw"
    other = tmp_path / "other.raw"
    path.write_bytes(b"")
    other.write_bytes(b"")
    existing = CountingParser(path)
    mass_spectra = MassSpectraBase(
        path,
        sample_name=path.stem,
        spectra_parser=existing,
    )

    mass_spectra.raw_file_location = other

    assert existing.closed == 1
    constructed = CountingParser.constructed
    opened = mass_spectra.spectra_parser
    assert CountingParser.constructed == constructed + 1
    assert opened is not existing
    assert opened.file_location == other
    again = mass_spectra.spectra_parser
    assert again is opened
    assert CountingParser.constructed == constructed + 1


def test_pickle_omits_parser_and_reopens_after_load(tmp_path):
    _reset_counts()
    path = tmp_path / "sample.raw"
    path.write_bytes(b"")
    existing = CountingParser(path)
    mass_spectra = MassSpectraBase(
        path,
        sample_name=path.stem,
        spectra_parser=existing,
    )

    restored = pickle.loads(pickle.dumps(mass_spectra))

    assert mass_spectra.spectra_parser is existing
    assert existing.closed == 0
    assert restored._spectra_parser is None
    constructed = CountingParser.constructed
    reopened = restored.spectra_parser
    assert CountingParser.constructed == constructed + 1
    assert reopened is not existing
    assert reopened.file_location == path


def test_owned_parser_closes_when_mass_spectra_is_collected(tmp_path):
    _reset_counts()
    path = tmp_path / "sample.raw"
    path.write_bytes(b"")

    def _open_and_drop():
        mass_spectra = MassSpectraBase(path)
        mass_spectra.spectra_parser_class = CountingParser
        mass_spectra.spectra_parser
        del mass_spectra
        gc.collect()

    _open_and_drop()
    gc.collect()
    assert CountingParser.constructed == 1
    assert CountingParser.closed == 1


def test_get_average_mass_spectrum_uses_one_parser():
    class Host(LCCalculations):
        def __init__(self):
            self.polarity = "positive"
            self.file_location = "average-test"
            self.builds = 0

        @property
        def spectra_parser(self):
            self.builds += 1
            return self

        def get_mass_spectrum_from_scan(self, scan, spectrum_mode="profile", auto_process=False):
            return SimpleNamespace(
                _mz_exp=np.array([100.0, 101.0]),
                _abundance=np.array([10.0, 5.0]),
            )

        def check_if_grid(self, data):
            return True

    host = Host()
    spectrum = host.get_average_mass_spectrum(
        [1, 2, 3],
        apex_scan=2,
        use_parser=True,
        auto_process=False,
    )
    assert host.builds == 1
    assert spectrum is not None


def test_load_raw_data_checks_parser_class_once():
    calls = []

    class ImportMassSpectraThermoMSFileReader:
        def get_ms_raw(self, spectra, scan_df, time_range):
            return {
                spectra: pd.DataFrame(
                    {"scan": [1], "mz": [100.0], "intensity": [10.0]}
                )
            }

    class Sample:
        def __init__(self):
            self._ms_unprocessed = {}
            self.scan_df = pd.DataFrame({"scan": [1]})
            self.spectra_parser_class = ImportMassSpectraThermoMSFileReader

        @property
        def spectra_parser(self):
            calls.append("open")
            return ImportMassSpectraThermoMSFileReader()

    sample = Sample()

    class Collection:
        samples = [sample]

        def __getitem__(self, index):
            return sample

    LCMSCollection.load_raw_data(Collection(), 0)
    assert calls == ["open"]
    assert 1 in sample._ms_unprocessed


def test_hdf_creation_time_reads_parser_once(tmp_path):
    calls = []

    class Spectra:
        def __iter__(self):
            return iter(())

        file_location = SimpleNamespace(name="a.raw", _str=str(tmp_path / "a.raw"))
        analyzer = "Unknown"
        instrument_label = "Unknown"
        sample_name = "a"
        polarity = "negative"
        spectra_parser_class = CountingParser

        @property
        def spectra_parser(self):
            calls.append("open")
            return SimpleNamespace(get_creation_time=lambda: datetime(2018, 1, 2))

    exporter = HighResMassSpectraExport(str(tmp_path / "out"), Spectra())
    exporter.to_hdf()
    assert calls == ["open"]


def test_close_file_disposes_both_handles_once():
    from corems.mass_spectra.input.rawFileReader import ThermoBaseClass

    events = []
    reader = ThermoBaseClass.__new__(ThermoBaseClass)
    reader.iRawDataPlus = _Disposable(events)
    reader.iFileHeader = _Disposable(events)

    reader.close_file()
    reader.close_file()

    assert events == ["dispose", "dispose"]
    assert reader.iRawDataPlus is None
    assert reader.iFileHeader is None


def test_close_file_skips_header_without_dispose():
    from corems.mass_spectra.input.rawFileReader import ThermoBaseClass

    events = []
    reader = ThermoBaseClass.__new__(ThermoBaseClass)
    reader.iRawDataPlus = _Disposable(events)
    reader.iFileHeader = _NoDispose()

    reader.close_file()

    assert events == ["dispose"]
    assert reader.iRawDataPlus is None
    assert reader.iFileHeader is None


def test_dropping_thermo_reader_closes_once():
    from corems.mass_spectra.input.rawFileReader import ThermoBaseClass

    events = []

    def _drop():
        reader = ThermoBaseClass.__new__(ThermoBaseClass)
        reader.iRawDataPlus = _Disposable(events)
        reader.iFileHeader = _Disposable(events)
        del reader
        gc.collect()

    _drop()
    gc.collect()
    assert events == ["dispose", "dispose"]
