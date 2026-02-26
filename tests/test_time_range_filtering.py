"""
Test time range filtering functionality for LC-MS parsers.

This module tests the time_range parameter across different parser implementations
to ensure efficient loading of targeted retention time windows.
"""
from pathlib import Path
import pytest
import os
import tempfile
import shutil
import time

from corems.mass_spectra.input.mzml import MZMLSpectraParser
from corems.mass_spectra.input import rawFileReader
from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectra
from corems.mass_spectra.output.export import LCMSExport
from corems.encapsulation.factory.parameters import LCMSParameters


# Module-level fixtures
@pytest.fixture(scope="module")
def mzml_file():
    """Path to test mzML file."""
    return (
        Path.cwd()
        / "tests/tests_data/lcms/"
        / "test_centroid_neg_RP_metab.mzML"
    )


@pytest.fixture(scope="module")
def thermo_file():
    """Path to test Thermo RAW file."""
    return (
        Path.cwd()
        / "tests/tests_data/lcms/"
        / "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.raw"
    )


@pytest.fixture
def mzml_parser(mzml_file):
    """Instantiate mzML parser."""
    return MZMLSpectraParser(mzml_file)


@pytest.fixture
def thermo_parser(thermo_file):
    """Instantiate Thermo parser."""
    return rawFileReader.ImportMassSpectraThermoMSFileReader(thermo_file)


@pytest.fixture
def hdf5_parser(mzml_file):
    """Create temporary HDF5 file and return parser."""
    temp_dir = tempfile.mkdtemp()
    base_name = "test_time_range"
    base_path = os.path.join(temp_dir, base_name)
    
    # Parse mzML and create HDF5 file
    parser = MZMLSpectraParser(mzml_file)
    lcms_obj = parser.get_lcms_obj(spectra="ms1")
    lcms_obj.parameters = LCMSParameters(use_defaults=True)
    
    # Save to HDF5 using LCMSExport
    exporter = LCMSExport(base_path, lcms_obj)
    exporter.to_hdf(overwrite=True)
    
    # The actual HDF5 file will be at base_path.corems/base_name.hdf5
    hdf5_path = os.path.join(temp_dir, f"{base_name}.corems", f"{base_name}.hdf5")
    
    hdf5_parser = ReadCoreMSHDFMassSpectra(hdf5_path)
    
    yield hdf5_parser
    
    # Cleanup
    try:
        corems_dir = os.path.join(temp_dir, f"{base_name}.corems")
        if os.path.exists(corems_dir):
            shutil.rmtree(corems_dir)
        os.rmdir(temp_dir)
    except Exception:
        pass


@pytest.fixture(params=['mzml', 'thermo', 'hdf5'])
def parser(request, mzml_parser, thermo_parser, hdf5_parser):
    """Parametrized fixture that provides all three parser types."""
    if request.param == 'mzml':
        return mzml_parser
    elif request.param == 'thermo':
        return thermo_parser
    else:
        return hdf5_parser


class TestTimeRangeFiltering:
    """Test time range filtering across all parser implementations."""
    
    def test_get_scans_in_time_range_single_range(self, parser):
        """Test getting scans within a single time range."""
        all_scan_df = parser.get_scan_df()
        time_range = (1.0, 2.0)
        scans_in_range = parser.get_scans_in_time_range(time_range)
        
        assert isinstance(scans_in_range, list)
        assert len(scans_in_range) > 0
        
        # Verify all scans are in range
        for scan_num in scans_in_range:
            scan_time = all_scan_df[all_scan_df.scan == scan_num].scan_time.values[0]
            assert time_range[0] <= scan_time <= time_range[1]
    
    def test_get_scans_in_time_range_multiple_ranges(self, parser):
        """Test getting scans within multiple time ranges."""
        all_scan_df = parser.get_scan_df()
        time_ranges = [(0.5, 1.5), (3.0, 4.0)]
        scans_in_range = parser.get_scans_in_time_range(time_ranges)
        
        assert isinstance(scans_in_range, list)
        assert len(scans_in_range) > 0
        
        # Verify all scans are in at least one range
        for scan_num in scans_in_range:
            scan_time = all_scan_df[all_scan_df.scan == scan_num].scan_time.values[0]
            assert any(start <= scan_time <= end for start, end in time_ranges)
    
    def test_get_scans_in_time_range_with_ms_level(self, parser):
        """Test filtering by both time range and MS level."""
        all_scan_df = parser.get_scan_df()
        time_range = (1.0, 3.0)
        ms1_scans = parser.get_scans_in_time_range(time_range, ms_level=1)
        
        assert len(ms1_scans) > 0
        
        # Verify all returned scans are MS1 and in time range
        for scan_num in ms1_scans:
            scan_info = all_scan_df[all_scan_df.scan == scan_num].iloc[0]
            assert scan_info.ms_level == 1
            assert time_range[0] <= scan_info.scan_time <= time_range[1]
    
    def test_get_scan_df_with_time_range(self, parser):
        """Test get_scan_df with time range filtering."""
        all_scan_df = parser.get_scan_df()
        time_range = (1.0, 2.0)
        filtered_scan_df = parser.get_scan_df(time_range=time_range)
        
        assert len(filtered_scan_df) < len(all_scan_df)
        assert len(filtered_scan_df) > 0
        assert all(time_range[0] <= t <= time_range[1] for t in filtered_scan_df.scan_time)
    
    def test_get_scan_df_with_multiple_time_ranges(self, parser):
        """Test get_scan_df with multiple time ranges."""
        time_ranges = [(0.5, 1.0), (2.0, 2.5)]
        filtered_scan_df = parser.get_scan_df(time_range=time_ranges)
        
        for scan_time in filtered_scan_df.scan_time:
            assert any(start <= scan_time <= end for start, end in time_ranges)
    
    def test_edge_cases(self, parser):
        """Test edge cases for time range filtering."""
        all_scan_df = parser.get_scan_df()
        min_time = all_scan_df.scan_time.min()
        max_time = all_scan_df.scan_time.max()
        
        # Range covering all data
        scans = parser.get_scans_in_time_range((0, 999))
        assert len(scans) == len(all_scan_df)
        
        # Range with no data
        scans = parser.get_scans_in_time_range((999, 1000))
        assert len(scans) == 0
        
        # Range at exact boundaries
        scans = parser.get_scans_in_time_range((min_time, max_time))
        assert len(scans) == len(all_scan_df)
    
    def test_normalize_time_range_helper(self):
        """Test the _normalize_time_range static helper method."""
        from corems.mass_spectra.input.parserbase import SpectraParserInterface
        
        assert SpectraParserInterface._normalize_time_range((1.0, 2.0)) == [(1.0, 2.0)]
        assert SpectraParserInterface._normalize_time_range([(1.0, 2.0), (3.0, 4.0)]) == [(1.0, 2.0), (3.0, 4.0)]
        assert SpectraParserInterface._normalize_time_range(None) is None


class TestParserSpecificFeatures:
    """Test parser-specific features and workflows."""
    
    def test_mzml_lcms_obj_with_time_range(self, mzml_parser):
        """Test loading mzML LCMS object with time range filtering."""
        all_scan_df = mzml_parser.get_scan_df()
        time_range = (1.0, 2.0)
        
        filtered_lcms = mzml_parser.get_lcms_obj(spectra="ms1", time_range=time_range)
        
        assert len(filtered_lcms.scan_df) < len(all_scan_df)
        assert len(filtered_lcms.scan_df) > 0
        assert all(time_range[0] <= t <= time_range[1] for t in filtered_lcms.scan_df.scan_time)
        assert len(filtered_lcms._scans_number_list) == len(filtered_lcms.scan_df)
        assert len(filtered_lcms._retention_time_list) == len(filtered_lcms.scan_df)
    
    def test_thermo_lcms_obj_with_time_range(self, thermo_parser):
        """Test loading Thermo LCMS object with time range filtering."""
        all_scan_df = thermo_parser.get_scan_df()
        time_range = (0.5, 1.0)
        
        filtered_lcms = thermo_parser.get_lcms_obj(spectra="ms1", time_range=time_range)
        
        assert len(filtered_lcms.scan_df) < len(all_scan_df)
        assert len(filtered_lcms.scan_df) > 0
        assert all(time_range[0] <= t <= time_range[1] for t in filtered_lcms.scan_df.scan_time)
        assert len(filtered_lcms._scans_number_list) == len(filtered_lcms.scan_df)
        assert len(filtered_lcms._retention_time_list) == len(filtered_lcms.scan_df)
    
    def test_hdf5_lcms_obj_with_time_range(self, hdf5_parser):
        """Test HDF5 LCMS object accepts time_range parameter."""
        all_scan_df = hdf5_parser.get_scan_df()
        time_range = (1.0, 2.0)
        
        # HDF5 parser accepts time_range for interface consistency
        _ = hdf5_parser.get_lcms_obj(
            load_raw=True,
            load_light=False,
            use_original_parser=False,
            time_range=time_range
        )
        
        # Verify scan_df can be filtered
        filtered_scan_df = hdf5_parser.get_scan_df(time_range=time_range)
        assert len(filtered_scan_df) < len(all_scan_df)
        assert len(filtered_scan_df) > 0
    
    def test_thermo_performance_benchmark(self, thermo_parser):
        """Test that time range filtering provides measurable performance improvement."""
        target_rt = 5.0
        rt_window = 2.0
        time_range = (target_rt - rt_window, target_rt + rt_window)
        
        # Time filtered load
        start_filtered = time.time()
        filtered_lcms = thermo_parser.get_lcms_obj(spectra="ms1", time_range=time_range)
        filtered_lcms.parameters = LCMSParameters(use_defaults=True)
        filtered_lcms.parameters.lc_ms.peak_picking_method = "persistent homology"
        filtered_lcms.parameters.lc_ms.ph_inten_min_rel = 0.01
        filtered_lcms.parameters.lc_ms.ph_persis_min_rel = 0.05
        filtered_lcms.find_mass_features()
        time_filtered = time.time() - start_filtered
        
        # Time full load
        start_full = time.time()
        full_lcms = thermo_parser.get_lcms_obj(spectra="ms1")
        full_lcms.parameters = LCMSParameters(use_defaults=True)
        full_lcms.parameters.lc_ms.peak_picking_method = "persistent homology"
        full_lcms.parameters.lc_ms.ph_inten_min_rel = 0.01
        full_lcms.parameters.lc_ms.ph_persis_min_rel = 0.05
        full_lcms.find_mass_features()
        time_full = time.time() - start_full
        
        speedup = time_full / time_filtered if time_filtered > 0 else 0
        
        print(f"\n{'='*60}")
        print("Time Range Filtering Performance Test")
        print(f"{'='*60}")
        print(f"Time range: {time_range[0]}-{time_range[1]} minutes")
        print(f"Filtered load: {len(filtered_lcms.scan_df)} scans in {time_filtered:.2f}s")
        print(f"Full load: {len(full_lcms.scan_df)} scans in {time_full:.2f}s")
        print(f"Speedup: {speedup:.2f}x")
        print(f"{'='*60}")
        
        assert len(filtered_lcms.scan_df) < len(full_lcms.scan_df)
        assert time_filtered < time_full
        assert speedup >= 1.2, f"Expected at least 1.2x speedup, got {speedup:.2f}x"


class TestMassFeatureIntegration:
    """Test time range filtering integration with mass feature detection."""
    
    def test_targeted_search_with_time_range(self, mzml_file):
        """Test targeted search with time-filtered data."""
        parser = MZMLSpectraParser(mzml_file)
        target_time_range = (1.0, 2.0)
        
        lcms_obj = parser.get_lcms_obj(spectra="ms1", time_range=target_time_range)
        lcms_obj.parameters = LCMSParameters(use_defaults=True)
        lcms_obj.parameters.lc_ms.peak_picking_method = "centroided_persistent_homology"
        lcms_obj.parameters.mass_spectrum["ms1"].mass_spectrum.noise_threshold_method = "relative_abundance"
        lcms_obj.find_mass_features()
        
        assert len(lcms_obj.mass_features) > 0
        
        # Verify all features are within target time range
        for mf_id, mf in lcms_obj.mass_features.items():
            assert target_time_range[0] <= mf.retention_time <= target_time_range[1]
    
    def test_multiple_time_windows_for_standards(self, mzml_file):
        """Test loading multiple time windows for internal standards."""
        parser = MZMLSpectraParser(mzml_file)
        standard_windows = [(0.3, 0.7), (1.8, 2.2), (3.3, 3.7)]
        
        lcms_obj = parser.get_lcms_obj(spectra="ms1", time_range=standard_windows)
        lcms_obj.parameters = LCMSParameters(use_defaults=True)
        lcms_obj.parameters.lc_ms.peak_picking_method = "centroided_persistent_homology"
        lcms_obj.parameters.mass_spectrum["ms1"].mass_spectrum.noise_threshold_method = "relative_abundance"
        lcms_obj.find_mass_features()
        
        # Verify scan times are only from specified windows
        for scan_time in lcms_obj.scan_df.scan_time:
            assert any(start <= scan_time <= end for start, end in standard_windows)


