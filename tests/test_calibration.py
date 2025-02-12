import sys
from pathlib import Path
import pytest

from corems.mass_spectrum.calc.AutoRecalibration import HighResRecalibration
from corems.mass_spectrum.calc.Calibration import MzDomainCalibration
from corems.mass_spectrum.calc.CalibrationCalc import FreqDomain_Calibration
from corems.mass_spectrum.input.massList import ReadCoremsMasslist
from corems.molecular_id.calc.ClusterFilter import ClusteringFilter
from corems.molecular_id.search.findOxygenPeaks import FindOxygenPeaks
from corems.encapsulation.factory.parameters import MSParameters, reset_ms_parameters


@pytest.fixture
def mass_spectrum_centroid():

    file_location = Path.cwd() / "tests/tests_data/ftms/ESI_NEG_SRFA_UnCal_Unassign.csv"
    MSParameters.mass_spectrum.noise_threshold_method = "relative_abundance"
    MSParameters.mass_spectrum.noise_threshold_min_relative_abundance = 0.1

    # load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    mass_list_reader = ReadCoremsMasslist(
        file_location, analyzer="ICR", instrument_label="12T"
    )

    mass_spectrum = mass_list_reader.get_mass_spectrum(loadSettings=False)

    # Return the MSParameters to the default values
    reset_ms_parameters()

    return mass_spectrum


def test_mz_domain_calibration(mass_spectrum_ftms, ref_file_location):
    print("test_mz_domain_calibration")
    mass_spectrum_ftms.settings.min_calib_ppm_error = -10
    mass_spectrum_ftms.settings.max_calib_ppm_error = 10
    mass_spectrum_ftms.filter_by_noise_threshold()

    # Check that mass_spectrum_ftms has not been calibrated
    assert set(mass_spectrum_ftms.mz_cal) == {None}

    MzDomainCalibration(mass_spectrum_ftms, ref_file_location).run()

    # Check that the calibration was successful
    assert mass_spectrum_ftms.calibration_RMS < 2


def test_autorecalibration(mass_spectrum_ftms, ref_file_location):
    mass_spectrum_ftms.filter_by_noise_threshold()

    # Check that mass_spectrum_ftms has not been calibrated
    assert set(mass_spectrum_ftms.mz_cal) == {None}

    auto_error_bounds = HighResRecalibration(
        mass_spectrum_ftms, plot=False, docker=False
    ).determine_error_boundaries()

    mass_spectrum_ftms.settings.min_calib_ppm_error = auto_error_bounds[-1][0]
    mass_spectrum_ftms.settings.max_calib_ppm_error = auto_error_bounds[-1][1]

    MzDomainCalibration(mass_spectrum_ftms, ref_file_location).run()

    # Check that the calibration was successful
    assert mass_spectrum_ftms.calibration_RMS < 2


def test_segmentedmzcalibration(mass_spectrum_ftms, ref_file_location):
    # Tests profile mode recalibration
    mass_spectrum_ftms.filter_by_noise_threshold()
    mass_spectrum_ftms.parameters.mass_spectrum.min_calib_ppm_error = -5
    mass_spectrum_ftms.parameters.mass_spectrum.max_calib_ppm_error = 5

    # Check that mass_spectrum_ftms has not been calibrated
    assert set(mass_spectrum_ftms.mz_cal) == {None}

    MzDomainCalibration(mass_spectrum_ftms, ref_file_location, mzsegment=(0, 300)).run()

    # Check that the calibration was successful
    assert mass_spectrum_ftms.calibration_RMS < 2


def test_old_calibration(mass_spectrum_ftms):
    usedatoms = {"C": (1, 100), "H": (4, 200), "O": (1, 10)}

    mass_spectrum_ftms.molecular_search_settings.url_database = "postgresql://coremsdb:coremsmolform@postgres:5432/molformula"
    mass_spectrum_ftms.molecular_search_settings.error_method = "None"
    mass_spectrum_ftms.molecular_search_settings.min_ppm_error = -5
    mass_spectrum_ftms.molecular_search_settings.max_ppm_error = 5
    mass_spectrum_ftms.molecular_search_settings.mz_error_range = 1
    mass_spectrum_ftms.molecular_search_settings.isProtonated = True
    mass_spectrum_ftms.molecular_search_settings.isRadical = True
    mass_spectrum_ftms.molecular_search_settings.usedAtoms = usedatoms

    # Check that mass_spectrum_ftms has not been calibrated by checking that mz_cal are all None
    assert set(mass_spectrum_ftms.mz_cal) == {None}

    find_formula_thread = FindOxygenPeaks(mass_spectrum_ftms)
    find_formula_thread.run()

    mspeaks_results = find_formula_thread.get_list_found_peaks()

    calibrate = FreqDomain_Calibration(mass_spectrum_ftms, mspeaks_results)
    calibrate.linear()
    calibrate.step_fit()
    calibrate.quadratic(iteration=True)
    calibrate.ledford_calibration()

    # Check that the calibration was successful
    assert set(mass_spectrum_ftms.mz_cal) != {None}

    mass_spectrum_ftms.molecular_search_settings.error_method = "symmetrical"
    mass_spectrum_ftms.molecular_search_settings.min_ppm_error = -3
    mass_spectrum_ftms.molecular_search_settings.max_ppm_error = 3
    mass_spectrum_ftms.molecular_search_settings.mz_error_range = 1
    mass_spectrum_ftms.molecular_search_settings.mz_error_average = 0
    mass_spectrum_ftms.molecular_search_settings.min_abun_error = -30  # percentage
    mass_spectrum_ftms.molecular_search_settings.max_abun_error = 70  # percentage
    mass_spectrum_ftms.molecular_search_settings.isProtonated = True
    mass_spectrum_ftms.molecular_search_settings.isRadical = True

    mass_spectrum_ftms.molecular_search_settings.usedAtoms = {
        "C": (1, 100),
        "H": (4, 200),
        "O": (0, 20),
        "N": (0, 1),
        "S": (0, 0),
        "P": (0, 0),
    }

    ClusteringFilter().filter_kendrick(mass_spectrum_ftms)


def test_mz_domain_calibration_centroid(mass_spectrum_centroid, ref_file_location):
    mass_spectrum_centroid.settings.min_calib_ppm_error = -10
    mass_spectrum_centroid.settings.max_calib_ppm_error = 10
    mass_spectrum_centroid.calib_pol_order = 1

    mass_spectrum_centroid.filter_by_noise_threshold()

    # Check that mass_spectrum_centroid has not been calibrated
    assert set(mass_spectrum_centroid.mz_cal) == {None}

    MzDomainCalibration(mass_spectrum_centroid, ref_file_location).run()

    # check there is an output
    assert mass_spectrum_centroid.calibration_points == 25
    assert round(mass_spectrum_centroid.calibration_RMS, 2) == round(0.591, 2)


def test_auto_calibration_centroid(mass_spectrum_centroid, ref_file_location):
    mass_spectrum_centroid.filter_by_noise_threshold()

    # Check that mass_spectrum_centroid has not been calibrated
    assert set(mass_spectrum_centroid.mz_cal) == {None}

    auto_error_bounds = HighResRecalibration(
        mass_spectrum_centroid, plot=False, docker=False
    ).determine_error_boundaries()

    mass_spectrum_centroid.settings.min_calib_ppm_error = auto_error_bounds[-1][0]
    mass_spectrum_centroid.settings.max_calib_ppm_error = auto_error_bounds[-1][1]

    MzDomainCalibration(mass_spectrum_centroid, ref_file_location).run()

    # Check that the calibration was successful
    assert mass_spectrum_centroid.calibration_RMS < 0.6


def test_segmentedmzcalibration_centroid(mass_spectrum_centroid, ref_file_location):
    mass_spectrum_centroid.filter_by_noise_threshold()
    mass_spectrum_centroid.settings.min_calib_ppm_error = -10
    mass_spectrum_centroid.settings.max_calib_ppm_error = 10

    # Check that mass_spectrum_centroid has not been calibrated
    assert set(mass_spectrum_centroid.mz_cal) == {None}

    MzDomainCalibration(
        mass_spectrum_centroid, ref_file_location, mzsegment=(0, 300)
    ).run()

    # Check that the calibration was successful
    assert mass_spectrum_centroid.calibration_RMS < 0.6
