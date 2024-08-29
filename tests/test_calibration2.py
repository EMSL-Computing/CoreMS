import sys
from pathlib import Path
import pytest

sys.path.append(".")

from corems.encapsulation.factory.parameters import MSParameters
from corems.mass_spectrum.calc.AutoRecalibration import HighResRecalibration
from corems.mass_spectrum.calc.Calibration import MzDomainCalibration
from corems.mass_spectrum.calc.CalibrationCalc import FreqDomain_Calibration
from corems.mass_spectrum.input.massList import ReadCoremsMasslist
from corems.molecular_id.calc.ClusterFilter import ClusteringFilter
from corems.molecular_id.search.findOxygenPeaks import FindOxygenPeaks
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.transient.input.brukerSolarix import ReadBrukerSolarix


@pytest.fixture
def mass_spectrum_ftms():
    """Creates a mass spectrum object to be used in the tests"""

    file_location = Path.cwd() / "tests/tests_data/ftms/ESI_NEG_SRFA.d/"
    bruker_reader = ReadBrukerSolarix(file_location)
    bruker_transient = bruker_reader.get_transient()

    # Instantiate the mass spectrum object
    mass_spectrum = bruker_transient.get_mass_spectrum(
        plot_result=False, auto_process=False, keep_profile=True
    )

    mass_spectrum.settings.noise_threshold_method = "log"
    mass_spectrum.settings.noise_threshold_log_nsigma = 12
    mass_spectrum.mspeaks_settings.peak_min_prominence_percent = 0.01

    # Process the mass spectrum
    mass_spectrum.process_mass_spec()

    return mass_spectrum

@pytest.fixture
def ref_file_location():
    return Path.cwd() / "tests/tests_data/ftms/SRFA.ref"

def test_mz_domain_calibration(mass_spectrum_ftms, ref_file_location):
    print("test_mz_domain_calibration")
    mass_spectrum_ftms.settings.min_calib_ppm_error = -10
    mass_spectrum_ftms.settings.max_calib_ppm_error = 10
    mass_spectrum_ftms.filter_by_noise_threshold()

    # Check that mass_spectrum_ftms has not been calibrated
    assert set(mass_spectrum_ftms.mz_cal) == {None}

    MzDomainCalibration(mass_spectrum_ftms, ref_file_location).run()

    # Check that the calibration was successful
    assert mass_spectrum_ftms.calibration_RMS < 0.6


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
    assert mass_spectrum_ftms.calibration_RMS < 0.6

def test_segmentedmzcalibration(mass_spectrum_ftms, ref_file_location):
    # Tests profile mode recalibration
    mass_spectrum_ftms.filter_by_noise_threshold()

    # Check that mass_spectrum_ftms has not been calibrated
    assert set(mass_spectrum_ftms.mz_cal) == {None}

    MzDomainCalibration(mass_spectrum_ftms, ref_file_location, mzsegment=(0, 300)).run()

    # Check that the calibration was successful
    assert mass_spectrum_ftms.calibration_RMS < 0.6

def test_old_calibration(mass_spectrum_ftms):
    usedatoms = {'C': (1,100) , 'H': (4,200), 'O': (1,10)}

    mass_spectrum_ftms.molecular_search_settings.url_database = ''
    mass_spectrum_ftms.molecular_search_settings.error_method = 'None'
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
