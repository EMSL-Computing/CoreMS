import pytest
from pathlib import Path

from corems.transient.input.brukerSolarix import ReadBrukerSolarix
from corems.encapsulation.factory.parameters import MSParameters

@pytest.fixture
def mass_spectrum_ftms(bruker_transient):
    """Creates a mass spectrum object to be used in the tests"""
    # Instantiate the mass spectrum object
    mass_spectrum = bruker_transient.get_mass_spectrum(
        plot_result=False, auto_process=False, keep_profile=True
    )
    mass_spectrum.parameters = MSParameters(use_defaults=True)
    # Process the mass spectrum
    mass_spectrum.process_mass_spec()

    return mass_spectrum

@pytest.fixture
def ref_file_location():
    """Returns the location of the reference file for calibration for the tests"""
    return Path.cwd() / "tests/tests_data/ftms/SRFA.ref"

@pytest.fixture
def ftms_file_location():
    """Returns the location of the FTMS file for the tests"""
    return Path.cwd() / "tests/tests_data/ftms/ESI_NEG_SRFA.d/"

@pytest.fixture
def bruker_transient(ftms_file_location):
    """Returns the transient object for the FTMS file"""
    bruker_reader = ReadBrukerSolarix(ftms_file_location)
    bruker_transient = bruker_reader.get_transient()

    return bruker_transient