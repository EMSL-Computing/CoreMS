import pytest
from pathlib import Path

from corems.transient.input.brukerSolarix import ReadBrukerSolarix
from corems.encapsulation.factory.parameters import MSParameters
from corems.mass_spectra.input.rawFileReader import ImportMassSpectraThermoMSFileReader


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


@pytest.fixture
def lcms_obj():
    """Returns an LCMS object for the tests"""
    file_raw = (
        Path.cwd()
        / "tests/tests_data/lcms/"
        / "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.raw"
    )
    parser = ImportMassSpectraThermoMSFileReader(file_raw)
    instrument_info = parser.get_instrument_info()
    assert instrument_info['model'] == "Orbitrap Velos Pro"
    creation_time = parser.get_creation_time()
    assert creation_time.year == 2018

    # Instatiate lc-ms data object using parser and pull in ms1 spectra into dataframe (without storing as MassSpectrum objects to save memory)
    myLCMSobj = parser.get_lcms_obj(spectra="ms1")

    return myLCMSobj


@pytest.fixture
def msp_file_location():
    """Returns the location of the msp file for the tests"""
    return Path.cwd() / "tests/tests_data/lcms/test_db.msp"


@pytest.fixture
def postgres_database():
    """Returns the location of the postgres database for the tests"""
    # Change this if running locally or the DB is running in a different location
    # Try to use postgres first (for CI/CD), fall back to sqlite3 (local)
    import socket
    try:
        # Test if postgres hostname is reachable
        socket.gethostbyname('postgres')
        return "postgresql://coremsdb:coremsmolform@postgres:5432/molformula"  ## Git CI/CD Build Pipeline
    except socket.gaierror:
        # Fall back to sqlite3 for local testing when postgres is not available
        return "" ## sqlite3 database (local)
