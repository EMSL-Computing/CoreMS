import pytest
import os
import shutil
from pathlib import Path
from urllib.request import urlopen

from corems.transient.input.brukerSolarix import ReadBrukerSolarix
from corems.encapsulation.factory.parameters import MSParameters
from corems.mass_spectra.input.rawFileReader import ImportMassSpectraThermoMSFileReader


LIPIDOMICS_SQLITE_URL = (
    "https://nmdcdemo.emsl.pnnl.gov/minio/lipidomics/parameter_files/"
    "202412_lipid_ref.sqlite"
)


def pytest_addoption(parser):
    parser.addoption(
        "--skip-lipidomics-db",
        action="store_true",
        default=False,
        help="Skip tests that require the lipidomics sqlite library.",
    )


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "lipidomics_db: mark test as requiring the lipidomics sqlite library",
    )


def pytest_collection_modifyitems(config, items):
    if not config.getoption("--skip-lipidomics-db"):
        return

    skip_marker = pytest.mark.skip(reason="skipped by --skip-lipidomics-db")
    for item in items:
        if "lipidomics_db" in item.keywords:
            item.add_marker(skip_marker)


def _download_lipidomics_db(destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(LIPIDOMICS_SQLITE_URL, timeout=300) as response:
        with open(destination, "wb") as out_file:
            shutil.copyfileobj(response, out_file)


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


@pytest.fixture(scope="session")
def lipidomics_sqlite_path(pytestconfig):
    """Returns a local sqlite path for lipidomics library searches.

    The fixture auto-downloads the sqlite file for local runs unless
    disabled by --skip-lipidomics-db or COREMS_LIPIDOMICS_AUTO_DOWNLOAD=0.
    """

    env_path = os.getenv("COREMS_LIPIDOMICS_SQLITE_PATH")
    sqlite_path = (
        Path(env_path).expanduser()
        if env_path
        else Path.cwd() / "tests/tests_data/lcms/202412_lipid_ref.sqlite"
    )

    if sqlite_path.exists():
        return sqlite_path

    if pytestconfig.getoption("--skip-lipidomics-db"):
        pytest.skip("lipidomics sqlite database unavailable and skip option set")

    auto_download = os.getenv("COREMS_LIPIDOMICS_AUTO_DOWNLOAD", "1").lower()
    if auto_download in {"0", "false", "no"}:
        pytest.skip(
            "lipidomics sqlite database missing and auto-download disabled. "
            f"Download from {LIPIDOMICS_SQLITE_URL}"
        )

    try:
        _download_lipidomics_db(sqlite_path)
    except Exception as exc:
        if os.getenv("CI"):
            raise RuntimeError(
                "Failed to download lipidomics sqlite database for CI run"
            ) from exc
        pytest.skip(
            "failed to download lipidomics sqlite database. "
            f"Download manually from {LIPIDOMICS_SQLITE_URL}. Error: {exc}"
        )

    return sqlite_path
