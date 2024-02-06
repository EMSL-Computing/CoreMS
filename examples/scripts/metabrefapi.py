import requests
import os
import numpy as np
import re
from ms_entropy import FlashEntropySearch
import copy
from corems.molecular_id.factory.EI_SQL import EI_LowRes_SQLite


# Globals
METABREF_GCMS_LIBRARY_URL = "https://metabref.emsl.pnnl.gov/api/mslevel/1"
METABREF_DI_MS2_LIBRARY_URL = "https://metabref.emsl.pnnl.gov/api/mslevel/2"
METABREF_FAMES_URL = "https://metabref.emsl.pnnl.gov/api/fames"


def set_token(path):
    """
    Set environment variable for MetabRef database token.

    Parameters
    ----------
    path : str
        Path to token.

    """

    # Read token from file
    with open(path, "r", encoding="utf-8") as f:
        token = f.readline().strip()

    # Set environment variable
    os.environ["METABREF_TOKEN"] = token


def get_token():
    """
    Get environment variable for MetabRef database token.

    Returns
    -------
    str
        Token string.

    """

    # Check for token
    if "METABREF_TOKEN" not in os.environ:
        raise ValueError("Must set METABREF_TOKEN environment variable.")

    # Get token from environment variables
    return os.environ.get("METABREF_TOKEN")


def get_header():
    """
    Access stored MetabRef database token and prepare as header.

    Returns
    -------
    str
        Header string.

    """

    # Get token
    token = get_token()

    # Pad header information
    header = {"Authorization": f"Bearer {token}",
              "Content-Type": "text/plain"}

    return header


def get_query(url):
    """
    Request payload from URL according to `get` protocol.

    Parameters
    ----------
    url : str
        URL for request.

    Returns
    -------
    dict
        Response as JSON.

    """

    # Query URL via `get`
    response = requests.get(url, headers=get_header())

    # Check response
    response.raise_for_status()

    # Return as JSON
    return response.json()


def post_query(url, variable, values, tolerance):
    """
    Request payload from URL according to `post` protocol.

    Parameters
    ----------
    url : str
        URL for request.

    Returns
    -------
    dict
        Response as JSON.

    """

    # Coerce to string
    if type(variable) is not str:
        variable = str(variable).replace(" ", "")

    if type(values) is not str:
        values = str(values).replace(" ", "")

    if type(tolerance) is not str:
        tolerance = str(tolerance).replace(" ", "")

    # Query URL via `post`
    response = requests.post(os.path.join(url, variable, tolerance),
                             data=values,
                             headers=get_header())

    # Check response
    response.raise_for_status()

    # Return as JSON
    return response.json()


def get_metabref_gcms_library():
    """
    Request MetabRef GC/MS library.

    Returns
    -------
    dict
        Response as JSON.

    """

    return get_query(METABREF_GCMS_LIBRARY_URL)['GC-MS']


def get_metabref_fames_library():
    """
    Request MetabRef GC/MS FAMEs library.

    Returns
    -------
    dict
        Response as JSON.

    """

    return get_query(METABREF_FAMES_URL)['GC-MS']


def get_metabref_di_ms2_library():
    """
    Request MetabRef DI MS2 library.

    Returns
    -------
    dict
        Response as JSON.

    """

    return get_query(METABREF_DI_MS2_LIBRARY_URL)['DI-MS']


def metabref_spectrum_to_array(spectrum, normalize=True):
    """
    Convert MetabRef-formatted spectrum to array.

    Parameters
    ----------
    spectrum : str
        MetabRef spectrum, i.e. list of (m/z,abundance) pairs.
    normalize : bool
        Normalize the spectrum by its magnitude.

    Returns
    -------
    :obj:`~numpy.array`
        Array of shape (N, 2).

    """

    # Convert parenthesis-delimited string to array
    arr = np.array(re.findall(r'\(([^,]+),([^)]+)\)',
                   spectrum), dtype=float).reshape(-1, 2)

    # Normalize the array
    if normalize:
        arr[:, -1] = arr[:, -1] / arr[:, -1].sum()

    return arr


def metabref_to_flashentropy(metabref_lib, normalize=True):
    """
    Convert metabref-formatted library to FlashEntropy library.

    Parameters
    ----------
    metabref_lib : dict
        MetabRef MS2 library in JSON format.
    normalize : bool
        Normalize each spectrum by its magnitude.

    Returns
    -------
    :obj:`~ms_entropy.FlashEntropySearch`
        MS2 library as FlashEntropy search instance.

    """

    # Deep copy CoreMS library
    fe_lib = copy.deepcopy(metabref_lib)

    # Enumerate spectra
    for i, source in enumerate(fe_lib):
        # Reorganize source dict
        spectrum = source['spectrum_data']

        # Convert CoreMS spectrum to array, store as `peaks`
        spectrum['peaks'] = metabref_spectrum_to_array(spectrum.pop('mz'),
                                                       normalize=normalize)

        # Rename for FlashEntropy
        spectrum['precursor_mz'] = spectrum.pop('precursor_ion')

        # Overwrite spectrum
        fe_lib[i] = spectrum

    # Initialize FlashEntropy
    fes = FlashEntropySearch()

    # Build index 
    fes.build_index(fe_lib)

    return fes


def metabref_to_corems(metabref_lib, normalize=True):
    """
    Convert MetabRef-formatted library to CoreMS-formatted library for local ingestion.

    Parameters
    ----------
    metabref_lib : dict
        MetabRef GC-MS library in JSON format.
    normalize : bool
        Normalize each spectrum by its magnitude.

    Returns
    -------
    list of dict
        List of each spectrum contained in dictionary.

    """

    # Known fields
    known_fields = ['ri', 'rt', 'comments', 'peak_count', 'mz', 'usi', 'id']

    # Local result container
    corems_lib = []

    # Enumerate spectra
    for i, source in enumerate(metabref_lib):
        # Reorganize source dict
        spectrum = source['spectrum_data']
        
        # Empty spectrum container
        target = {}
        
        # Convert CoreMS spectrum to array
        arr =  metabref_spectrum_to_array(spectrum['mz'],
                                          normalize=normalize)

        # Store as mz, abundance
        target['mz'] = arr[:, 0]
        target['abundance'] = arr[:, 1]

        # Other known fields
        target['NAME'] = source.get('molecule_name', None)
        target['COMMENT'] = spectrum.get('comments', None)
        target['RI'] = spectrum.get('ri', None)
        target['RT'] = spectrum.get('rt', None)
        target['NUM PEAKS'] = spectrum.get('peak_count', None)
        target['metabref_id'] = spectrum.get('id', None)

        # Extra meta data
        for k in spectrum:
            if k not in known_fields:
                target[k] = spectrum[k]

        # Add to CoreMS list
        corems_lib.append(target)

    return corems_lib


def corems_to_sqlite(corems_lib, url='sqlite://', normalize=True):
    """
    Convert CoreMS-formatted library to SQLite database for local ingestion.

    Parameters
    ----------
    corems_lib : dict
        CoreMS GC-MS library.
    url : str
        URL to SQLite prefix.
    normalize : bool
        Normalize each spectrum by its magnitude.

    Returns
    -------
    sqlite database
        Spectra contained in SQLite database.

    """

    # Initialize SQLite object
    sqlite_obj = EI_LowRes_SQLite(url=url)

    # Iterate spectra
    for data_dict in corems_lib:
        # Parse number of peaks
        if not data_dict.get('NUM PEAKS'):
            data_dict['NUM PEAKS'] = len(data_dict.get('mz'))

        # Parse CAS number
        if not data_dict.get('CASNO'):
            data_dict['CASNO'] = data_dict.get('CAS')

        if not data_dict['CASNO']:
            data_dict['CASNO'] = 0

        # Attempt addition to sqlite
        try:
            sqlite_obj.add_compound(data_dict)
        except:
            print(data_dict['NAME'])
    
    return sqlite_obj
