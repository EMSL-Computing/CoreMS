import copy
import os
import re

import numpy as np
import requests
from ms_entropy import FlashEntropySearch

from corems.molecular_id.factory.EI_SQL import EI_LowRes_SQLite, Metadatar

# Globals
METABREF_GCMS_LIBRARY_URL = "https://metabref.emsl.pnnl.gov/api/mslevel/1"
METABREF_DI_MS2_LIBRARY_URL = "https://metabref.emsl.pnnl.gov/api/mslevel/2"
METABREF_FAMES_URL = "https://metabref.emsl.pnnl.gov/api/fames"
METABREF_MOLECULARDATA_URL = "https://metabref.emsl.pnnl.gov/moleculardata/{}"


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

    # All below key:value lookups are based on CoreMS class definitions
    # NOT MetabRef content. For example, MetabRef has keys for PubChem,
    # USI, etc. that are not considered below.

    # Dictionary to map metabref keys to corems keys
    metadatar_cols = {'casno': 'cas',
                      'inchikey': 'inchikey',
                      'inchi': 'inchi',
                      'chebi': 'chebi',
                      'smiles': 'smiles',
                      'kegg': 'kegg',
                      'iupac_name': 'iupac_name',
                      'traditional_name': 'traditional_name', # Not present in metabref
                      'common_name': 'common_name' # Not present in metabref
                     }

    # Dictionary to map metabref keys to corems keys
    lowres_ei_compound_cols = {'id': 'metabref_id',
                               'molecule_name': 'name', # Is this correct?
                               'classify': 'classify', # Not present in metabref
                               'formula': 'formula',
                               'ri': 'ri',
                               'rt': 'retention_time',
                               'source': 'source', # Not present in metabref
                               'casno': 'casno',
                               'comments': 'comment',
                               'source_temp_c': 'source_temp_c', # Not present in metabref
                               'ev': 'ev', # Not present in metabref
                               'peak_count': 'peaks_count',
                               'mz': 'mz',
                               'abundance': 'abundance'
                              }

    # Local result container
    corems_lib = []

    # Enumerate spectra
    for i, source_ in enumerate(metabref_lib):
        # Copy source to prevent modification
        source = source_.copy()
        
        # Flatten source dict
        source =  source.pop('spectrum_data') | source

        # Parse target data
        target = {lowres_ei_compound_cols[k]: v for k, v in source.items() if k in lowres_ei_compound_cols}
        
        # Explicitly add this to connect with LowResCompoundRef later 
        target['rt'] = source['rt']

        # Parse (mz, abundance)
        arr = metabref_spectrum_to_array(target['mz'],
                                         normalize=normalize)
        target['mz'] = arr[:, 0]
        target['abundance'] = arr[:, 1]

        # Parse meta data
        target['metadata'] = {metadatar_cols[k]: v for k, v in source.items() if k in metadatar_cols}

        # Add anything else
        for k in source:
            if k not in lowres_ei_compound_cols:
                target[k] = source[k]

        # Add to CoreMS list
        corems_lib.append(target)

    return corems_lib


def corems_to_sqlite(corems_lib, url='sqlite://'):
    """
    Convert CoreMS-formatted library to SQLite database for local ingestion.

    Parameters
    ----------
    corems_lib : dict
        CoreMS GC-MS library.
    url : str
        URL to SQLite prefix.

    Returns
    -------
    sqlite database
        Spectra contained in SQLite database.

    """
    # Dictionary to map corems keys to all-caps keys
    capped_cols = {'name': 'NAME',
                   'formula': 'FORM',
                   'ri': 'RI',
                   'retention_time': 'RT',
                   'source': 'SOURCE',
                   'casno': 'CASNO',
                   'comment': 'COMMENT',
                   'peaks_count': 'NUM PEAKS'
                  }

    # Initialize SQLite object
    sqlite_obj = EI_LowRes_SQLite(url=url)

    # Iterate spectra
    for data_dict_ in corems_lib:
        # Copy source to prevent modification
        data_dict = data_dict_.copy()
        
        # Add missing capped values
        for k, v in capped_cols.items():
            # Key exists
            if k in data_dict:
                # # This will replace the key
                # data_dict[v] = data_dict.pop(k)

                # This will keep both keys
                data_dict[v] = data_dict[k]

        # Parse number of peaks
        if not data_dict.get('NUM PEAKS'):
            data_dict['NUM PEAKS'] = len(data_dict.get('mz'))

        # Parse CAS number
        if not data_dict.get('CASNO'):
            data_dict['CASNO'] = data_dict.get('CAS')

        if not data_dict['CASNO']:
            data_dict['CASNO'] = 0

        # Build linked metadata table
        if 'metadata' in data_dict:
            if len(data_dict['metadata']) > 0:
                data_dict['metadatar'] = Metadatar(**data_dict.pop('metadata'))
            else:
                data_dict.pop('metadata')

        # Attempt addition to sqlite
        try:
            sqlite_obj.add_compound(data_dict)
        except:
            print(data_dict['NAME'])
    
    return sqlite_obj
