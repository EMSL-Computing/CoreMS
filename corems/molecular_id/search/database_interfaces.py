import os
import re
from abc import ABC
from io import StringIO
from pathlib import Path

import numpy as np
import requests
import pandas as pd
from ms_entropy import FlashEntropySearch

from corems.molecular_id.factory.EI_SQL import (
    EI_LowRes_SQLite,
    Metadatar,
    MetaboliteMetadata,
)
from corems.molecular_id.factory.lipid_molecular_metadata import LipidMetadata
from corems.mass_spectra.calc.lc_calc import find_closest


class SpectralDatabaseInterface(ABC):
    """
    Base class that facilitates connection to spectral reference databases,
    such as EMSL's Metabolomics Reference Database (MetabRef).

    """

    def __init__(self, key=None):
        """
        Initialize instance.

        Parameters
        ----------
        key : str
            Token key.

        """

        self.key = key

    def set_token(self, path):
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
        os.environ[self.key] = token

    def get_token(self):
        """
        Get environment variable for database token.

        Returns
        -------
        str
            Token string.

        """

        # Check for token
        if self.key not in os.environ:
            raise ValueError("Must set {} environment variable.".format(self.key))

        # Get token from environment variables
        return os.environ.get(self.key)

    def get_header(self):
        """
        Access stored database token and prepare as header.

        Returns
        -------
        str
            Header string.

        """

        # Get token
        token = self.get_token()

        # Pad header information
        header = {"Authorization": f"Bearer {token}", "Content-Type": "text/plain"}

        return header

    def get_query(self, url, use_header=True):
        """
        Request payload from URL according to `get` protocol.

        Parameters
        ----------
        url : str
            URL for request.
        use_header: bool
            Whether or not the query should include the header

        Returns
        -------
        dict
            Response as JSON.

        """

        # Query URL via `get`
        if use_header:
            response = requests.get(url, headers=self.get_header())
        else:
            response = requests.get(url)

        # Check response
        response.raise_for_status()

        # Return as JSON
        return response.json()

    def post_query(self, url, variable, values, tolerance):
        """
        Request payload from URL according to `post` protocol.

        Parameters
        ----------
        url : str
            URL for request.
        variable : str
            Variable to query.
        values : str
            Specific values of `variable` to query.
        tolerance : str
            Query tolerance relative to `values`.

        Returns
        -------
        dict
            Response as JSON.

        """

        # Coerce to string
        if not isinstance(variable, str):
            variable = str(variable).replace(" ", "")

        if not isinstance(values, str):
            values = str(values).replace(" ", "")

        if not isinstance(tolerance, str):
            tolerance = str(tolerance).replace(" ", "")

        # Query URL via `post`
        response = requests.post(
            os.path.join(url, variable, tolerance),
            data=values,
            headers=self.get_header(),
        )

        # Check response
        response.raise_for_status()

        # Return as JSON
        return response.json()

    def _check_flash_entropy_kwargs(self, fe_kwargs):
        """
        Check FlashEntropy keyword arguments.

        Parameters
        ----------
        fe_kwargs : dict
            Keyword arguments for FlashEntropy search.


        Raises
        ------
        ValueError
            If "min_ms2_difference_in_da" or "max_ms2_tolerance_in_da" are present in `fe_kwargs` and they
            are not equal.

        """
        # If "min_ms2_difference_in_da" in fe_kwargs, check that "max_ms2_tolerance_in_da" is also present and that min_ms2_difference_in_da = 2xmax_ms2_tolerance_in_da
        if (
            "min_ms2_difference_in_da" in fe_kwargs
            or "max_ms2_tolerance_in_da" in fe_kwargs
        ):
            if (
                "min_ms2_difference_in_da" not in fe_kwargs
                or "max_ms2_tolerance_in_da" not in fe_kwargs
            ):
                raise ValueError(
                    "Both 'min_ms2_difference_in_da' and 'max_ms2_tolerance_in_da' must be specified."
                )
            if (
                fe_kwargs["min_ms2_difference_in_da"]
                != 2 * fe_kwargs["max_ms2_tolerance_in_da"]
            ):
                raise ValueError(
                    "The values of 'min_ms2_difference_in_da' must be exactly 2x 'max_ms2_tolerance_in_da'."
                )

    def _get_format_func(self, format):
        """
        Obtain format function by key.

        Returns
        -------
        func
            Formatting function.
        """

        if format.lower() in self.format_map.keys():
            return self.format_map[format.lower()]

        raise ValueError(("{} not a supported format.").format(format))

    def _dict_to_dataclass(self, metabref_lib, data_class):
        """
        Convert dictionary to dataclass.

        Notes
        -----
        This function will pull the attributes a dataclass and its parent class
        and convert the dictionary to a dataclass instance with the appropriate
        attributes.

        Parameters
        ----------
        data_class : :obj:`~dataclasses.dataclass`
            Dataclass to convert to.
        metabref_lib : dict
            Metabref dictionary object to convert to dataclass.

        Returns
        -------
        :obj:`~dataclasses.dataclass`
            Dataclass instance.

        """

        # Get list of expected attributes of data_class
        data_class_keys = list(data_class.__annotations__.keys())

        # Does the data_class inherit from another class, if so, get the attributes of the parent class as well
        if len(data_class.__mro__) > 2:
            parent_class_keys = list(data_class.__bases__[0].__annotations__.keys())
            data_class_keys = list(set(data_class_keys + parent_class_keys))

        # Remove keys that are not in the data_class from the input dictionary
        input_dict = {k: v for k, v in metabref_lib.items() if k in data_class_keys}

        # Add keys that are in the data class but not in the input dictionary as None
        for key in data_class_keys:
            if key not in input_dict.keys():
                input_dict[key] = None
        return data_class(**input_dict)

    @staticmethod
    def normalize_peaks(arr):
        """
        Normalize peaks in an array.

        Parameters
        ----------
        arr : :obj:`~numpy.array`
            Array of shape (N, 2), with m/z in the first column and abundance in
            the second.

        Returns
        -------
        :obj:`~numpy.array`
            Normalized array of shape (N, 2), with m/z in the first column and
            normalized abundance in the second.
        """
        # Normalize the array
        arr[:, -1] = arr[:, -1] / arr[:, -1].sum()

        return arr

    @staticmethod
    def _build_flash_entropy_index(fe_lib, fe_kwargs={}, clean_spectra=True):
        """
        Build FlashEntropy index.

        Parameters
        ----------
        fe_lib : list
            List of spectra to build index from. Can be a list of dictionaries or
            a FlashEntropy search instance.
        fe_kwargs : dict, optional
            Keyword arguments for FlashEntropy search.
        clean_spectra : bool, optional
            Clean spectra before building index. Default is True.

        Returns
        -------
        :obj:`~ms_entropy.FlashEntropySearch`
            FlashEntropy search instance.

        """
        # Initialize FlashEntropy
        fe_init_kws = [
            "max_ms2_tolerance_in_da",
            "mz_index_step",
            "low_memory",
            "path_data",
        ]
        fe_init_kws = {k: v for k, v in fe_kwargs.items() if k in fe_init_kws}
        fes = FlashEntropySearch(**fe_init_kws)

        # Build FlashEntropy index
        fe_index_kws = [
            "max_indexed_mz",
            "precursor_ions_removal_da",
            "noise_threshold",
            "min_ms2_difference_in_da",
            "max_peak_num",
        ]
        fe_index_kws = {k: v for k, v in fe_kwargs.items() if k in fe_index_kws}
        fes.build_index(fe_lib, **fe_index_kws, clean_spectra=clean_spectra)

        return fes


class MetabRefInterface(SpectralDatabaseInterface):
    """
    Interface to the Metabolomics Reference Database.
    """

    def __init__(self):
        """
        Initialize instance.

        """

        super().__init__(key=None)

    def spectrum_to_array(self, spectrum, normalize=True):
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
            Array of shape (N, 2), with m/z in the first column and abundance in
            the second.

        """

        # Convert parenthesis-delimited string to array
        arr = np.array(
            re.findall(r"\(([^,]+),([^)]+)\)", spectrum), dtype=float
        ).reshape(-1, 2)

        if normalize:
            arr = self.normalize_peaks(arr)

        return arr

    def _to_flashentropy(self, metabref_lib, normalize=True, fe_kwargs={}):
        """
        Convert metabref-formatted library to FlashEntropy library.

        Parameters
        ----------
        metabref_lib : dict
            MetabRef MS2 library in JSON format or FlashEntropy search instance (for reformatting at different MS2 separation).
        normalize : bool
            Normalize each spectrum by its magnitude.
        fe_kwargs : dict, optional
            Keyword arguments for instantiation of FlashEntropy search and building index for FlashEntropy search;
            any keys not recognized will be ignored. By default, all parameters set to defaults.

        Returns
        -------
        :obj:`~ms_entropy.FlashEntropySearch`
            MS2 library as FlashEntropy search instance.

        Raises
        ------
        ValueError
            If "min_ms2_difference_in_da" or "max_ms2_tolerance_in_da" are present in `fe_kwargs` and they are not equal.

        """
        self._check_flash_entropy_kwargs(fe_kwargs)

        # Initialize empty library
        fe_lib = []

        # Enumerate spectra
        for i, source in enumerate(metabref_lib):
            # Reorganize source dict, if necessary
            if "spectrum_data" in source.keys():
                spectrum = source["spectrum_data"]
            else:
                spectrum = source

            # Rename precursor_mz key for FlashEntropy
            if "precursor_mz" not in spectrum.keys():
                spectrum["precursor_mz"] = spectrum.pop("precursor_ion")

            # Convert CoreMS spectrum to array and clean, store as `peaks`
            spectrum["peaks"] = self.spectrum_to_array(
                spectrum["mz"], normalize=normalize
            )

            # Add spectrum to library
            fe_lib.append(spectrum)

        # Build FlashEntropy index
        fe_search = self._build_flash_entropy_index(fe_lib, fe_kwargs=fe_kwargs)

        return fe_search

    def get_query(self, url, use_header=False):
        """Overwrites the get_query method on the parent class to default to not use a header

        Notes
        -----
        As of January 2025, the metabref database no longer requires a token and therefore no header is needed

        """
        return super().get_query(url, use_header)


class MetabRefGCInterface(MetabRefInterface):
    """
    Interface to the Metabolomics Reference Database.
    """

    def __init__(self):
        """
        Initialize instance.

        """

        super().__init__()
        self.GCMS_LIBRARY_URL = "https://metabref.emsl.pnnl.gov/api/mslevel/1"
        self.FAMES_URL = "https://metabref.emsl.pnnl.gov/api/fames"

        self.__init_format_map__()

    def __init_format_map__(self):
        """
        Initialize database format mapper, enabling multiple format requests.

        """

        # Define format workflows
        self.format_map = {
            "json": lambda x, normalize, fe_kwargs: x,
            "dict": lambda x,
            normalize,
            fe_kwargs: self._to_LowResolutionEICompound_dict(x, normalize),
            "sql": lambda x,
            normalize,
            fe_kwargs: self._LowResolutionEICompound_dict_to_sqlite(
                self._to_LowResolutionEICompound_dict(x, normalize)
            ),
        }

        # Add aliases
        self.format_map["metabref"] = self.format_map["json"]
        self.format_map["datadict"] = self.format_map["dict"]
        self.format_map["data-dict"] = self.format_map["dict"]
        self.format_map["lowreseicompound"] = self.format_map["dict"]
        self.format_map["lowres"] = self.format_map["dict"]
        self.format_map["lowresgc"] = self.format_map["dict"]
        self.format_map["sqlite"] = self.format_map["sql"]

    def available_formats(self):
        """
        View list of available formats.

        Returns
        -------
        list
            Format map keys.
        """

        return list(self.format_map.keys())

    def get_library(self, format="json", normalize=False):
        """
        Request MetabRef GC/MS library.

        Parameters
        ----------
        format : str
            Format of requested library, i.e. "json", "sql", "flashentropy".
            See `available_formats` method for aliases.
        normalize : bool
            Normalize the spectrum by its magnitude.

        Returns
        -------
        Library in requested format.

        """

        # Init format function
        format_func = self._get_format_func(format)

        return format_func(
            self.get_query(self.GCMS_LIBRARY_URL)["GC-MS"], normalize, {}
        )

    def get_fames(self, format="json", normalize=False):
        """
        Request MetabRef GC/MS FAMEs library.

        Parameters
        ----------
        format : str
            Format of requested library, i.e. "json", "sql", "flashentropy".
            See `available_formats` method for aliases.
        normalize : bool
            Normalize the spectrum by its magnitude.

        Returns
        -------
        Library in requested format.

        """

        # Init format function
        format_func = self._get_format_func(format)

        return format_func(self.get_query(self.FAMES_URL)["GC-MS"], normalize, {})

    def _to_LowResolutionEICompound_dict(self, metabref_lib, normalize=False):
        """
        Convert MetabRef-formatted library to CoreMS LowResolutionEICompound-formatted
        dictionary for local ingestion.

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
        metadatar_cols = {
            "casno": "cas",
            "inchikey": "inchikey",
            "inchi": "inchi",
            "chebi": "chebi",
            "smiles": "smiles",
            "kegg": "kegg",
            "iupac_name": "iupac_name",
            "traditional_name": "traditional_name",  # Not present in metabref
            "common_name": "common_name",  # Not present in metabref
        }

        # Dictionary to map metabref keys to corems keys
        lowres_ei_compound_cols = {
            "id": "metabref_id",
            "molecule_name": "name",  # Is this correct?
            "classify": "classify",  # Not present in metabref
            "formula": "formula",
            "ri": "ri",
            "rt": "retention_time",
            "source": "source",  # Not present in metabref
            "casno": "casno",
            "comments": "comment",
            "source_temp_c": "source_temp_c",  # Not present in metabref
            "ev": "ev",  # Not present in metabref
            "peak_count": "peaks_count",
            "mz": "mz",
            "abundance": "abundance",
        }

        # Local result container
        corems_lib = []

        # Enumerate spectra
        for i, source_ in enumerate(metabref_lib):
            # Copy source to prevent modification
            source = source_.copy()

            # Flatten source dict
            source = source.pop("spectrum_data") | source

            # Parse target data
            target = {
                lowres_ei_compound_cols[k]: v
                for k, v in source.items()
                if k in lowres_ei_compound_cols
            }

            # Explicitly add this to connect with LowResCompoundRef later
            target["rt"] = source["rt"]

            # Parse (mz, abundance)
            arr = self.spectrum_to_array(target["mz"], normalize=normalize)
            target["mz"] = arr[:, 0]
            target["abundance"] = arr[:, 1]

            # Parse meta data
            target["metadata"] = {
                metadatar_cols[k]: v for k, v in source.items() if k in metadatar_cols
            }

            # Add anything else
            for k in source:
                if k not in lowres_ei_compound_cols:
                    target[k] = source[k]

            # Add to CoreMS list
            corems_lib.append(target)

        return corems_lib

    def _LowResolutionEICompound_dict_to_sqlite(
        self, lowres_ei_compound_dict, url="sqlite://"
    ):
        """
        Convert CoreMS LowResolutionEICompound-formatted dictionary to SQLite
        database for local ingestion.

        Parameters
        ----------
        lowres_ei_compound_dict : dict
            CoreMS GC-MS library formatted for LowResolutionEICompound.
        url : str
            URL to SQLite prefix.

        Returns
        -------
        sqlite database
            Spectra contained in SQLite database.

        """

        # Dictionary to map corems keys to all-caps keys
        capped_cols = {
            "name": "NAME",
            "formula": "FORM",
            "ri": "RI",
            "retention_time": "RT",
            "source": "SOURCE",
            "casno": "CASNO",
            "comment": "COMMENT",
            "peaks_count": "NUM PEAKS",
        }

        # Initialize SQLite object
        sqlite_obj = EI_LowRes_SQLite(url=url)

        # Iterate spectra
        for _data_dict in lowres_ei_compound_dict:
            # Copy source to prevent modification
            data_dict = _data_dict.copy()

            # Add missing capped values
            for k, v in capped_cols.items():
                # Key exists
                if k in data_dict:
                    # # This will replace the key
                    # data_dict[v] = data_dict.pop(k)

                    # This will keep both keys
                    data_dict[v] = data_dict[k]

            # Parse number of peaks
            if not data_dict.get("NUM PEAKS"):
                data_dict["NUM PEAKS"] = len(data_dict.get("mz"))

            # Parse CAS number
            if not data_dict.get("CASNO"):
                data_dict["CASNO"] = data_dict.get("CAS")

            if not data_dict["CASNO"]:
                data_dict["CASNO"] = 0

            # Build linked metadata table
            if "metadata" in data_dict:
                if len(data_dict["metadata"]) > 0:
                    data_dict["metadatar"] = Metadatar(**data_dict.pop("metadata"))
                else:
                    data_dict.pop("metadata")

            # Attempt addition to sqlite
            try:
                sqlite_obj.add_compound(data_dict)
            except:
                print(data_dict["NAME"])

        return sqlite_obj


class MetabRefLCInterface(MetabRefInterface):
    """
    Interface to the Metabolomics Reference Database for LC-MS data.
    """

    def __init__(self):
        """
        Initialize instance.

        """

        super().__init__()

        # API endpoint for precursor m/z search
        # inputs = mz, tolerance (in Da), polarity, page_no, per_page
        self.PRECURSOR_MZ_URL = "https://metabref.emsl.pnnl.gov/api/precursors/m/{}/t/{}/{}?page={}&per_page={}"

        # API endpoint for returning full list of precursor m/z values in database
        # inputs = polarity, page_no, per_page
        self.PRECURSOR_MZ_ALL_URL = (
            "https://metabref.emsl.pnnl.gov/api/precursors/{}?page={}&per_page={}"
        )

        self.__init_format_map__()

    def __init_format_map__(self):
        """
        Initialize database format mapper, enabling multiple format requests.

        """

        # Define format workflows
        self.format_map = {
            "json": lambda x, normalize, fe_kwargs: x,
            "flashentropy": lambda x, normalize, fe_kwargs: self._to_flashentropy(
                x, normalize, fe_kwargs
            ),
        }

        # Add aliases
        self.format_map["metabref"] = self.format_map["json"]
        self.format_map["fe"] = self.format_map["flashentropy"]
        self.format_map["flash-entropy"] = self.format_map["flashentropy"]

    def query_by_precursor(
        self, mz_list, polarity, mz_tol_ppm, mz_tol_da_api=0.2, max_per_page=50
    ):
        """
        Query MetabRef by precursor m/z values.

        Parameters
        ----------
        mz_list : list
            List of precursor m/z values.
        polarity : str
            Ionization polarity, either "positive" or "negative".
        mz_tol_ppm : float
            Tolerance in ppm for each precursor m/z value.
            Used for retrieving from a potential match from database.
        mz_tol_da_api : float, optional
            Maximum tolerance between precursor m/z values for API search, in daltons.
            Used to group similar mzs into a single API query for speed. Default is 0.2.
        max_per_page : int, optional
            Maximum records to return from MetabRef API query at a time.  Default is 50.

        Returns
        -------
        list
            List of library entries in original JSON format.
        """

        # If polarity is anything other than positive or negative, raise error
        if polarity not in ["positive", "negative"]:
            raise ValueError("Polarity must be 'positive' or 'negative'")

        # Cluster groups of mz according to mz_tol_da_api for precursor query
        mz_list.sort()
        mz_groups = [[mz_list[0]]]
        for x in mz_list[1:]:
            if abs(x - mz_groups[-1][0]) <= mz_tol_da_api:
                mz_groups[-1].append(x)
            else:
                mz_groups.append([x])

        # Query MetabRef for each mz group
        lib = []
        for mz_group in mz_groups:
            mz = np.mean(mz_group)
            if len(mz_group) == 1:
                mz = mz_group[0]
                tol = mz_tol_ppm * 10**-6 * mz
            else:
                mz = (max(mz_group) - min(mz_group)) / 2 + min(mz_group)
                tol = (max(mz_group) - min(mz_group)) / 2 + mz_tol_ppm**-6 * max(
                    mz_group
                )

            # Get first page of results
            response = self.get_query(
                self.PRECURSOR_MZ_URL.format(
                    str(mz), str(tol), polarity, 1, max_per_page
                )
            )
            lib = lib + response["results"]

            # If there are more pages of results, get them
            if response["total_pages"] > 1:
                for i in np.arange(2, response["total_pages"] + 1):
                    lib = (
                        lib
                        + self.get_query(
                            self.PRECURSOR_MZ_URL.format(
                                str(mz), str(tol), polarity, i, max_per_page
                            )
                        )["results"]
                    )

        return lib

    def request_all_precursors(self, polarity, per_page=50000):
        """
        Request all precursor m/z values for MS2 spectra from MetabRef.

        Parameters
        ----------
        polarity : str
            Ionization polarity, either "positive" or "negative".
        per_page : int, optional
            Number of records to fetch per call. Default is 50000

        Returns
        -------
        list
            List of all precursor m/z values, sorted.
        """
        # If polarity is anything other than positive or negative, raise error
        if polarity not in ["positive", "negative"]:
            raise ValueError("Polarity must be 'positive' or 'negative'")

        precursors = []

        # Get first page of results and total number of pages of results
        response = self.get_query(
            self.PRECURSOR_MZ_ALL_URL.format(polarity, str(1), str(per_page))
        )
        total_pages = response["total_pages"]
        precursors.extend([x["precursor_ion"] for x in response["results"]])

        # Go through remaining pages of results
        for i in np.arange(2, total_pages + 1):
            response = self.get_query(
                self.PRECURSOR_MZ_ALL_URL.format(polarity, str(i), str(per_page))
            )
            precursors.extend([x["precursor_ion"] for x in response["results"]])

        # Sort precursors from smallest to largest and remove duplicates
        precursors = list(set(precursors))
        precursors.sort()

        return precursors

    def get_lipid_library(
        self,
        mz_list,
        polarity,
        mz_tol_ppm,
        mz_tol_da_api=0.2,
        format="json",
        normalize=True,
        fe_kwargs={},
    ):
        """
        Request MetabRef lipid library.

        Parameters
        ----------
        mz_list : list
            List of precursor m/z values.
        polarity : str
            Ionization polarity, either "positive" or "negative".
        mz_tol_ppm : float
            Tolerance in ppm for each precursor m/z value.
            Used for retrieving from a potential match from database.
        mz_tol_da_api : float, optional
            Maximum tolerance between precursor m/z values for API search, in daltons.
            Used to group similar mzs into a single API query for speed. Default is 0.2.
        format : str, optional
            Format of requested library, i.e. "json", "sql", "flashentropy".
            See `available_formats` method for aliases. Default is "json".
        normalize : bool, optional
            Normalize the spectrum by its magnitude. Default is True.
        fe_kwargs : dict, optional
            Keyword arguments for FlashEntropy search. Default is {}.

        Returns
        -------
        tuple
            Library in requested format and lipid metadata as a LipidMetadata dataclass.

        """
        mz_list.sort()
        mz_list = np.array(mz_list)

        # Get all precursors in the library matching the polarity
        precusors_in_lib = self.request_all_precursors(polarity=polarity)
        precusors_in_lib = np.array(precusors_in_lib)

        # Compare the mz_list with the precursors in the library, keep any mzs that are within mz_tol of any precursor in the library
        lib_mz_df = pd.DataFrame(precusors_in_lib, columns=["lib_mz"])
        lib_mz_df["closest_obs_mz"] = mz_list[
            find_closest(mz_list, lib_mz_df.lib_mz.values)
        ]
        lib_mz_df["mz_diff_ppm"] = np.abs(
            (lib_mz_df["lib_mz"] - lib_mz_df["closest_obs_mz"])
            / lib_mz_df["lib_mz"]
            * 1e6
        )
        lib_mz_sub = lib_mz_df[lib_mz_df["mz_diff_ppm"] <= mz_tol_ppm]

        # Do the same in the opposite direction
        mz_df = pd.DataFrame(mz_list, columns=["mass_feature_mz"])
        mz_df["closest_lib_pre_mz"] = precusors_in_lib[
            find_closest(precusors_in_lib, mz_df.mass_feature_mz.values)
        ]
        mz_df["mz_diff_ppm"] = np.abs(
            (mz_df["mass_feature_mz"] - mz_df["closest_lib_pre_mz"])
            / mz_df["mass_feature_mz"]
            * 1e6
        )
        mz_df_sub = mz_df[mz_df["mz_diff_ppm"] <= mz_tol_ppm]

        # Evaluate which is fewer mzs - lib_mz_sub or mz_df_sub and use that as the input for next step
        if len(lib_mz_sub) < len(mz_df_sub):
            mzs_to_query = lib_mz_sub.lib_mz.values
        else:
            mzs_to_query = mz_df_sub.mass_feature_mz.values

        # Query the library for the precursors in the mz_list that are in the library to retrieve the spectra and metadata
        lib = self.query_by_precursor(
            mz_list=mzs_to_query,
            polarity=polarity,
            mz_tol_ppm=mz_tol_ppm,
            mz_tol_da_api=mz_tol_da_api,
        )

        # Pull out lipid metadata from the metabref library and convert to LipidMetadata dataclass
        mol_data_dict = {x["id"]: x["Molecular Data"] for x in lib}
        lipid_lib = {x["id"]: x["Lipid Tree"] for x in lib if "Lipid Tree" in x.keys()}
        mol_data_dict = {k: {**v, **lipid_lib[k]} for k, v in mol_data_dict.items()}
        mol_data_dict = {
            k: self._dict_to_dataclass(v, LipidMetadata)
            for k, v in mol_data_dict.items()
        }

        # Remove lipid metadata from the metabref library
        lib = [
            {k: v for k, v in x.items() if k not in ["Molecular Data", "Lipid Tree"]}
            for x in lib
        ]
        # Unpack the 'Lipid Fragments' key and the 'MSO Data" key from each entry
        for x in lib:
            if "Lipid Fragments" in x.keys():
                x.update(x.pop("Lipid Fragments"))
            if "MSO Data" in x.keys():
                x.update(x.pop("MSO Data"))

        # Format the spectral library
        format_func = self._get_format_func(format)
        lib = format_func(lib, normalize=normalize, fe_kwargs=fe_kwargs)
        return (lib, mol_data_dict)


class MSPInterface(SpectralDatabaseInterface):
    """
    Interface to parse NIST MSP files
    """

    def __init__(self, file_path):
        """
        Initialize instance.

        Parameters
        ----------
        file_path : str
            Path to a local MSP file.

        Attributes
        ----------
        file_path : str
            Path to the MSP file.
        _file_content : str
            Content of the MSP file.
        _data_frame : :obj:`~pandas.DataFrame`
            DataFrame of spectra from the MSP file with unaltered content.
        """
        super().__init__(key=None)

        self.file_path = file_path
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(
                f"File {self.file_path} does not exist. Please check the file path."
            )
        with open(self.file_path, "r") as f:
            self._file_content = f.read()

        self._data_frame = self._read_msp_file()
        self.__init_format_map__()

    def __init_format_map__(self):
        """
        Initialize database format mapper, enabling multiple format requests.

        """

        # x is a pandas dataframe similar to self._data_frame format
        # Define format workflows
        self.format_map = {
            "msp": lambda x, normalize, fe_kwargs: self._to_msp(x, normalize),
            "flashentropy": lambda x, normalize, fe_kwargs: self._to_flashentropy(
                x, normalize, fe_kwargs
            ),
            "df": lambda x, normalize, fe_kwargs: self._to_df(x, normalize),
        }

        # Add aliases
        self.format_map["fe"] = self.format_map["flashentropy"]
        self.format_map["flash-entropy"] = self.format_map["flashentropy"]
        self.format_map["dataframe"] = self.format_map["df"]
        self.format_map["data-frame"] = self.format_map["df"]

    def _read_msp_file(self):
        """
        Reads the MSP files into the pandas dataframe, and sort/remove zero intensity ions in MS/MS spectra.

        Returns
        -------
        :obj:`~pandas.DataFrame`
            DataFrame of spectra from the MSP file, exacly as it is in the file (no sorting, filtering etc)
        """
        # If input_dataframe is provided, return it it
        spectra = []
        spectrum = {}

        f = StringIO(self._file_content)
        for line in f:
            line = line.strip()
            if not line:
                continue  # Skip empty lines

            # Handle metadata
            if ":" in line:
                key, value = line.split(":", 1)
                key = key.strip().lower()
                value = value.strip()

                if key == "name":
                    # Save current spectrum and start a new one
                    if spectrum:
                        spectra.append(spectrum)
                    spectrum = {"name": value, "peaks": []}
                else:
                    spectrum[key] = value

            # Handle peak data (assumed to start with a number)
            elif line[0].isdigit():
                peaks = line.split()
                m_z = float(peaks[0])
                intensity = float(peaks[1])
                spectrum["peaks"].append(([m_z, intensity]))
        # Save the last spectrum
        if spectrum:
            spectra.append(spectrum)

        df = pd.DataFrame(spectra)
        for column in df.columns:
            if column != "peaks":  # Skip 'peaks' column
                try:
                    df[column] = pd.to_numeric(df[column], errors="raise")
                except:
                    pass
        return df

    def _to_df(self, input_dataframe, normalize=True):
        """
        Convert MSP-derived library to FlashEntropy library. 

        Parameters
        ----------
        input_dataframe : :obj:`~pandas.DataFrame`
            Input DataFrame containing MSP-formatted spectra.
        normalize : bool, optional
            Normalize each spectrum by its magnitude.
            Default is True.

        Returns
        -------
        :obj:`~pandas.DataFrame`
            DataFrame of with desired normalization
        """
        if not normalize:
            return input_dataframe
        else:
            # Convert to dictionary
            db_dict = input_dataframe.to_dict(orient="records")

            # Initialize empty library
            lib = []

            # Enumerate spectra
            for i, source in enumerate(db_dict):
                spectrum = source
                # Check that spectrum["peaks"] exists
                if "peaks" not in spectrum.keys():
                    raise KeyError(
                        "MSP not interpretted correctly, 'peaks' key not found in spectrum, check _dataframe attribute."
                    )

                # Convert spectrum["peaks"] to numpy array
                if not isinstance(spectrum["peaks"], np.ndarray):
                    spectrum["peaks"] = np.array(spectrum["peaks"])

                # Normalize peaks, if requested
                if normalize:
                    spectrum["peaks"] = self.normalize_peaks(spectrum["peaks"])
                    spectrum["num peaks"] = len(spectrum["peaks"])

                # Add spectrum to library
                lib.append(spectrum)
            
            # Convert to DataFrame
            df = pd.DataFrame(lib)
            return df
    
    def _to_flashentropy(self, input_dataframe, normalize=True, fe_kwargs={}):
        """
        Convert MSP-derived library to FlashEntropy library.

        Parameters
        ----------
        input_dataframe : :obj:`~pandas.DataFrame`
            Input DataFrame containing MSP-formatted spectra.
        normalize : bool
            Normalize each spectrum by its magnitude.
        fe_kwargs : dict, optional
            Keyword arguments for instantiation of FlashEntropy search and building index for FlashEntropy search;
            any keys not recognized will be ignored. By default, all parameters set to defaults.

        Returns
        -------
        :obj:`~ms_entropy.FlashEntropySearch`
            MS2 library as FlashEntropy search instance.

        Raises
        ------
        ValueError
            If "min_ms2_difference_in_da" or "max_ms2_tolerance_in_da" are present in `fe_kwargs` and they
        """
        self._check_flash_entropy_kwargs(fe_kwargs)

        db_df = input_dataframe

        # Convert to dictionary
        db_dict = db_df.to_dict(orient="records")

        # Initialize empty library
        fe_lib = []

        # Enumerate spectra
        for i, source in enumerate(db_dict):
            # Reorganize source dict, if necessary
            if "spectrum_data" in source.keys():
                spectrum = source["spectrum_data"]
            else:
                spectrum = source

            # Rename precursor_mz key for FlashEntropy
            if "precursor_mz" not in spectrum.keys():
                if "precursormz" in spectrum:
                    spectrum["precursor_mz"] = spectrum.pop("precursormz")
                elif "precursor_ion" in spectrum:
                    spectrum["precursor_mz"] = spectrum.pop("precursor_ion")
                else:
                    raise KeyError(
                        "MSP must have either 'precursormz' or 'precursor_ion' key to be converted to FlashEntropy format."
                    )

            # Check that spectrum["peaks"] exists
            if "peaks" not in spectrum.keys():
                raise KeyError(
                    "MSP not interpretted correctly, 'peaks' key not found in spectrum, check _dataframe attribute."
                )

            # Convert spectrum["peaks"] to numpy array
            if not isinstance(spectrum["peaks"], np.ndarray):
                spectrum["peaks"] = np.array(spectrum["peaks"])

            # Normalize peaks, if requested
            if normalize:
                spectrum["peaks"] = self.normalize_peaks(spectrum["peaks"])

            # Add spectrum to library
            fe_lib.append(spectrum)

        # Build FlashEntropy index
        fe_search = self._build_flash_entropy_index(fe_lib, fe_kwargs=fe_kwargs)

        return fe_search

    def _to_msp(self, input_dataframe, normalize=True):
        #TODO KRH: Write this functionality or remove before merging into master branch

        raise NotImplementedError(
            "MSP writing functionality not yet implemented."
        )
    
    def _check_msp_compatibility(self):
        """
        Check if the MSP file is compatible with the get_metabolomics_spectra_library method and provide feedback if it is not.
        """
        # Check polarity
        if (
            "polarity" not in self._data_frame.columns
            and "ionmode" not in self._data_frame.columns
        ):
            raise ValueError(
                "Neither 'polarity' nor 'ionmode' columns found in the input MSP metadata. Please check the file."
            )
        polarity_column = (
            "polarity" if "polarity" in self._data_frame.columns else "ionmode"
        )

        # Check if polarity_column contents is either "positive" or "negative"
        if not all(self._data_frame[polarity_column].isin(["positive", "negative"])):
            raise ValueError(
                f"Input field on MSP '{polarity_column}' must contain only 'positive' or 'negative' values."
            )

        # Check if the MSP file contains the required columns for metabolite metadata
        # inchikey, by name, not null
        # either formula or molecular_formula, not null
        if not all(self._data_frame["inchikey"].notnull()):
            raise ValueError(
                "Input field on MSP 'inchikey' must contain only non-null values."
            )
        if (
            "formula" not in self._data_frame.columns
            and "molecular_formula" not in self._data_frame.columns
        ):
            raise ValueError(
                "Input field on MSP must contain either 'formula' or 'molecular_formula' columns."
            )
        molecular_formula_column = (
            "formula" if "formula" in self._data_frame.columns else "molecular_formula"
        )
        if not all(self._data_frame[molecular_formula_column].notnull()):
            raise ValueError(
                f"Input field on MSP '{molecular_formula_column}' must contain only non-null values."
            )

    def get_metabolomics_spectra_library(
        self,
        polarity,
        metabolite_metadata_mapping={},
        format="fe",
        normalize=True,
        fe_kwargs={},
    ):
        """
        Prepare metabolomics spectra library and associated metabolite metadata

        Note: this uses the inchikey as the index for the metabolite metadata dataframe and for connecting to the spectra, so it must be in the input

        """
        # Check if the MSP file is compatible with the get_metabolomics_spectra_library method
        self._check_msp_compatibility()

        # Check if the polarity parameter is valid and if a polarity column exists in the dataframe
        if polarity not in ["positive", "negative"]:
            raise ValueError("Polarity must be 'positive' or 'negative'")
        polarity_column = (
            "polarity" if "polarity" in self._data_frame.columns else "ionmode"
        )

        # Get a subset of the initial dataframea by polarity
        db_df = self._data_frame[self._data_frame[polarity_column] == polarity].copy()

        # Rename the columns of the db_df to match the MetaboliteMetadata dataclass using the metabolite_metadata_mapping
        # If the mapping is not provided, use the default mapping
        if not metabolite_metadata_mapping:
            metabolite_metadata_mapping = {
                "chebi_id": "chebi",
                "kegg_id": "kegg",
                "refmet_name": "common_name",
                "molecular_formula": "formula",
                "gnps_spectra_id":"id",
                "precursormz": "precursor_mz",
                "precursortype":"ion_type"
            }
        db_df.rename(columns=metabolite_metadata_mapping, inplace=True)
        db_df["molecular_data_id"] = db_df["inchikey"]



        # Check if the resulting dataframe has the required columns for the flash entropy search
        required_columns = ["molecular_data_id", "precursor_mz", "ion_type", "id"]
        for col in required_columns:
            if col not in db_df.columns:
                raise ValueError(
                    f"Input field on MSP must contain '{col}' column for FlashEntropy search."
                )

        # Pull out the metabolite metadata from the dataframe and put it into a different dataframe
        # First get a list of the possible attributes of the MetaboliteMetadata dataclass
        metabolite_metadata_keys = list(MetaboliteMetadata.__annotations__.keys())
        # Replace id with molecular_data_id in metabolite_metadata_keys
        metabolite_metadata_keys = [
            "molecular_data_id" if x == "id" else x for x in metabolite_metadata_keys
        ]
        metabolite_metadata_df = db_df[
            db_df.columns[db_df.columns.isin(metabolite_metadata_keys)]
        ].copy()

        # Make unique and recast the id column for metabolite metadata
        metabolite_metadata_df.drop_duplicates(subset=["molecular_data_id"], inplace=True)
        metabolite_metadata_df["id"] = metabolite_metadata_df["molecular_data_id"]

        # Convert to a dictionary using the inchikey as the key
        metabolite_metadata_dict = metabolite_metadata_df.to_dict(
            orient="records"
        )
        metabolite_metadata_dict = {
            v["id"]: self._dict_to_dataclass(v, MetaboliteMetadata)
            for v in metabolite_metadata_dict
        }

        # Remove the metabolite metadata columns from the original dataframe
        for key in metabolite_metadata_keys:
            if key != "molecular_data_id":
                if key in db_df.columns:
                    db_df.drop(columns=key, inplace=True)

        # Format the spectral library
        format_func = self._get_format_func(format)
        lib = format_func(db_df, normalize=normalize, fe_kwargs=fe_kwargs)
        return (lib, metabolite_metadata_dict)
