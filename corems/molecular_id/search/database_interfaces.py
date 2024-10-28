import os
import re
from abc import ABC

import numpy as np
import requests
import pandas as pd
from ms_entropy import FlashEntropySearch

from corems.molecular_id.factory.EI_SQL import EI_LowRes_SQLite, Metadatar
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

        if self.key is None:
            raise ValueError(
                "Must specify environment variable key for token associatedwith this database interface."
            )

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

    def get_query(self, url):
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
        response = requests.get(url, headers=self.get_header())

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


class MetabRefInterface(SpectralDatabaseInterface):
    """
    Interface to the Metabolomics Reference Database.
    """

    def __init__(self):
        """
        Initialize instance.

        """

        super().__init__(key="METABREF_TOKEN")

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

        # Normalize the array
        if normalize:
            arr[:, -1] = arr[:, -1] / arr[:, -1].sum()

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
        fes.build_index(fe_lib, **fe_index_kws, clean_spectra=True)

        return fes

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
        self.PRECURSOR_MZ_URL = (
            "https://metabref.emsl.pnnl.gov/api/precursors/m/{}/t/{}/{}"
        )

        # API endpoint for returning full list of precursor m/z values in database
        self.PRECURSOR_MZ_ALL_URL = "https://metabref.emsl.pnnl.gov/api/precursors/{}"

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

    def query_by_precursor(self, mz_list, polarity, mz_tol_ppm, mz_tol_da_api=0.2):
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
            if abs(x - mz_groups[-1][-1]) <= mz_tol_da_api:
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
            lib = lib + self.get_query(
                self.PRECURSOR_MZ_URL.format(str(mz), str(tol), polarity)
            )

        return lib

    def request_all_precursors(self, polarity):
        """
        Request all precursor m/z values from MetabRef.

        Parameters
        ----------
        polarity : str
            Ionization polarity, either "positive" or "negative".

        Returns
        -------
        list
            List of all precursor m/z values.
        """
        # If polarity is anything other than positive or negative, raise error
        if polarity not in ["positive", "negative"]:
            raise ValueError("Polarity must be 'positive' or 'negative'")

        # Query MetabRef for all precursor m/z values
        return self.get_query(self.PRECURSOR_MZ_ALL_URL.format(polarity))

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

        # Get all precursors in the library matching the polarity
        precusors_in_lib = self.request_all_precursors(polarity=polarity)
        precusors_in_lib.sort()
        precusors_in_lib = np.array(precusors_in_lib)

        # Compare the mz_list with the precursors in the library, keep any mzs that are within mz_tol of any precursor in the library
        mz_list = np.array(mz_list)
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

        # Query the library for the precursors in the mz_list that are in the library to retrieve the spectra and metadata
        lib = self.query_by_precursor(
            mz_list=mz_df_sub.mass_feature_mz.values,
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

        # Format the spectral library
        format_func = self._get_format_func(format)
        lib = format_func(lib, normalize=normalize, fe_kwargs=fe_kwargs)
        return (lib, mol_data_dict)
