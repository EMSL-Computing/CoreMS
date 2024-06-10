class SpectrumSearchResults:
    """Class for storing Search Results for a single Spectrum Query

    Parameters
    ----------
    query_spectrum : MassSpectrum
        The queried mass spectrum
    precursor_mz : float, optional
        The queried precursor_mz. None is interpreted as an open query, i.e. no precursor_mz
    spectral_similarity_search_results : dict
        The search results for the queried spectrum, which will be unpacked into attributes

    Attributes
    ----------
    query_spectrum : MassSpectrum
        The queried mass spectrum
    query_spectrum_id : int
        The id of the queried spectrum (the scan number within an MassSpectra object)
    precursor_mz : float
        The precursor m/z of the queried spectrum

    Other Possible Attributes
    -------------------------
    ref_mol_id : str
        The id of the molecule associated with the query spectrum in reference database
    ref_ms_id : str
        The id of the query spectrum in reference database
    ref_precursor_mz : float
        The precursor mass of the query spectrum
    precursor_mz_error_ppm : float
        The ppm error between the query spectrum and the reference spectrum
    entropy_similarity : float
        The entropy similarity between the query spectrum and the reference spectrum
    ref_ion_type : str
        The ion type of the reference spectrum, i.e. [M+H]+, [M+Na]+, etc.
    query_mz_in_ref_n : list
        The number of query m/z peaks that are in the reference spectrum
    query_mz_in_ref_fract : float
        The fraction of query m/z peaks that are in the reference spectrum
    query_frag_types : list
        The fragment types of the query spectrum that are in the reference spectrum,
        i.e. LSF (lipid species fragments) or MSF (molecular species fragments),
        generally used for only for lipidomics
    ref_mz_in_query_n : list
        The number of reference m/z peaks that are in the query spectrum
    ref_mz_in_query_fract : float
        The fraction of reference m/z peaks that are in the query spectrum
    ref_frag_types : list
        The fragment types of the reference spectrum,
        i.e. LSF (lipid species fragments) or MSF (molecular species fragments),
        generally used for only for lipidomics

    Methods
    -------
    *to_dataframe().
        Convert the SpectrumSearchResults to a pandas DataFrame

    """

    def __init__(
        self, query_spectrum, precursor_mz, spectral_similarity_search_results
    ):
        self.query_spectrum = query_spectrum
        self.precursor_mz = precursor_mz
        if query_spectrum is not None:
            if query_spectrum.scan_number is not None:
                self.query_spectrum_id = query_spectrum.scan_number
        attribute_keys = [
            "ref_mol_id",
            "ref_ms_id",
            "ref_precursor_mz",
            "precursor_mz_error_ppm",
            "ref_ion_type",
            "entropy_similarity",
            "query_mz_in_ref_n",
            "query_mz_in_ref_fract",
            "query_frag_types",
            "ref_mz_in_query_n",
            "ref_mz_in_query_fract",
            "ref_frag_types",
        ]
        for key in spectral_similarity_search_results.keys():
            if key in attribute_keys:
                setattr(self, key, spectral_similarity_search_results[key])

    def to_dataframe(self, cols_to_drop=None):
        """Convert the SpectrumSearchResults to a pandas DataFrame

        Parameters
        ----------
        cols_to_drop : list, optional
            A list of columns to drop from the DataFrame. Default is None.

        Returns
        -------
        pandas.DataFrame
            A DataFrame with the SpectrumSearchResults attributes as columns

        """
        import pandas as pd

        # get initial dict
        df = pd.DataFrame.from_dict(self.__dict__).copy()

        # remove ms2_spectrum column
        df = df.drop(columns=["query_spectrum"])

        # drop additional columns
        if cols_to_drop is not None:
            df = df.drop(columns=cols_to_drop)

        # reorder to high to low entropy similarity
        df = df.sort_values(by=["entropy_similarity"], ascending=False)

        # rename id to query_scan_number
        return pd.DataFrame(df)
