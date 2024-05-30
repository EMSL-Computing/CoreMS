__author__ = "Katherine R. Heal"
__date__ = "2024-01-24"

from dataclasses import dataclass

from .EI_SQL import MetaboliteMetadata


@dataclass
class LipidMetadata(MetaboliteMetadata):
    """Dataclass for the Lipid Metadata
    
    Parameters
    ----------
    name : str
        The name of the lipid, using the LIPID MAPS nomenclature
    casno : str
        The CAS number of the lipid
    formula : str
        The molecular formula of the lipid
    pubchem_id : str
        The PubChem ID of the lipid
    structure_level : str
        The structure level of the lipid, following the LIPID MAPS classification
    lipid_summed_name : str
        The summed name of the lipid, aka lipid species, 
        following the LIPID MAPS classification
    lipid_subclass : str
        The subclass of the lipid, following the LIPID MAPS classification
    lipid_class : str
        The class of the lipid, following the LIPID MAPS classification
    lipid_category : str
        The category of the lipid, following the LIPID MAPS classification
    """

    name: str
    casno: str
    formula: str
    pubchem_id: str
    structure_level: str

    lipid_summed_name: str
    lipid_subclass: str
    lipid_class: str
    lipid_category: str


class MS2SearchResults:
    """Class for storing MS2 Search Results for a single MS2 Spectrum

    Parameters
    ----------
    ms2_spectrum : MS2Spectrum
        The queried MS2 Spectrum
    precursor_mz : float, optional
        The queried precursor_mz.  Default is None, which implies an open search.
    ms2_search_results : dict
        The MS2 Search Results for the queried MS2 Spectrum, which will be unpacked into attributes

    Attributes
    ----------
    ms2_spectrum : MS2Spectrum
        The queried MS2 Spectrum
    ms2_spectrum_id : int
        The id of the queried MS2 Spectrum within the LCMSObject (generally the scan number)
    precursor_mz : float
        The precursor m/z of the queried MS2 Spectrum

    Other Attributes Possible
    -------------------------
    metabref_mol_id : str
        The mol_id of the molecule associated with the query spectrum in metabref database
    metabref_ms_id : str
        The ms_id of the query spectrum in metabref database
    metabref_precursor_mz : float
        The precursor mass of the query spectrum
    precursor_mz_error_ppm : float
        The ppm error between the query spectrum and the mass feature
    entropy_similarity : float
        The entropy similarity between the query spectrum and the reference spectrum
    fe_lib_id : str
        The id of the reference spectrum in emphermeral flash entorpy library
    query_mz_in_lib : list
        The number of query m/z peaks that are in the reference spectrum
    query_mz_in_lib_fract : float
        The fraction of query m/z peaks that are in the reference spectrum
    query_frags : list
        The number of query fragments that are in the reference spectrum
    lib_mz_in_query : list
        The number of reference m/z peaks that are in the query spectrum

    Methods
    -------
    *to_dataframe().
        Convert the MS2SeachResults to a pandas DataFrame

    """

    def __init__(self, ms2_spectrum, precursor_mz, ms2_search_results):
        self.ms2_spectrum = ms2_spectrum
        self.precursor_mz = precursor_mz
        if ms2_spectrum is not None:
            if ms2_spectrum.scan_number is not None:
                self.ms2_spectrum_id = ms2_spectrum.scan_number
        attribute_keys = [
            "metabref_mol_id",
            "metabref_ms_id",
            "metabref_precursor_mz",
            "precursor_mz_error_ppm",
            "ion_type",
            "entropy_similarity",
            "query_mz_in_lib",
            "query_mz_in_lib_fract",
            "query_frags",
            "lib_mz_in_query",
            "lib_mz_in_query_fract",
            "lib_frags",
        ]
        for key in ms2_search_results.keys():
            if key in attribute_keys:
                setattr(self, key, ms2_search_results[key])

    def to_dataframe(self):
        """Convert the MS2SearchResults to a pandas DataFrame

        Returns
        -------
        pandas.DataFrame
            A DataFrame with the MS2SearchResults attributes as columns

        """
        import pandas as pd

        # get initial dict
        df = pd.DataFrame.from_dict(self.__dict__).copy()

        # remove ms2_spectrum column
        df = df.drop(columns=["ms2_spectrum", "query_mz_in_lib", "lib_mz_in_query"])

        # reorder to high to low entropy similarity
        df = df.sort_values(by=["entropy_similarity"], ascending=False)

        # rename id to query_scan_number
        df = df.rename(columns={"ms2_spectrum_id": "query_scan_number"})
        return pd.DataFrame(df)