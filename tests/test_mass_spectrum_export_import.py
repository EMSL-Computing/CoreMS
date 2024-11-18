import os
import pytest

from corems.mass_spectrum.input.numpyArray import ms_from_array_centroid
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.mass_spectrum.output.export import HighResMassSpecExport
from corems.mass_spectrum.input.coremsHDF5 import ReadCoreMSHDF_MassSpectrum


@pytest.fixture
def mass_spectrum_silico():
    # Test for generating accurate molecular formula from a single mass using the local sql database
    # Now also tests that it is handling isotopes correctly (for non-adducts)
    mz = [760.58156938877, 761.58548]
    abundance = [1000, 400]
    rp, s2n = [[1, 1], [10, 10]]

    mass_spectrum_obj = ms_from_array_centroid(
        mz, abundance, rp, s2n, "single mf search", polarity=1, auto_process=False
    )

    # Set the settings for the molecular search on the mass spectrum object
    mass_spectrum_obj.settings.noise_threshold_method = "relative_abundance"
    mass_spectrum_obj.settings.noise_threshold_absolute_abundance = 0

    mass_spectrum_obj.molecular_search_settings.url_database = "postgresql://coremsdb:coremsmolform@postgres:5432/molformula"
    mass_spectrum_obj.molecular_search_settings.error_method = "None"
    mass_spectrum_obj.molecular_search_settings.min_ppm_error = -5
    mass_spectrum_obj.molecular_search_settings.max_ppm_error = 5
    mass_spectrum_obj.molecular_search_settings.mz_error_range = 1
    mass_spectrum_obj.molecular_search_settings.isProtonated = True
    mass_spectrum_obj.molecular_search_settings.isRadical = False
    mass_spectrum_obj.molecular_search_settings.isAdduct = False
    mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = False
    mass_spectrum_obj.molecular_search_settings.use_isotopologue_filter = False

    usedatoms = {"C": (1, 57), "H": (4, 200), "N": (0, 1)}
    mass_spectrum_obj.molecular_search_settings.usedAtoms = usedatoms

    mass_spectrum_obj.process_mass_spec()

    return mass_spectrum_obj


def test_molecular_formula_search(mass_spectrum_silico):
    SearchMolecularFormulas(
        mass_spectrum_silico, find_isotopologues=True
    ).run_worker_ms_peaks([mass_spectrum_silico[0]])

    ms_df1 = mass_spectrum_silico.to_dataframe()
    assert mass_spectrum_silico[0][0].string == "C56 H73 N1"
    assert ms_df1.shape == (2, 26)
    assert mass_spectrum_silico[1][0].string == "C55 H73 N1 13C1"


def test_mass_spec_export_import_with_annote(mass_spectrum_silico):
    SearchMolecularFormulas(
        mass_spectrum_silico, find_isotopologues=True
    ).run_worker_ms_peaks([mass_spectrum_silico[0]])

    exportMS = HighResMassSpecExport("my_mass_spec", mass_spectrum_silico)
    exportMS._output_type = "hdf5"
    exportMS.save()

    parser = ReadCoreMSHDF_MassSpectrum("my_mass_spec.hdf5")
    mass_spectrum_obj2 = parser.get_mass_spectrum(auto_process=True, load_settings=True)

    ms_df2 = mass_spectrum_obj2.to_dataframe()
    assert mass_spectrum_obj2[0][0].string == "C56 H73 N1"
    assert ms_df2.shape == (2, 26)
    assert mass_spectrum_obj2[1][0].string == "C55 H73 N1 13C1"
    assert mass_spectrum_obj2._mz_exp[0] == 760.58156938877

    # Remove the file
    os.remove("my_mass_spec.hdf5")