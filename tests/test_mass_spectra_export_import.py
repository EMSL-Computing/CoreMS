import shutil

from corems.mass_spectrum.input.numpyArray import ms_from_array_centroid
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters
from corems.mass_spectrum.output.export import HighResMassSpecExport
from corems.mass_spectrum.input.coremsHDF5 import ReadCoreMSHDF_MassSpectrum


def prep_mass_spec_obj():
    # Test for generating accurate molecular formula from a single mass using the local sql database
    # Now also tests that it is handling isotopes correctly (for non-adducts)
    mz = [760.58156938877, 761.58548]
    abundance = [1000, 400]
    rp, s2n = [[1, 1], [10, 10]]

    MSParameters.mass_spectrum.noise_threshold_method = "relative_abundance"
    MSParameters.mass_spectrum.noise_threshold_absolute_abundance = 0

    MSParameters.molecular_search.url_database = ""
    MSParameters.molecular_search.error_method = "None"
    MSParameters.molecular_search.min_ppm_error = -5
    MSParameters.molecular_search.max_ppm_error = 5
    MSParameters.molecular_search.mz_error_range = 1
    MSParameters.molecular_search.isProtonated = True
    MSParameters.molecular_search.isRadical = False
    MSParameters.molecular_search.isAdduct = False

    usedatoms = {"C": (1, 57), "H": (4, 200), "N": (0, 1)}
    MSParameters.molecular_search.usedAtoms = usedatoms
    mass_spectrum_obj = ms_from_array_centroid(
        mz, abundance, rp, s2n, "single mf search", polarity=1, auto_process=True
    )
    return mass_spectrum_obj


def run_molecular_formula_search(mass_spectrum_obj):
    mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = False
    mass_spectrum_obj.molecular_search_settings.use_isotopologue_filter = False
    SearchMolecularFormulas(
        mass_spectrum_obj, find_isotopologues=True
    ).run_worker_ms_peaks([mass_spectrum_obj[0]])
    return mass_spectrum_obj


def test_mass_spec_export_import_with_annote():
    mass_spectrum_obj = prep_mass_spec_obj()
    mass_spectrum_obj = run_molecular_formula_search(mass_spectrum_obj)
    ms_df1 = mass_spectrum_obj.to_dataframe()
    assert mass_spectrum_obj[0][0].string == "C56 H73 N1"
    assert ms_df1.shape == (2, 26)
    assert mass_spectrum_obj[1][0].string == "C55 H73 N1 13C1"

    exportMS = HighResMassSpecExport("my_mass_spec", mass_spectrum_obj)
    exportMS._output_type = "hdf5"
    exportMS.save()

    parser = ReadCoreMSHDF_MassSpectrum("my_mass_spec.hdf5")
    mass_spectrum_obj2 = parser.get_mass_spectrum(auto_process=True, load_settings=True)

    ms_df2 = mass_spectrum_obj2.to_dataframe()
    assert mass_spectrum_obj2[0][0].string == "C56 H73 N1"
    assert ms_df2.shape == (2, 26)
    assert mass_spectrum_obj2[1][0].string == "C55 H73 N1 13C1"

    # Remove the file
    shutil.rmtree(
        "my_mass_spec.hdf5",
        ignore_errors=True,
    )


if __name__ == "__main__":
    test_mass_spec_export_import_with_annote()
