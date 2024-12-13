import sys

from corems.molecular_formula.input.masslist_ref import ImportMassListRef
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas

def test_search_imported_ref_files(mass_spectrum_ftms, ref_file_location):
    mass_spectrum_obj = mass_spectrum_ftms
    mass_spectrum_obj.molecular_search_settings.url_database = "postgresql://coremsdb:coremsmolform@postgres:5432/molformula"
    mf_references_list = ImportMassListRef(ref_file_location).from_bruker_ref_file()
    assert len(mf_references_list) == 60
    assert round(mf_references_list[0].mz_calc, 2) == 149.06
    assert mf_references_list[0].class_label == "O2"

    ion_type = "unknown"

    ms_peaks_assigned = SearchMolecularFormulas(mass_spectrum_obj).search_mol_formulas(
        mf_references_list, ion_type, neutral_molform=False, find_isotopologues=False
    )

    assert (len(ms_peaks_assigned)) > 10
