
import os, sys
from pathlib import Path
sys.path.append(".")

from corems.molecular_formula.input.masslist_ref import ImportMassListRef
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.constant import Labels

import pytest

__author__ = "Yuri E. Corilo"
__date__ = "Jul 02, 2019"

def get_mass_spectrum():

    from corems.mass_spectrum.input.massList import ReadMassList

    file_location = Path.cwd() / "tests/tests_data/ftms/" / "ESI_NEG_ESFA.ascii"

    #polarity needs to be set or read from the file

    polarity = -1
   
    return ReadMassList(file_location).get_mass_spectrum(polarity, auto_process=True)
    
def test_search_imported_ref_files():

    mass_spectrum_obj = get_mass_spectrum()
    
    ref_file_location = os.path.join(os.getcwd(),  os.path.normcase("tests/tests_data/ftms/")) + "SRFA.ref"

    mf_references_list = ImportMassListRef(ref_file_location).from_bruker_ref_file()

    for mf in mf_references_list:
        
        print(mf.mz_calc, mf.class_label)
    
    ion_type = 'unknown'

    ms_peaks_assigned = SearchMolecularFormulas(mass_spectrum_obj).search_mol_formulas( mf_references_list, ion_type, neutral_molform=False, find_isotopologues=False)

    assert (len(ms_peaks_assigned)) > 0

if __name__ == '__main__':
    
    test_search_imported_ref_files()