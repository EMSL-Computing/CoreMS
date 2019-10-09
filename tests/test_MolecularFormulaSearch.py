__author__ = "Yuri E. Corilo"
__date__ = "Jul 25, 2019"

import os
import sys
import time
import pytest
sys.path.append('.')
from enviroms.mass_spectrum.factory.MSPeakClasses import MSPeak
from enviroms.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaSearchSettings, MoleculaLookupDictSettings
from enviroms.mass_spectrum.input.textMassList import Read_MassList
from enviroms.molecular_id.search.MolecularFormulaSearch import SearchMolecularFormulas
from enviroms.molecular_id.calc.ClusterFilter import ClusteringFilter
from enviroms.transient.input.BrukerSolarix import ReadBrukerSolarix

def creat_mass_spectrum():

    directory = os.path.join(os.getcwd(), "tests/tests_data/")
    
    file_name = os.path.normcase("ESI_NEG_SRFA.d/")

    file_location = directory + file_name

    bruker_reader = ReadBrukerSolarix(file_location)

    bruker_transient = bruker_reader.get_transient()
    
    mass_spectrum_obj = bruker_transient.get_mass_spectrum(
        plot_result=False, auto_process=True)


    # polariy need to be set if reading a text file
    #polariy = -1
    # load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    #mass_list_reader = Read_MassList(file_location, polariy, delimiter="  ")
    #mass_spectrum_obj  = mass_list_reader.get_mass_spectrum(auto_process=True)

    return mass_spectrum_obj

def test_mspeak_search():

    mass_spec_obj = creat_mass_spectrum()

    mspeak_obj = mass_spec_obj.most_abundant_mspeak
    
    SearchMolecularFormulas().run_worker_ms_peak(mspeak_obj, mass_spec_obj)

    if mspeak_obj.is_assigned:

        print(mspeak_obj[0].mz_error, mspeak_obj[0].to_string_formated)

def test_molecular_formula_search_db():
    
    mass_spec_obj = creat_mass_spectrum()
    
    time1 = time.time()
    
    SearchMolecularFormulas(first_hit=False).run_worker_mass_spectrum(mass_spec_obj)
    
    print('searching molecular formulas took %.3f seconds' % (time.time() - time1))
    
    i = 0
    j = 0
    error = list()
    mass = list()
    abundance = list()
    
    for mspeak in mass_spec_obj.sort_by_abundance():
        
        if mspeak.is_assigned:
            i += 1
            for mformula in mspeak:
                mass.append(mspeak.mz_exp)
                error.append(mformula._calc_assigment_mass_error(mspeak.mz_exp))
                abundance.append(mspeak.abundance)
        else:
            j += 1
            pass
    
    print('%i peaks assigned and %i peaks not assigned' % (i, j))
    
if __name__ == "__main__":

    #test_molecular_formula_search_db()
    test_mspeak_search()