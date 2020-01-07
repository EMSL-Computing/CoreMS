__author__ = "Yuri E. Corilo"
__date__ = "Jul 25, 2019"


import sys
import time
from pathlib import Path
sys.path.append('.')

import pytest

from corems.mass_spectrum.input.numpyArray import ms_from_array_centroid
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.transient.input.BrukerSolarix import ReadBrukerSolarix
from corems.encapsulation.settings.molecular_id.MolecularIDSettings import MolecularSearchSettings

def create_mass_spectrum():
    
    file_location = Path.cwd() /  "ESI_NEG_SRFA.d"

    bruker_reader = ReadBrukerSolarix(file_location)

    bruker_transient = bruker_reader.get_transient()
    
    mass_spectrum_obj = bruker_transient.get_mass_spectrum(
        plot_result=False, auto_process=True)


    # polariy need to be set if reading a text file
    #polariy = -1
    # load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    #mass_list_reader = Read_MassList(file_location,  )
    #mass_spectrum_obj  = mass_list_reader.get_mass_spectrum(polarity, auto_process=True)

    return mass_spectrum_obj


def test_run_molecular_formula_search():

    MolecularSearchSettings.usedAtoms['F'] = (0,2)
    MolecularSearchSettings.usedAtoms['P'] = (0,2)
    MolecularSearchSettings.usedAtoms['Cl'] = (0,1)
    MolecularSearchSettings.isAdduct = True

    MolecularSearchSettings.used_atom_valences['P'] = 3
    MolecularSearchSettings.used_atom_valences['F'] = 1
    MolecularSearchSettings.used_atom_valences['Cl'] = 0

    mz = [215.09269]
    abundance = [1]
    rp, s2n = [1,1]
    dataname = 'one peak'
    mass_spectrum_obj = ms_from_array_centroid(mz, abundance, rp, s2n, dataname)

    SearchMolecularFormulas().run_worker_ms_peak(mass_spectrum_obj[0], mass_spectrum_obj)
    ms_peak = mass_spectrum_obj[0]
    print(ms_peak.mz_exp)
    if ms_peak.is_assigned:
        for formula in ms_peak:
            print(formula.to_string_formated, formula.mz_error)

    MolecularSearchSettings.usedAtoms['F'] = (0,0)
    MolecularSearchSettings.usedAtoms['P'] = (0,0)
    MolecularSearchSettings.usedAtoms['Cl'] = (0,0)
    MolecularSearchSettings.isAdduct = False    

def test_mspeak_search():

    
    mass_spec_obj = create_mass_spectrum()
    
    print("OK")

    mspeak_obj = mass_spec_obj.most_abundant_mspeak
    
    SearchMolecularFormulas().run_worker_ms_peak(mspeak_obj, mass_spec_obj)

    print("OK2")
    if mspeak_obj.is_assigned:
        
        print(mspeak_obj.molecular_formula_earth_filter().to_string)
        print(mspeak_obj.molecular_formula_water_filter().to_string)
        print(mspeak_obj.molecular_formula_air_filter().to_string)
        print(mspeak_obj.cia_score_S_P_error().to_string)
        print(mspeak_obj.cia_score_N_S_P_error().to_string)

        print(mspeak_obj[0].mz_error, mspeak_obj[0].to_string_formated)

def test_molecular_formula_search_db():
    
    mass_spec_obj = create_mass_spectrum()
    
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
                error.append(mformula._calc_assignment_mass_error(mspeak.mz_exp))
                abundance.append(mspeak.abundance)
        else:
            j += 1
            pass
    
    print('%i peaks assigned and %i peaks not assigned' % (i, j))

def test_priorityAssignment():
    
    
    MolecularSearchSettings.error_method = 'None'
    MolecularSearchSettings.min_mz_error = -3
    MolecularSearchSettings.max_mz_error = 3
    MolecularSearchSettings.mz_error_range = 1
    MolecularSearchSettings.isProtonated = True 
    MolecularSearchSettings.isRadical= False 
    MolecularSearchSettings.isAdduct= False 

    mass_spec_obj = create_mass_spectrum()
    
    assignOx = OxygenPriorityAssignment(mass_spec_obj) 

    assignOx.run() 
 
    

if __name__ == "__main__":

    #test_priorityAssignment()
    #test_molecular_formula_search_db()
    #test_run_molecular_formula_search()
    test_mspeak_search()
