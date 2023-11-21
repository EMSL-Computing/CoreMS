__author__ = "Yuri E. Corilo"
__date__ = "Jul 25, 2019"


import sys
sys.path.append('.')

import time
from pathlib import Path

import pytest

from corems.molecular_id.factory.classification import  HeteroatomsClassification
from corems.mass_spectrum.input.numpyArray import ms_from_array_centroid
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.transient.input.brukerSolarix import ReadBrukerSolarix
from corems.encapsulation.factory.parameters import MSParameters

def create_mass_spectrum():
    
    file_location = Path.cwd() / "tests/tests_data/ftms/ESI_NEG_SRFA.d/"

    bruker_reader = ReadBrukerSolarix(file_location)

    MSParameters.mass_spectrum.noise_threshold_method = 'log'
    MSParameters.mass_spectrum.noise_threshold_log_nsigma = 10
    MSParameters.ms_peak.peak_min_prominence_percent = 1

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
    
    #MSParameters.molecular_search.url_database = ''
    MSParameters.molecular_search.usedAtoms['F'] = (0,0)
    MSParameters.molecular_search.usedAtoms['P'] = (0,0)
    MSParameters.molecular_search.usedAtoms['Cl'] = (0,0)
    MSParameters.molecular_search.isAdduct = False
    MSParameters.molecular_search.isRadical = False
    
    MSParameters.molecular_search.used_atom_valences['P'] = 0
    MSParameters.molecular_search.used_atom_valences['F'] = 0
    MSParameters.molecular_search.used_atom_valences['Cl'] = 0

    MSParameters.mass_spectrum.noise_threshold_method = 'relative_abundance'
    MSParameters.mass_spectrum.noise_threshold_min_relative_abundance = 0.1
    
    mz = [215.09269]
    abundance = [1]
    rp, s2n = [1] ,[1]
    dataname = 'one peak'
    MSParameters.mass_spectrum.noise_threshold_method = 'relative_abundance'
    mass_spectrum_obj = ms_from_array_centroid(mz, abundance, rp, s2n, dataname, auto_process=False)

    SearchMolecularFormulas(mass_spectrum_obj).run_worker_ms_peaks([mass_spectrum_obj[0]])
    
    ms_peak = mass_spectrum_obj[0]
    print(ms_peak.mz_exp)
    if ms_peak.is_assigned:
        for formula in ms_peak:
            print(formula.string_formated, formula.mz_error)

def test_mspeak_search():

    mass_spec_obj = create_mass_spectrum()
    
    print("OK")

    mspeak_obj = mass_spec_obj.most_abundant_mspeak
    
    SearchMolecularFormulas(mass_spec_obj).run_worker_ms_peaks([mspeak_obj])

    print("OK2")
    if mspeak_obj.is_assigned:
        
        print(mspeak_obj.molecular_formula_earth_filter().string)
        print(mspeak_obj.molecular_formula_water_filter().string)
        print(mspeak_obj.molecular_formula_air_filter().string)
        print(mspeak_obj.cia_score_S_P_error().string)
        print(mspeak_obj.cia_score_N_S_P_error().string)
        print(mspeak_obj.best_molecular_formula_candidate.string)
        print(mspeak_obj[0].mz_error, mspeak_obj[0].string_formated)

def test_molecular_formula_search_db():
    
    MSParameters.molecular_search.isAdduct = False
    MSParameters.molecular_search.isRadical = False

    mass_spec_obj = create_mass_spectrum()
    
    time1 = time.time()
    
    SearchMolecularFormulas(mass_spec_obj, first_hit=True).run_worker_mass_spectrum()
    
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
                error.append(mformula.mz_error)
                abundance.append(mspeak.abundance)
        else:
            j += 1
            pass
    
    print('%i peaks assigned and %i peaks not assigned' % (i, j))

def test_priorityAssignment():
    
    MSParameters.molecular_search.error_method = 'None'
    MSParameters.molecular_search.min_ppm_error  = -3
    MSParameters.molecular_search.max_ppm_error = 5
    MSParameters.molecular_search.mz_error_range = 1
    MSParameters.molecular_search.isProtonated = True 
    MSParameters.molecular_search.isRadical= True 
    MSParameters.molecular_search.isAdduct= False 
    usedatoms = {'C': (1,100) , 'H': (4,200), 'O': (1,10)}
    MSParameters.molecular_search.usedAtoms = usedatoms
    
    mass_spec_obj = create_mass_spectrum()
    mass_spec_obj.process_mass_spec()

    assignOx = OxygenPriorityAssignment(mass_spec_obj) 

    assignOx.run() 

    #test classification 
    mass_spec_obj.percentile_assigned()

    mass_spectrum_by_classes = HeteroatomsClassification(mass_spec_obj)

    mass_spectrum_by_classes.plot_ms_assigned_unassigned()
    
    mass_spectrum_by_classes.atoms_ratio_all("H", "C")

    mass_spectrum_by_classes.atoms_ratio_all("H", "C")

if __name__ == "__main__":

    #test_priorityAssignment()
    #()
    #test_molecular_formula_search_db()
    test_run_molecular_formula_search()
    #test_mspeak_search()
