
import os
import sys
sys.path.append(".")

import pytest
from matplotlib import colors as mcolors
from matplotlib import pyplot

from enviroms.encapsulation.settings.molecular_id.MolecularIDSettings import  MoleculaLookupTableSettings, MoleculaSearchSettings
from enviroms.mass_spectrum.calc.CalibrationCalc import FreqDomain_Calibration, MZDomain_Calibration
from enviroms.mass_spectrum.output.MassSpecExport import MassSpecExport 
from enviroms.molecular_id.calc.FindOxigenPeaks import FindOxygenPeaks
from enviroms.molecular_id.calc.PrioriryAssignment import OxigenPriorityAssignment
from enviroms.molecular_id.calc.MolecularFormulaSearch import SearchMolecularFormulas
from enviroms.transient.input.BrukerSolarix import ReadBrukerSolarix


def calibrate(mass_spectrum_obj):
    
    MoleculaSearchSettings.error_method = 'average'
    MoleculaSearchSettings.min_mz_error = -5
    MoleculaSearchSettings.max_mz_error = 5
    MoleculaSearchSettings.mz_error_range = 1

    find_formula_thread = FindOxygenPeaks(mass_spectrum_obj, lookupTableSettings)
    find_formula_thread.run()
    mspeaks_results = find_formula_thread.get_list_found_peaks()
    
    calibrate = FreqDomain_Calibration(mass_spectrum_obj, mspeaks_results)
    calibrate.step_fit()
    mass_spectrum_obj.clear_molecular_formulas()

def assign_mf(mass_spectrum_obj):
    
    MoleculaSearchSettings.error_method = 'symmetrical'
    MoleculaSearchSettings.min_mz_error = -1
    MoleculaSearchSettings.max_mz_error = 1
    MoleculaSearchSettings.mz_error_range = 1
    MoleculaSearchSettings.mz_error_average = 0
    MoleculaSearchSettings.min_abun_error = -30 # percentage
    MoleculaSearchSettings.max_abun_error = 70 # percentage
    MoleculaSearchSettings.isProtonated = True
    MoleculaSearchSettings.isRadical = False
    MoleculaSearchSettings.isAdduct = True

    assignOx = OxigenPriorityAssignment(mass_spectrum_obj, lookupTableSettings)
    assignOx.start()
    assignOx.join()

def search_mf(mass_spectrum_obj):

    MoleculaSearchSettings.error_method = 'symmetrical'
    MoleculaSearchSettings.min_mz_error = -5
    MoleculaSearchSettings.max_mz_error = 5
    MoleculaSearchSettings.mz_error_range = 2
    MoleculaSearchSettings.mz_error_average = 0
    MoleculaSearchSettings.min_abun_error = -30 # percentage
    MoleculaSearchSettings.max_abun_error = 70 # percentage
    MoleculaSearchSettings.isProtonated = False
    MoleculaSearchSettings.isRadical = False
    MoleculaSearchSettings.isAdduct = True

    SearchMolecularFormulas().run_worker_mass_spectrum(mass_spectrum_obj, lookupTableSettings)

def plot():
    #colors = list(mcolors.XKCD_COLORS.keys())
    #oxigens = range(6,21)
    
    #for o in oxigens:
        #o_c = list()
    
    #pyplot.plot(mass_spectrum_obj.mz_exp_profile, mass_spectrum_obj.abundance_profile)
    
    pyplot.scatter(mass_spectrum_obj.mz_exp, mass_spectrum_obj.resolving_power,
                                         s=mass_spectrum_obj.signal_to_noise, 
                                         cmap='seismic')
    print(max(mass_spectrum_obj.signal_to_noise), min(mass_spectrum_obj.signal_to_noise))
    for mspeak in mass_spectrum_obj:
        
        if mspeak:
            #molecular_formula = mspeak.molecular_formula_lowest_error
            off_set = 0 
            for molecular_formula in mspeak:
                
                if not molecular_formula.is_isotopologue:
                    if molecular_formula.mspeak_indexes_isotopologues:
                        #pyplot.annotate(molecular_formula.to_string, (mspeak.mz_exp + off_set, mspeak.abundance ))
                        
                        off_set +=  0.1
                    #if molecular_formula['O'] == o:
                    #    if  not molecular_formula.is_isotopologue:
                    #        pyplot.plot(molecular_formula['C'], molecular_formula.dbe, "o",   color=colors[molecular_formula['O']])
                    #        pyplot.plot(molecular_formula['C'], molecular_formula.dbe, "o",   color=colors[molecular_formula['O']])
                    #        pyplot.annotate(molecular_formula.class_label, (molecular_formula['C']+0.5, molecular_formula.dbe+0.5))

    pyplot.show()                
    
if __name__ == "__main__":

    #file_location = os.path.join(os.getcwd(), "data/") + os.path.normcase("20190315_WK_CADY_O68_S1_PCP55A_H2O_000001.d")

    file_location = os.path.join(os.getcwd(), "tests/tests_data/") + os.path.normcase("ESI_NEG_SRFA.d/")
    
    bruker_reader = ReadBrukerSolarix(file_location)

    bruker_transient = bruker_reader.get_transient()

    mass_spectrum_obj = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=True)

    #mass_spectrum_obj.plot_mz_domain_profile_and_noise_threshold()

    lookupTableSettings = MoleculaLookupTableSettings()
    
    lookupTableSettings.usedAtoms['O'] = (0, 30)
    lookupTableSettings.usedAtoms['N'] = (0, 3)
    lookupTableSettings.usedAtoms['S'] = (0, 3)
    lookupTableSettings.usedAtoms['Cl'] = (1, 1)

    #calibrate(mass_spectrum_obj)

    #assign_mf(mass_spectrum_obj)

    #search_mf(mass_spectrum_obj)

    MassSpecExport(mass_spectrum_obj.filename, mass_spectrum_obj, 'excel').start()


    plot()
