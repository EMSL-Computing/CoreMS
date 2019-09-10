
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
from enviroms.transient.input.BrukerSolarix import ReadBrukerSolarix


def calibrate(mass_spectrum_obj):
    
    MoleculaSearchSettings.error_method = 'average'
    MoleculaSearchSettings.min_mz_error = -5
    MoleculaSearchSettings.max_mz_error = 1
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

def plot():
    colors = list(mcolors.XKCD_COLORS.keys())
    oxigens = range(6,21)
    for o in oxigens:
        #o_c = list()
        for mspeak in mass_spectrum_obj:
            if mspeak:
                #molecular_formula = mspeak.molecular_formula_lowest_error
                for molecular_formula in mspeak:
                    if molecular_formula['O'] == o:
                        if  not molecular_formula.is_isotopologue:
                            pyplot.plot(molecular_formula['C'], molecular_formula.dbe, "o",   color=colors[molecular_formula['O']])
                            pyplot.plot(molecular_formula['C'], molecular_formula.dbe, "o",   color=colors[molecular_formula['O']])
                            pyplot.annotate(molecular_formula.class_label, (molecular_formula['C']+0.5, molecular_formula.dbe+0.5))

        pyplot.show()                
    
if __name__ == "__main__":

    file_location = os.path.join(os.getcwd(), "tests/tests_data/") + os.path.normcase("ESI_NEG_SRFA.d/")
    
    bruker_reader = ReadBrukerSolarix(file_location)

    bruker_transient = bruker_reader.get_transient()

    mass_spectrum_obj = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=True)

    lookupTableSettings = MoleculaLookupTableSettings()
    
    lookupTableSettings.usedAtoms['O'] = (1, 22)
    
    calibrate(mass_spectrum_obj)

    assign_mf(mass_spectrum_obj)

    MassSpecExport('neg_esi_srfa_1ppm_test', mass_spectrum_obj, 'excel').start()

    plot()
