
import os
import sys
import time
sys.path.append(".")
import numpy as np
from matplotlib import pyplot as pylab
import matplotlib
from numpy import average

from enviroms.emsl.yec.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaSearchSettings
from enviroms.emsl.yec.mass_spectrum.calc.CalibrationCalc import MZDomain_Calibration, FreqDomain_Calibration
from enviroms.emsl.yec.mass_spectrum.input.TextMassList import Read_MassList
from enviroms.emsl.yec.molecular_id.calc.FindOxigenPeaks import FindOxygenPeaks
from enviroms.emsl.yec.transient.input.BrukerSolarix import ReadBrukerSolarix
from enviroms.emsl.yec.molecular_id.calc.MolecularFormulaSearch import SearchMolecularFormulas
from enviroms.emsl.yec.molecular_id.calc.ClusterFilter import ClusteringFilter

def creat_mass_spectrum(file_location):

    bruker_reader = ReadBrukerSolarix(file_location)

    bruker_transient = bruker_reader.get_transient()

    mass_spectrum_obj = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=True)
    #mass_spectrum_obj.plot_mz_domain_profile_and_noise_threshold()
    
    # polariy need to be set if reading a text file
    #polariy = -1
    # load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    #mass_list_reader = Read_MassList(file_location, polariy, delimiter="  ")
    #mass_spectrum_obj  = mass_list_reader.get_mass_spectrum(auto_process=True)

    return mass_spectrum_obj

if __name__ == "__main__":
    
    # MoleculaLookupTableSettings and MoleculaSearchSettings at
    # enviroms\emsl\yec\encapsulation\settings\molecular_id\MolecularIDSettings.py
    # for changing settings of the lookup table and searching algorithms

    directory = os.path.join(os.getcwd(), "data/")

    #file_name = os.path.normcase("20190616_WK_ESFA_0pt2mgml_ESI_Neg_0pt8sFID_000001.d/")

    file_name = os.path.normcase("20190205_WK_SRFA_opt_000001.d/")

    #file_name = "20190616_WK_ESFA_0pt2mgml_ESI_Neg_1pt4sFID_000001.ascii"

    file_location = directory + file_name

    mass_spectrum = creat_mass_spectrum(file_location)
    
    time1 = time.time()

    find_formula_thread = FindOxygenPeaks(mass_spectrum)
    find_formula_thread.run()
    mspeaks_results = find_formula_thread.get_list_found_peaks()
    
    calibrate = FreqDomain_Calibration(mass_spectrum, mspeaks_results)
    #calibrate.ledford_calibration()
    calibrate.step_fit()
    mass_spectrum.clear_molecular_formulas()

    MoleculaSearchSettings.error_method = 'symmetrical'
    MoleculaSearchSettings.min_mz_error = -3
    MoleculaSearchSettings.max_mz_error = 1
    MoleculaSearchSettings.mz_error_range = 1
    MoleculaSearchSettings.mz_error_average = 0
    
    #ClusteringFilter().filter_kendrick(mass_spectrum)
    SearchMolecularFormulas().run_worker_mass_spectrum(mass_spectrum)
    
    error = list()
    mass = list()
    abundance = list()
    freq_exp = list()
    mz_theo = list()

    for mspeak in mass_spectrum:
        
        if mspeak:
            
            #molecular_formula = mspeak.molecular_formula_lowest_error
            for molecular_formula in mspeak:
                
            #if not molecular_formula.is_isotopologue:
                freq_exp.append(mspeak.freq_exp)
                mass.append(mspeak.mz_exp)
                error.append(molecular_formula._calc_assigment_mass_error(mspeak.mz_exp))
                abundance.append(mspeak.abundance)
                mz_theo.append(molecular_formula.mz_theor)

    #print(np.average(error), np.std(error))
    pylab.plot(mass_spectrum.mz_exp, mass_spectrum.abundance) 
    pylab.plot(mass, abundance, "o")  
    pylab.show()  
    pylab.plot(mass, error, "o")  
    pylab.show()  
