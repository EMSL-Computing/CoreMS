__author__ = "Yuri E. Corilo"
__date__ = "Aug 26, 2019"

import os, sys, time, pytest
sys.path.append(".")
import numpy as np
from matplotlib import pyplot as pylab
import matplotlib
from numpy import average

from enviroms.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaSearchSettings, MoleculaLookupTableSettings
from enviroms.mass_spectrum.calc.CalibrationCalc import MZDomain_Calibration, FreqDomain_Calibration
from enviroms.mass_spectrum.input.TextMassList import Read_MassList
from enviroms.molecular_id.calc.FindOxigenPeaks import FindOxygenPeaks
from enviroms.transient.input.BrukerSolarix import ReadBrukerSolarix
from enviroms.molecular_id.calc.MolecularFormulaSearch import SearchMolecularFormulas
from enviroms.molecular_id.calc.ClusterFilter import ClusteringFilter

def creat_mass_spectrum(file_location):
    '''parse transient data from Bruker into a mass spectrum class object

        Parameters
        ----------
        file_location: str
            The full path of the *.d data folder 
        
        Returns
        -------
        MassSpecfromFreq() class
           (See MassSpecfromFreq class for more details)  
        '''

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

def test_calibration():
    
    ''' Mass calibration test module: 
            - creates a mass spectrum objec
            - find oxygen most abundant peaks separated by 14Da
            - calibrate on frequency domain using ledford equation
            - filter data based on kendrick mass with CH2O base
            - search for all molecular formula candidates 

        Returns
        -------
        Nothing
            
            Store the results inside the mass spectrum class 
            (See Docs for the structural details)  
    '''
    
    directory = os.path.join(os.getcwd(), "tests/tests_data/")

    file_name = os.path.normcase("ESI_NEG_SRFA.d/")

    #file_name = "20190616_WK_ESFA_0pt2mgml_ESI_Neg_1pt4sFID_000001.ascii"

    file_location = directory + file_name

    mass_spectrum = creat_mass_spectrum(file_location)
    
    LookupTableSettings = MoleculaLookupTableSettings()

    MoleculaSearchSettings.error_method = 'average'
    MoleculaSearchSettings.min_mz_error = -5
    MoleculaSearchSettings.max_mz_error = 1
    MoleculaSearchSettings.mz_error_range = 1

    find_formula_thread = FindOxygenPeaks(mass_spectrum, LookupTableSettings)
    find_formula_thread.run()
    mspeaks_results = find_formula_thread.get_list_found_peaks()
    
    calibrate = FreqDomain_Calibration(mass_spectrum, mspeaks_results)
    calibrate.ledford_calibration()
    #calibrate.step_fit()
    mass_spectrum.clear_molecular_formulas()

    MoleculaSearchSettings.error_method = 'symmetrical'
    MoleculaSearchSettings.min_mz_error = -1
    MoleculaSearchSettings.max_mz_error = 1
    MoleculaSearchSettings.mz_error_range = 2
    MoleculaSearchSettings.mz_error_average = 0
    MoleculaSearchSettings.min_abun_error = -30 # percentage 
    MoleculaSearchSettings.max_abun_error = 70 # percentage 
    MoleculaSearchSettings.isProtonated = True 
    MoleculaSearchSettings.isRadical= True 
    
    LookupTableSettings.usedAtoms = {'C': (1, 100),
                 'H': (4, 200),
                 'O': (1, 20),
                 'N': (0, 0),
                 'S': (0, 0),
                 'P': (0, 0),
                 }
    
    #print(len(mass_spectrum))
    ClusteringFilter().filter_kendrick(mass_spectrum)
    #print(len(mass_spectrum))
    SearchMolecularFormulas().run_worker_mass_spectrum(mass_spectrum, LookupTableSettings)
    
    error = list()
    error_iso = list()
    mass_iso = list()
    mass = list()
    abundance = list()
    abundance_iso = list()
    freq_exp = list()
    mz_theo = list()
    o_c = list()
    h_c = list()
    colors = list(matplotlib.colors.XKCD_COLORS.keys())
    colors_oxigen = []
    oxigens = range(8,21)
    for o in oxigens:
        o_c = list()
        for mspeak in mass_spectrum:
            
            if mspeak:
                
                #molecular_formula = mspeak.molecular_formula_lowest_error
                for molecular_formula in mspeak:
                    if molecular_formula['O'] == o:
                        if  not molecular_formula.is_isotopologue:
                            freq_exp.append(mspeak.freq_exp)
                            mass.append(mspeak.mz_exp)
                            error.append(molecular_formula._calc_assigment_mass_error(mspeak.mz_exp))
                            abundance.append(mspeak.abundance)
                            mz_theo.append(molecular_formula.mz_theor)
                            o_c.append(molecular_formula.O_C)
                            h_c.append(molecular_formula.H_C)
                            pylab.plot(molecular_formula['C'], molecular_formula.dbe, "o",   color=colors[molecular_formula['O']])
                            pylab.annotate(molecular_formula['O'], (molecular_formula['C']+0.5, molecular_formula.dbe+0.5))

                        else:
                        
                            mass_iso.append(mspeak.mz_exp)
                            abundance_iso.append(mspeak.abundance)
                            error_iso.append(molecular_formula._calc_assigment_mass_error(mspeak.mz_exp))
        print(max(o_c), min(o_c))
        pylab.show()                
        
    if __name__ == "__main__":
        #don not plot if running as unit test
        
        print(np.average(error), np.std(error), (len(error)+len(error_iso))/len(mass_spectrum)*100)
        #pylab.plot(mass_spectrum.mz_exp, mass_spectrum.abundance)
        #pylab.plot(mass, abundance, "o") 
        #pylab.plot(mass_iso, abundance_iso, "o", color=colors_oxigen)  
        pylab.show()  
        pylab.plot(mass, error, "o")  
        pylab.plot(mass_iso, error_iso, "o", color='red')  
        pylab.show()  
        pylab.plot( o_c, h_c, "o", color='red')  
        pylab.show()  

if __name__ == "__main__":
    
    test_calibration()
     