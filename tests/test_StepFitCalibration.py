__author__ = "Yuri E. Corilo"
__date__ = "Aug 26, 2019"


import os, sys, time, pytest, matplotlib
sys.path.append(".")

import numpy as np
from matplotlib import pyplot as pylab

from corems.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaSearchSettings, MoleculaLookupDictSettings
from corems.mass_spectrum.calc.CalibrationCalc import FreqDomain_Calibration
#from corems.mass_spectrum.input.textMassList import Read_MassList
from corems.molecular_id.search.FindOxigenPeaks import FindOxygenPeaks
from corems.transient.input.BrukerSolarix import ReadBrukerSolarix
from corems.molecular_id.search.MolecularFormulaSearch import SearchMolecularFormulas
from corems.molecular_id.calc.ClusterFilter import ClusteringFilter
from corems.mass_spectrum.output.export import MassSpecExport

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

    file_location = directory + file_name

    mass_spectrum = creat_mass_spectrum(file_location)
    
    MoleculaSearchSettings.error_method = 'None'
    MoleculaSearchSettings.min_mz_error = -5
    MoleculaSearchSettings.max_mz_error = 5
    MoleculaSearchSettings.mz_error_range = 1

    find_formula_thread = FindOxygenPeaks(mass_spectrum)
    find_formula_thread.start()
    find_formula_thread.join()
    
    mspeaks_results = find_formula_thread.get_list_found_peaks()
    
    calibrate = FreqDomain_Calibration(mass_spectrum, mspeaks_results)
    calibrate.ledford_calibration()
    #calibrate.step_fit()
    mass_spectrum.clear_molecular_formulas()

    MoleculaSearchSettings.error_method = 'symmetrical'
    MoleculaSearchSettings.min_mz_error = -3
    MoleculaSearchSettings.max_mz_error = 3
    MoleculaSearchSettings.mz_error_range = 1
    MoleculaSearchSettings.mz_error_average = 0
    MoleculaSearchSettings.min_abun_error = -30 # percentage 
    MoleculaSearchSettings.max_abun_error = 70 # percentage 
    MoleculaSearchSettings.isProtonated = True 
    MoleculaSearchSettings.isRadical= True 
    
    MoleculaSearchSettings.usedAtoms = {'C': (1, 100),
                 'H': (4, 100),
                 'O': (0, 20),
                 'N': (0, 1),
                 'S': (0, 0),
                 'P': (0, 0),
                 }
    
    #print(len(mass_spectrum))
    ClusteringFilter().filter_kendrick(mass_spectrum)
    #print(len(mass_spectrum))
    time0 = time.time()
    print('started')
    SearchMolecularFormulas().run_worker_mass_spectrum(mass_spectrum)
    print(time.time()-time0)
    
    exportMS= MassSpecExport('neg_esi_srfa_1ppm_test', mass_spectrum)
    exportMS.run()
    
    exportMS.output_type = 'csv'
    exportMS.run()
    
    exportMS.output_type = 'pandas'
    exportMS.run()
    
    exportMS.get_pandas_df()

if __name__ == "__main__":
    
    test_calibration()
     