__author__ = "Yuri E. Corilo"
__date__ = "Aug 26, 2019"


import  sys, time, pytest, matplotlib
from pathlib import Path
sys.path.append(".")

import numpy as np
from matplotlib import pyplot 

from corems.encapsulation.factory.parameters import MSParameters
from corems.mass_spectrum.calc.CalibrationCalc import FreqDomain_Calibration, MZDomain_Calibration
#from corems.mass_spectrum.input.massList import Read_MassList
from corems.molecular_id.search.findOxygenPeaks import FindOxygenPeaks
from corems.transient.input.brukerSolarix import ReadBrukerSolarix
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.molecular_id.calc.ClusterFilter import ClusteringFilter

def create_mass_spectrum(file_location):
    
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
    
    #mass_spectrum_obj.plot_profile_and_noise_threshold()
    
    # polarity need to be set if reading a text file
    #polarity = -1
    # load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    #mass_list_reader = Read_MassList(file_location, polarity,  )
    #mass_spectrum_obj  = mass_list_reader.get_mass_spectrum(auto_process=True)

    return mass_spectrum_obj



def test_calibration():
    
    ''' Mass calibration test module: 
            - creates a mass spectrum object
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
    
    MSParameters.molecular_search.error_method = 'None'
    MSParameters.molecular_search.min_ppm_error  = -5
    MSParameters.molecular_search.max_ppm_error = 5
    MSParameters.molecular_search.mz_error_range = 1
    MSParameters.molecular_search.isProtonated = True 
    MSParameters.molecular_search.isRadical= True 

    file_location = Path.cwd() / "tests/tests_data/ESI_NEG_SRFA.d/"

    mass_spectrum = create_mass_spectrum(file_location)

    find_formula_thread = FindOxygenPeaks(mass_spectrum)
    find_formula_thread.run()
    #find_formula_thread.join()
    
    mspeaks_results = find_formula_thread.get_list_found_peaks()
    
    calibrate = FreqDomain_Calibration(mass_spectrum, mspeaks_results)
    calibrate.linear()
    calibrate.step_fit()
    calibrate.quadratic(iteration=True)
    calibrate.ledford_calibration()
    
    calibrate = MZDomain_Calibration(mass_spectrum, mspeaks_results, include_isotopologue=True)
    calibrate.linear(iteration=True)
    calibrate.ledford_inverted_calibration(iteration=True)
    calibrate.quadratic(iteration=True)
    
    mass_spectrum.clear_molecular_formulas()

    MSParameters.molecular_search.error_method = 'symmetrical'
    MSParameters.molecular_search.min_ppm_error  = -3
    MSParameters.molecular_search.max_ppm_error = 3
    MSParameters.molecular_search.mz_error_range = 1
    MSParameters.molecular_search.mz_error_average = 0
    MSParameters.molecular_search.min_abun_error = -30 # percentage 
    MSParameters.molecular_search.max_abun_error = 70 # percentage 
    MSParameters.molecular_search.isProtonated = True 
    MSParameters.molecular_search.isRadical= True 
    
    MSParameters.molecular_search.usedAtoms = {'C': (1, 100),
                 'H': (4, 200),
                 'O': (0, 20),
                 'N': (0, 1),
                 'S': (0, 0),
                 'P': (0, 0),
                 }
    
    #print(len(mass_spectrum))
    ClusteringFilter().filter_kendrick(mass_spectrum)
    #print(len(mass_spectrum))
   
    SearchMolecularFormulas(mass_spectrum).run_worker_mass_spectrum()
    ClusteringFilter().remove_assignment_by_mass_error(mass_spectrum)  

def test_import_ref_list():
    pass    

if __name__ == "__main__":
    
    test_calibration()
   