__author__ = "Yuri E. Corilo"
__date__ = "Aug 26, 2019"


import  sys, time, pytest, matplotlib
from pathlib import Path
sys.path.append(".")

import numpy as np
from matplotlib import pyplot 

from corems.encapsulation.settings.molecular_id.MolecularIDSettings import MolecularSearchSettings
from corems.mass_spectrum.calc.CalibrationCalc import FreqDomain_Calibration, MZDomain_Calibration
#from corems.mass_spectrum.input.massList import Read_MassList
from corems.molecular_id.search.findOxygenPeaks import FindOxygenPeaks
from corems.transient.input.BrukerSolarix import ReadBrukerSolarix
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
    
    #mass_spectrum_obj.plot_mz_domain_profile_and_noise_threshold()
    
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
    
    
    
    
    MolecularSearchSettings.error_method = 'None'
    MolecularSearchSettings.min_mz_error = -5
    MolecularSearchSettings.max_mz_error = 5
    MolecularSearchSettings.mz_error_range = 1
    MolecularSearchSettings.isProtonated = True 
    MolecularSearchSettings.isRadical= True 

    file_location = Path.cwd() /  "ESI_NEG_SRFA.d/"

    mass_spectrum = create_mass_spectrum(file_location)

    find_formula_thread = FindOxygenPeaks(mass_spectrum)
    find_formula_thread.start()
    find_formula_thread.join()
    
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

    MolecularSearchSettings.error_method = 'symmetrical'
    MolecularSearchSettings.min_mz_error = -3
    MolecularSearchSettings.max_mz_error = 3
    MolecularSearchSettings.mz_error_range = 1
    MolecularSearchSettings.mz_error_average = 0
    MolecularSearchSettings.min_abun_error = -30 # percentage 
    MolecularSearchSettings.max_abun_error = 70 # percentage 
    MolecularSearchSettings.isProtonated = True 
    MolecularSearchSettings.isRadical= True 
    
    MolecularSearchSettings.usedAtoms = {'C': (1, 100),
                 'H': (4, 200),
                 'O': (0, 20),
                 'N': (0, 1),
                 'S': (0, 0),
                 'P': (0, 0),
                 }
    
    #print(len(mass_spectrum))
    ClusteringFilter().filter_kendrick(mass_spectrum)
    #print(len(mass_spectrum))
   
    SearchMolecularFormulas().run_worker_mass_spectrum(mass_spectrum)
    ClusteringFilter().remove_assigment_by_mass_error(mass_spectrum)  

def test_import_ref_list():
    pass    

if __name__ == "__main__":
    
    from corems.encapsulation.settings.io import settings_parsers  

    MolecularSearchSettings.error_method = 'None'
    MolecularSearchSettings.min_mz_error = -7
    MolecularSearchSettings.max_mz_error = 0
    MolecularSearchSettings.mz_error_range = 1
    MolecularSearchSettings.isProtonated = True 
    MolecularSearchSettings.isRadical= False 
    MolecularSearchSettings.isAdduct= False 
    MolecularSearchSettings.usedAtoms['O'] = (1,20)
    
    print(MolecularSearchSettings.usedAtoms)

    #settings_parsers.load_search_setting_json(settings_path="SettingsCoreMS.json")    

    file_location = Path.cwd() /  "ESI_NEG_SRFA.d/"

    #file_location = Path("C:\\Users\\eber373\\OneDrive - PNNL\\Trabalhos\\Mayes\\Mayes_V1D76Alt_ICR_23Sept19_Alder_Infuse_p05_1_01_48741.d")

    mass_spectrum = create_mass_spectrum(file_location)
    
    print(mass_spectrum.polarity)
    
    find_formula_thread = FindOxygenPeaks(mass_spectrum)
    find_formula_thread.run()
    mspeaks_results = find_formula_thread.get_list_found_peaks()
    
    
    #mass_spectrum.plot_mz_domain_profile_and_noise_threshold()
    calibrate = FreqDomain_Calibration(mass_spectrum, mspeaks_results)
    calibrate.step_fit()

    mass_spectrum.molecular_search_settings.min_mz_error = -1
    mass_spectrum.molecular_search_settings.max_mz_error = 1
    
    find_formula_thread = FindOxygenPeaks(mass_spectrum)
    find_formula_thread.run()
    mspeaks_results = find_formula_thread.get_list_found_peaks()
    
    fig, ax = pyplot.subplots()
    

    ax.plot(mass_spectrum.mz_exp_profile, mass_spectrum.abundance_profile)
    for mspeak in mspeaks_results:
        if mspeak:
            for mf in mspeak:
                
                ax.plot(mspeak.mz_exp, mspeak.abundance, 'o', c='g')     
                ax.annotate(mspeak[0].to_string_formated,  (mspeak.mz_exp, mspeak.abundance + 100), fontsize=22)     
                #print(peak.mz_exp,mf.to_string, mf.mz_error )
    
   
    ax.label_outer()
    ax.get_yaxis().set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=16)

    pyplot.xlabel("m/z",fontsize=16)
    #pyplot.ylabel("DBE")
    pyplot.show()
    
    #test_calibration()
   