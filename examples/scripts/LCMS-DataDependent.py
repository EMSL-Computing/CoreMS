



from logging import warn
from typing import Dict, List
import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.append("./")

from pathlib import Path
import glob

import matplotlib.pyplot as plt
import pandas as pd 
# from PySide2.QtWidgets import QFileDialog, QApplication
# from PySide2.QtCore import Qt

from corems.encapsulation.factory.parameters import LCMSParameters, MSParameters
from corems.molecular_formula.factory.MolecularFormulaFactory import MolecularFormula
from corems.molecular_formula.input.masslist_ref import ImportMassListRef
from corems.encapsulation.constant import Labels
from corems.mass_spectra.input import rawFileReader
from corems import get_dirname, get_filename

def run_thermo(file_location, target_mzs: List[float]):
    
    print(file_location)
    
    LCMSParameters.lc_ms.start_scan = -1
    LCMSParameters.lc_ms.end_scan = -1

    LCMSParameters.lc_ms.smooth_window = 5
    
    LCMSParameters.lc_ms.min_peak_datapoints = 3
    LCMSParameters.lc_ms.peak_height_min_percent = 0.1

    LCMSParameters.lc_ms.eic_signal_threshold = 0.1
    LCMSParameters.lc_ms.eic_tolerance_ppm = 5
    LCMSParameters.lc_ms.enforce_target_ms2 = False
    LCMSParameters.lc_ms.average_target_mz = False
    
    parser = rawFileReader.ImportDataDependentThermoMSFileReader(file_location, target_mzs)

    tic_data, ax_tic = parser.get_tic(ms_type='MS', peak_detection=True, 
                                      smooth=True, plot=False)

    ms2_tic, ax_ms2_tic = parser.get_tic(ms_type='MS2', peak_detection=False, plot=False)

    #print(data)

    # get selected data dependent mzs 
    target_mzs = parser.selected_mzs

    eics_data, ax_eic = parser.get_eics(target_mzs,
                                        tic_data,
                                        smooth=True,
                                        plot=False,
                                        legend=False,
                                        peak_detection=True,
                                        ax=ax_tic)
    

    #ax_eic.plot(tic_data.time, tic_data.tic, c='black')
    return eics_data
    

def read_lib(ref_filepath:Path):

    ion_charge = -1   
    iontype = Labels.protonated_de_ion
    
    mf_references_dict = ImportMassListRef(ref_filepath).from_lcms_lib_file(ion_charge, iontype)

    return mf_references_dict

def auto_process(mf_references_dict: Dict[str, Dict[float, List[MolecularFormula]]], datadir: Path):

    error = lambda x, ref: (x - ref)/ref * 1000000 
    
    #((self._mspeak_parent.mz_exp - self.mz_calc)/self.mz_calc)*multi_factor
    import os
    rootdir = datadir

    for mix_name in mf_references_dict.keys():
        
        file_locations = glob.glob(str(rootdir) + "/*.raw")
        file_paths = []
        
        for file_path in file_locations:
            if mix_name in file_path:
                file_paths.append(file_path)
        
        if file_paths:
            for file_path in file_paths:
                
                dict_tarrget_mzs = mf_references_dict.get(mix_name)    
                #mf_list = [mf for mf in mf_references_dict.get(mix_name)]

                target_mzs = dict_tarrget_mzs.keys()
                
                #print(target_mzs)
                eics_data = run_thermo(file_path, target_mzs)

                for mz, eic_data in eics_data.items():

                    if eic_data.apexes:
                        
                        possible_mf = dict_tarrget_mzs.get(mz)
                            
                        print("m/z =  {}, formulas = {}, names = {},  peaks indexes = {}, retention times = {}, abundance = {}".format(mz,
                                    [mf_obj.string for mf_obj in possible_mf],
                                    [mf_obj.name for mf_obj in possible_mf],
                                    eic_data.apexes,
                                    [eic_data.time[apex[1]] for apex in eic_data.apexes],
                                    [eic_data.eic[apex[1]] for apex in eic_data.apexes]) )
                            #print(mz, eic_data.apexes, [eic_data.time[apex[1]] for apex in eic_data.apexes])
                print()
                print()
                #print(parser.get_all_filters())
                #plt.show()

        else:
            warn("Could not find a any raw data with mix name: {}".format(mix_name))

if __name__ == "__main__":
    
    dirpath = "/Users/eber373/OneDrive - PNNL/Documents/Data/LCMS/RAW Files/C18/1st run/NEG/"
    ref_dirpath = "/Users/eber373/OneDrive - PNNL/Documents/Data/LCMS/"
    filename = "LCMS_5191_CapDev_C18_Mix1_NEG_28Apr2021.raw"
    ref_filename = "LCMS_StantardLibrary.csv"
    # run_multiprocess()
    # cpu_percents = monitor(target=run_multiprocess)
    # print(cpu_percents)
    # file_location = get_filename()
    file_location = Path(dirpath + filename)
    ref_file_location = Path(ref_dirpath + ref_filename)
    
    if file_location:
        
        #get the list of molecular formulas objects for each stardard compound maped by mix name
        mf_references_dict = read_lib(ref_file_location)
        
        #loop trought a directory, find and match mix name to raw file, 
        # search eic, do peak picking, for target compounds, currently enforcing parent ion selected to MS2, 
        # to change settings, chage  LCMSParameters.lc_ms parameters inside run_thermo() function:
        # TODO, search correspondent MS1 for m/z on peaks found, (sum 3 datapoints?)
        #       do molecular formula search on MS1 level, than get candidates and MF assiggment scores,
        # calculates fragment masses, and based on the highest score molecular formuyla calculate fragment formula 
        auto_process(mf_references_dict, dirpath)
        #run_thermo(file_location)
