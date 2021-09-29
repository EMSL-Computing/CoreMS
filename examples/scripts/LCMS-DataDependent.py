



from ast import parse
from logging import warn
import re
from typing import Dict, List, Tuple
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

from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import LCMSParameters, MSParameters
from corems.molecular_formula.factory.MolecularFormulaFactory import MolecularFormula
from corems.molecular_formula.input.masslist_ref import ImportMassListRef
from corems.encapsulation.constant import Labels
from corems.mass_spectra.input import rawFileReader
from corems import get_dirname, get_filename

def run_thermo(file_location, target_mzs: List[float]) -> Tuple[Dict[float, rawFileReader.EIC_Data], rawFileReader.ImportDataDependentThermoMSFileReader]:
    
    print(file_location)
    
    LCMSParameters.lc_ms.start_scan = -1
    LCMSParameters.lc_ms.end_scan = -1

    LCMSParameters.lc_ms.smooth_window = 5
    
    LCMSParameters.lc_ms.min_peak_datapoints = 3
    LCMSParameters.lc_ms.peak_height_min_percent = 0.1

    LCMSParameters.lc_ms.eic_signal_threshold = 0.1
    LCMSParameters.lc_ms.eic_tolerance_ppm = 5
    LCMSParameters.lc_ms.enforce_target_ms2 = True
    LCMSParameters.lc_ms.average_target_mz = False
    
    parser = rawFileReader.ImportDataDependentThermoMSFileReader(file_location, target_mzs)

    tic_data, ax_tic = parser.get_tic(ms_type='MS !d', peak_detection=True, 
                                      smooth=True, plot=True)
    
    #ms2_tic, ax_ms2_tic = parser.get_tic(ms_type='MS2', peak_detection=False, plot=False)

    #print(data)

    # get selected data dependent mzs 
    target_mzs = parser.selected_mzs
    
    print(target_mzs)

    eics_data, ax_eic = parser.get_eics(target_mzs,
                                        tic_data,
                                        smooth=True,
                                        plot=True,
                                        legend=False,
                                        peak_detection=True,
                                        ax=ax_tic)
    
    plt.show()
    #ax_eic.plot(tic_data.time, tic_data.tic, c='black')
    return eics_data, parser
    

def read_lib(ref_filepath:Path):

    ion_charge = -1   

    iontypes = [Labels.protonated_de_ion]
    
    mf_references_dict = ImportMassListRef(ref_filepath).from_lcms_lib_file(ion_charge, iontypes)

    return mf_references_dict

def single_process(mf_references_dict: Dict[str, Dict[float, List[MolecularFormula]]], datapath: Path):

    #get mix name from filename
    current_mix = (re.findall(r'Mix[0-9]{1,2}', str(datapath)))[0]

    #get target compounds mz and molecular formulas
    dict_tarrget_mzs = mf_references_dict.get(current_mix)   
    
    target_mzs = dict_tarrget_mzs.keys()

    eics_data, parser = run_thermo(datapath, target_mzs)

    #need to convert this to a lcms object
    scan_number_mass_spectrum = {}
    
    # mz is from calculate mz
    for mz, eic_data in eics_data.items():
        
        #all possible m/z from the same mix, should be one per m/z as per current lib
        possible_mf = dict_tarrget_mzs.get(mz)

        if eic_data.apexes:                    
           
            print("m/z =  {}, formulas = {}, names = {},  peaks indexes = {}, retention times = {}, abundance = {}".format(mz,
                                                                                    [mf_obj.string for mf_obj in possible_mf],
                                                                                    [mf_obj.name for mf_obj in possible_mf],
                                                                                    eic_data.apexes,
                                                                                    [eic_data.time[apex[1]] for apex in eic_data.apexes],
                                                                                    [eic_data.eic[apex[1]] for apex in eic_data.apexes]) )

            for peak_index in eic_data.apexes:
           
                
                apex_index = peak_index[1]
                retention_time = eic_data.time[apex_index]
                original_scan = eic_data.scans[apex_index]
                
                parser.chromatogram_settings.start_scan = original_scan
                parser.chromatogram_settings.end_scan = original_scan
                
                mass_spec = parser.get_average_mass_spectrum_in_scan_range()
                if mass_spec:
                    
                    if original_scan not in scan_number_mass_spectrum.keys():
                        mass_spec.retention_time = retention_time
                        scan_number_mass_spectrum[original_scan] = [mass_spec, possible_mf]
                        mass_spec.plot_mz_domain_profile()
                        #plt.show()
                    
                    else:
                        scan_number_mass_spectrum[original_scan][1].extend(possible_mf)
    
    # TODO: create lcms and add dependent scans based on scan number 
    # Search molecular formulas on the mass spectrum, might need to use ProxyObject?
    # Add Adducts search, right now only working for de or protonated species
    
    ion_type = Labels.protonated_de_ion
    
    precision_decimals = 4

    for scan, ms_mf in scan_number_mass_spectrum.items():
        
        dependent_scans = parser.iRawDataPlus.GetScanDependents(scan, precision_decimals)

        mass_spectcrum_obj = ms_mf[0]
        mf_references_list = ms_mf[1]

        print()       

        ms_peaks_assigned = SearchMolecularFormulas(mass_spectcrum_obj).search_mol_formulas( mf_references_list, ion_type, find_isotopologues=True)
        
        for peak in ms_peaks_assigned:
            
            for mf in peak:
                if not mf.is_isotopologue:
                    print(mass_spectcrum_obj.retention_time, mf.name, mf.mz_calc, mf.mz_error, mf.confidence_score, mf.isotopologue_similarity)  
        
        if  ms_peaks_assigned:
            
            print([(i.ScanIndex, [i for i in i.PrecursorMassArray]) for i in dependent_scans.ScanDependentDetailArray])
        
        print()       
        #print(ms_peaks_assigned)
        

def auto_process(mf_references_dict: Dict[str, Dict[float, List[MolecularFormula]]], datadir: Path):

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
                eics_data, parser = run_thermo(file_path, target_mzs)

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
        single_process(mf_references_dict, file_location)
        #loop trought a directory, find and match mix name to raw file, 
        # search eic, do peak picking, for target compounds, currently enforcing parent ion selected to MS2, 
        # to change settings, chage  LCMSParameters.lc_ms parameters inside run_thermo() function:
        # TODO, search correspondent MS1 for m/z on peaks found, (sum 3 datapoints?)
        #       do molecular formula search on MS1 level, than get candidates and MF assiggment scores,
        # calculates fragment masses, and based on the highest score molecular formuyla calculate fragment formula 
        #auto_process(mf_references_dict, dirpath)
        #run_thermo(file_location)
