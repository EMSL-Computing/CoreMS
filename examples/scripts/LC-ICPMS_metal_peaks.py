##### CWD 2021-10-27
#### (1) import .csv with LC-ICPMS data; (2) for each metal, pick peaks; (3) return dictionary with RT as key, metal as value

from dataclasses import dataclass, field
from typing import Dict, List, Tuple
import warnings
import sys

from scipy.ndimage.measurements import label



sys.path.append("./")
warnings.filterwarnings("ignore")


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from corems.encapsulation.factory.parameters import LCMSParameters
from corems.mass_spectra.input import rawFileReader
from corems.mass_spectra.calc import SignalProcessing as sp

@dataclass
class TIC_Data:
     '''
    Scans: [int]
        original thermo scan numbers
    Time: [floats]
        list of retention times
    TIC: [floats]
        total ion chromatogram
    Apexes: [int]    
        original thermo apex scan number after peak picking 
     '''
     
     scans : List[int] = field(default_factory=list)
     time : List[float] = field(default_factory=list)
     tic : List[float] = field(default_factory=list)
     apexes : List[int] = field(default_factory=list)

@dataclass
class EIC_Data:
     '''
    Scans: [int]
        original thermo scan numbers
    Time: [floats]
        list of retention times
    EIC: [floats]
        extracted ion chromatogram
    Apexes: [int]    
        original thermo apex scan number after peak picking 
    
     '''
     metal: str
     scans : List[int] = field(default_factory=list)
     time : List[float] = field(default_factory=list)
     eic : List[float] = field(default_factory=list)
     apexes : List[int] = field(default_factory=list)

def eic_centroid_detector(max_tic, eic_data:EIC_Data, parameters:LCMSParameters, smooth=True):
    # Do peak picking and store results inside EIC_Data class    
	
	max_prominence = parameters.lc_ms.peak_max_prominence_percent

	max_height = parameters.lc_ms.peak_height_max_percent

	signal_threshold = parameters.lc_ms.eic_signal_threshold

	min_peak_datapoints = parameters.lc_ms.min_peak_datapoints

	correct_baseline = parameters.lc_ms.correct_eic_baseline
	
	peak_derivative_threshold = parameters.lc_ms.peak_derivative_threshold

	if smooth:
		
		eic_signal = smooth_signal(eic_data.eic, parameters)
	
	else:
	
		eic_signal = eic_data.eic

	peak_indexes_generator = sp.peak_picking_first_derivative(eic_data.time, eic_signal, max_height, max_prominence, 
															  max(eic_signal), min_peak_datapoints,
															  peak_derivative_threshold, 
															  signal_threshold=signal_threshold, 
															  correct_baseline=correct_baseline,
															  check_abundance=True,
															  plot_res=False,)
	eic_data.apexes = [i for i in peak_indexes_generator]
	
	#plt.plot(eic_data.time, eic_signal, label=eic_data.metal)
	
	#for peak_index in eic_data.apexes:
	#	plt.plot(eic_data.time[peak_index[1]], eic_signal[peak_index[1]], marker='x')
	#plt.legend()
	#plt.show()

def centroid_detector(tic_data:TIC_Data, parameters:LCMSParameters):
    # Do peak picking and store results inside TIC_Data class    
	noise_std = parameters.lc_ms.std_noise_threshold

	method = parameters.lc_ms.noise_threshold_method
	
	#peak picking
	min_height = parameters.lc_ms.peak_height_min_percent 
	min_datapoints = parameters.lc_ms.min_peak_datapoints   
	
	# baseline detection
	max_prominence = parameters.lc_ms.peak_max_prominence_percent 
	max_height = parameters.lc_ms.peak_height_max_percent 
	
	peak_indexes_generator = sp.peak_detector_generator(tic_data.tic, noise_std, method, tic_data.time, max_height, min_height, max_prominence, min_datapoints)

	tic_data.apexes = [i for i in peak_indexes_generator]
	

def smooth_signal(signal, parameters:LCMSParameters):
            
	implemented_smooth_method = parameters.lc_ms.implemented_smooth_method 
	
	pol_order = parameters.lc_ms.savgol_pol_order

	smooth_method = parameters.lc_ms.smooth_method

	window_len = parameters.lc_ms.smooth_window

	return sp.smooth_signal(signal, window_len, smooth_method, pol_order, implemented_smooth_method)

def get_data(icpdata: pd.DataFrame, parameters:LCMSParameters) -> Tuple[TIC_Data, Dict[str, EIC_Data]]:

	eic_metal_dict = {}

	tic_data = TIC_Data(time= [], scans= [], tic= [], apexes= [])
	
	tic = np.zeros(len(icpdata))

	scans = icpdata.Number

	for columns_label in icpdata.columns[1:]:
		
		if columns_label[0:4] == 'Time':
			
			rts = icpdata[columns_label]
			
			metal_label = columns_label.replace('Time ', '')

			eic_signal = icpdata[metal_label] 
			
			tic = tic + eic_signal

			eic_data = EIC_Data(metal= metal_label, time= rts, 
								scans= scans, eic= eic_signal, apexes= [])

			
			eic_metal_dict[metal_label] = eic_data

	tic_data.scans = scans
	tic_data.tic = list(tic)

	#fig, ax = plt.subplots()
	#ax.plot(tic_data.scans, smooth_signal(tic_data.tic, parameters))
	#plt.show()

	return tic_data, eic_metal_dict

def get_metal_data(icpfile):
	
	parameters = LCMSParameters()
	parameters.lc_ms.smooth_window = 301
	parameters.lc_ms.eic_signal_threshold = 1 #0-1
	parameters.lc_ms.min_peak_datapoints = 5
	parameters.lc_ms.correct_eic_baseline = False
	parameters.lc_ms.peak_max_prominence_percent = 0.1
	parameters.lc_ms.peak_height_max_percent = 1

	icpdata = pd.read_csv(icpfile, sep=',')
	
	tic_data, dict_metal_eicdata = get_data(icpdata, parameters)

	max_tic = max(tic_data.tic)
	
	for metal, eic_data in dict_metal_eicdata.items():

		eic_centroid_detector(max_tic, eic_data, parameters)

		print(metal, eic_data.apexes)
	
	return eic_data

def search_ms1_data(icrfile:str, eic_data:EIC_Data):
	
	'''place holder for parsing and search LC FT-MS data'''
	def run_thermo(file_location) -> Tuple[rawFileReader.DataDependentLCMS, rawFileReader.ImportDataDependentThermoMSFileReader]:
    
	
		LCMSParameters.lc_ms.start_scan = -1
		LCMSParameters.lc_ms.end_scan = -1

		#LCMSParameters.lc_ms.smooth_window = 3
		
		#LCMSParameters.lc_ms.min_peak_datapoints = 5
		#LCMSParameters.lc_ms.peak_height_min_percent = 1
		#LCMSParameters.lc_ms.peak_derivative_threshold = 1

		#LCMSParameters.lc_ms.peak_max_prominence_percent = 1
		#LCMSParameters.lc_ms.peak_height_max_percent = 1
		
		parser = rawFileReader.ImportDataDependentThermoMSFileReader(file_location)
		
		lcms_obj = parser.get_lcms_obj()
		
		return lcms_obj, parser

	lcms_obj, parser = run_thermo(icrfile)	
	
	tic_data, ax_tic = lcms_obj.get_tic(ms_type='MS !d', peak_detection=True, 
                                      smooth=False, plot=True)
    
	print(tic_data.apexes)

	plt.show()

	#eics_data, ax_eic = lcms_obj.get_eics(tic_data,
    #                                   smooth=True,
    #                                    plot=False,
    #                                    legend=False,
    #                                    peak_detection=True,
    #                                    ax=ax_tic)

if __name__ == '__main__':
	
	icpfile = "tests/tests_data/icpms/cwd_211018_day7_8_c18_1uMcobalamin_10uL.csv"
	
	icrfile = "tests/tests_data/icpms/rmb_161221_kansas_h2o_2.raw"
	
	eic_data = get_metal_data(icpfile)
	
	icr_data = search_ms1_data(icrfile, eic_data)

