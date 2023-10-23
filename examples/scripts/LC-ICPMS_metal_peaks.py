##### CWD 2021-10-27
#### (1) import .csv with LC-ICPMS data; (2) for each metal, pick peaks; (3) return dictionary with RT as key, metal as value

from dataclasses import dataclass, field
from typing import Dict, List, Tuple
import warnings
import sys

sys.path.append("./")
warnings.filterwarnings("ignore")

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from corems.encapsulation.factory.parameters import LCMSParameters
from corems.mass_spectra.input import rawFileReader
from corems.mass_spectra.calc import SignalProcessing as sp
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas

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
	
	plt.plot(eic_data.time/60, eic_signal, label=eic_data.metal)
	
	for peak_index in eic_data.apexes:
		plt.plot(eic_data.time[peak_index[1]]/60, eic_signal[peak_index[1]], marker='x')
	
	plt.xlabel('Retention Time (min)')
	plt.ylabel('Intensity (cps)')
	plt.legend()
	plt.show()

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
			
			rts = icpdata[columns_label].to_numpy()
			
			metal_label = columns_label.replace('Time ', '')

			eic_signal = icpdata[metal_label] 
			
			tic = tic + eic_signal

			eic_data = EIC_Data(metal= metal_label, time= rts, 
								scans= scans, eic= eic_signal, apexes= [])

			
			eic_metal_dict[metal_label] = eic_data

	tic_data.scans = scans
	tic_data.tic = list(tic)

	#fig, ax = plt.subplots()
	#ax.plot(rts, smooth_signal(tic_data.tic, parameters))
	plt.xlabel('Retention Time')
	plt.ylabel('Intensity (cps)')
	plt.show()

	return tic_data, eic_metal_dict

def get_metal_data(icpfile, parameters:LCMSParameters):
	
	icpdata = pd.read_csv(icpfile, sep=',')
	
	tic_data, dict_metal_eicdata = get_data(icpdata, parameters)

	max_tic = max(tic_data.tic)
	
	for metal, eic_data in dict_metal_eicdata.items():

		eic_centroid_detector(max_tic, eic_data, parameters)

		#print(metal, eic_data.apexes)
	
	return dict_metal_eicdata

def search_ms1_data(icrfile:str, dict_metal_eicdata:  Dict[str, EIC_Data], parameters:LCMSParameters):
	
	'''place holder for parsing and search LC FT-MS data'''
	
	lcms_obj, parser = run_thermo(icrfile, parameters)	
	
	tic_data, ax_tic = lcms_obj.get_tic(ms_type='MS !d', peak_detection=True, 
                                      smooth=False, plot=True)
	
	plt.show()

	for metal, eic_data in dict_metal_eicdata.items():
		
		print(metal, eic_data.apexes)
		
		for peak_indexex in eic_data.apexes:
			
			ftms_scans_index = ([find_nearest_scan(eic_data.time[i], tic_data) for i in peak_indexex])
			ftms_scans = [tic_data.scans[i] for i in ftms_scans_index]
			ftms_times = [tic_data.time[i] for i in ftms_scans_index]
			
			retention_time = tic_data.time[ftms_scans_index[1]]

			print(ftms_scans)
			print(ftms_times)

			parser.chromatogram_settings.scans = (ftms_scans[0], ftms_scans[-1])
			
			mass_spec = parser.get_average_mass_spectrum(auto_process=False)
			mass_spec.retention_time = retention_time

			mass_spec.settings = parameters.mass_spectrum
			mass_spec.molecular_search_settings = parameters.ms1_molecular_search
			mass_spec.mspeaks_settings = parameters.ms_peak
			mass_spec.process_mass_spec()

			metal_atom = ''.join(i for i in metal if not i.isdigit())
			mass_spec.molecular_search_settings.usedAtoms[metal_atom] = (1,1)
			
			ax = mass_spec.plot_profile_and_noise_threshold()

			SearchMolecularFormulas(mass_spec, first_hit=False).run_worker_mass_spectrum()
			mass_spec.molecular_search_settings.usedAtoms[metal_atom] = (0,0)

			mass_spec.percentile_assigned(report_error=True)
			print(metal)
			filename = '{}_rt{}_{}'.format(metal, retention_time, mass_spec.sample_name).replace(".", "_")
			print(filename)
			mass_spec.to_csv(filename, write_metadata=False)
			
			for peak in mass_spec:
            
				for mf in peak:
					is_assigned = True
					
					annotation = "Mol. Form = {}\nm\z = {:.4f}\nerror = {:.4f}\nconfidence score = {:.2f}\nisotopologue score = {:.2f}".format(mf.string_formated, peak.mz_exp, mf.mz_error, mf.confidence_score, mf.isotopologue_similarity)
					
					ax.annotate(annotation , xy=(peak.mz_exp, peak.abundance),
												xytext=(+3, np.sign(peak.abundance)*-40), textcoords="offset points",
												horizontalalignment="left",
												verticalalignment="bottom" if peak.abundance > 0 else "top")
			plt.show()
			#time_range.append([eic_data.time[i]/60 for i in peak_indexex])



def run_thermo(file_location, parameters:LCMSParameters) -> Tuple[rawFileReader.DataDependentLCMS, rawFileReader.ImportDataDependentThermoMSFileReader]:
    
	#LCMSParameters.lc_ms.smooth_window = 3
	
	#LCMSParameters.lc_ms.min_peak_datapoints = 5
	#LCMSParameters.lc_ms.peak_height_min_percent = 1
	#LCMSParameters.lc_ms.peak_derivative_threshold = 1

	#LCMSParameters.lc_ms.peak_max_prominence_percent = 1
	#LCMSParameters.lc_ms.peak_height_max_percent = 1
	
	parser = rawFileReader.ImportDataDependentThermoMSFileReader(file_location)
	parser.parameters = parameters
	
	parser.find_nearest_scan()
	
	lcms_obj = parser.get_lcms_obj()
	
	return lcms_obj, parser


def find_nearest_scan(rt, ftms_data):

	lst = np.asarray(ftms_data.time)
    
	idx = (np.abs(lst - (rt/60))).argmin()
	
	return idx

if __name__ == '__main__':
	
	#icpfile = "tests/tests_data/icpms/cwd_211018_day7_8_c18_1uMcobalamin_10uL.csv"
	icpfile = 'tests/tests_data/icpms/161220_soils_hypercarb_3_kansas_qH2O.csv'
	icrfile = "tests/tests_data/icpms/rmb_161221_kansas_h2o_2.raw"
	
	parameters = LCMSParameters()

	parameters.lc_ms.smooth_window = 301
	parameters.lc_ms.eic_signal_threshold = 1 #0-1
	parameters.lc_ms.min_peak_datapoints = 5
	parameters.lc_ms.correct_eic_baseline = False
	parameters.lc_ms.peak_max_prominence_percent = 0.1
	parameters.lc_ms.peak_height_max_percent = 1

	parameters.mass_spectrum.threshold_method = 'log'
	parameters.mass_spectrum.noise_threshold_std = 1

	parameters.ms1_molecular_search.error_method = 'None'
	parameters.ms1_molecular_search.min_ppm_error = -1
	parameters.ms1_molecular_search.max_ppm_error = 1
	
	parameters.ms1_molecular_search.usedAtoms['C'] = (10, 100)
	parameters.ms1_molecular_search.usedAtoms['H'] = (4, 200)
	parameters.ms1_molecular_search.usedAtoms['O'] = (1, 30)
	parameters.ms1_molecular_search.usedAtoms['N'] = (1, 12)
	parameters.ms1_molecular_search.usedAtoms['S'] = (0, 0)

	parameters.ms1_molecular_search.isProtonated = True
	parameters.ms1_molecular_search.isRadical = False
	parameters.ms1_molecular_search.isAdduct = False

	dict_metal_eicdata = get_metal_data(icpfile, parameters)
	
	icr_data = search_ms1_data(icrfile, dict_metal_eicdata, parameters)

