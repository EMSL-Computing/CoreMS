__author__ = "Yuri E. Corilo"
__date__ = "Feb 12, 2020"

from pathlib import Path

from netCDF4 import Dataset

from corems.encapsulation.constant import Labels
from corems.encapsulation.settings.input import InputSetting
from corems.mass_spectra.factory.GC_Class import GCMSBase

from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecCentroid

class ReadAndiNetCDF:

	def __init__(self, file_location, analyzer='Quadruple', instrument_label='GCMS-Agilent', auto_process=True):
		
		file_location = Path(file_location)

		if not file_location.exists:

			raise FileNotFoundError("File does not exist at %s", file_location)

		self.file_location = file_location

		self.net_cdf_obj = Dataset(file_location, "r+", format='NETCDF3_CLASSIC')

		self.ionization_type = self.net_cdf_obj.test_ionization_mode

		self.experiment_type = self.net_cdf_obj.experiment_type

		self.list_scans = self.net_cdf_obj.variables.get("actual_scan_number")[:]

		self.initial_scan_number = self.list_scans[0]

		self.final_scan_number = self.list_scans[-1]

		self.analyzer = analyzer

		self.instrument_label = instrument_label
		
		self.gcms = GCMSBase(self.file_location, analyzer, instrument_label)

		#other unused variables
		#a_d_sampling_rate = self.net_cdf_obj.variables.get("a_d_sampling_rate")[:]
		#a_d_coaddition_factor = self.net_cdf_obj.variables.get("a_d_coaddition_factor")[:]
		#inter_scan_time = self.net_cdf_obj.variables.get("inter_scan_time")[:]
		#flag_count = self.net_cdf_obj.variables.get("flag_count")[:]
		#instrument_name = ''.join(self.net_cdf_obj.variables.get("instrument_name")[:].tolist())
		#instrument_id = self.net_cdf_obj.variables.get("instrument_id")[:]
		#instrument_mfr = self.net_cdf_obj.variables.get("instrument_mfr")[:]
		#instrument_model = self.net_cdf_obj.variables.get("instrument_model")[:]
		#instrument_serial_no = self.net_cdf_obj.variables.get("instrument_serial_no")[:]
		#instrument_sw_version = self.net_cdf_obj.variables.get("instrument_sw_version")[:]
		#instrument_fw_version = self.net_cdf_obj.variables.get("instrument_fw_version")[:]
		#instrument_os_version = self.net_cdf_obj.variables.get("instrument_os_version")[:]
		#instrument_app_version = self.net_cdf_obj.variables.get("instrument_app_version")[:]
		#instrument_comments = self.net_cdf_obj.variables.get("instrument_comments")[:]
		#scan_duration = self.net_cdf_obj.variables.get("scan_duration")[:]
		#time_range_min = self.net_cdf_obj.variables.get("time_range_min")[:]
		#time_range_max = self.net_cdf_obj.variables.get("time_range_max")[:]
		#scan_index_list = self.net_cdf_obj.variables.get("scan_index")[:]
		#time_values = self.net_cdf_obj.variables.get("time_values")[:]
		
	@property

	def polarity(self):

		polarity = str(self.net_cdf_obj.test_ionization_polarity)

		if polarity == 'Positive Polarity': return +1
			
		else: return -1    

	def get_mass_spectrum(self, mz, abun, rp, d_params):
				
		data_dict = {Labels.mz: mz,
				Labels.abundance: abun,
				Labels.rp: rp,
				Labels.s2n: None,
		}
            
		mass_spec = MassSpecCentroid(data_dict, d_params)

		return mass_spec

	def run(self):
        
		'''creates the gcms obj'''

		d_parameters = InputSetting.d_params(self.file_location)
		self.import_mass_spectra(d_parameters)

	def import_mass_spectra(self, d_params):
        
		ms_datapoints_per_scans = self.net_cdf_obj.variables.get("point_count")[:]

		list_tic = self.net_cdf_obj.variables.get("total_intensity")[:]	

		list_rt = self.net_cdf_obj.variables.get("scan_acquisition_time")[:]

		mass_values = self.net_cdf_obj.variables.get("mass_values")[:]
		intensity_values = self.net_cdf_obj.variables.get("intensity_values")[:]
		resolution = self.net_cdf_obj.variables.get("resolution")[:]

		finish_location = -1
		
		for scan_index in self.list_scans:

			datapoints = ms_datapoints_per_scans[scan_index]

			finish_location += datapoints 

			start_location = finish_location - datapoints + 1

			#finish_location = ms_datapoints_per_scans[start_location:last_pos]
						
			d_params["rt"] = list_rt[scan_index]

			d_params["scan_number"] = scan_index

			d_params['label'] = Labels.gcms_centroid

			d_params["polarity"] = self.polarity

			d_params['analyzer'] = self.analyzer

			d_params['instrument_label'] = self.instrument_label

			mz = mass_values[start_location:finish_location]
			
			abundance = intensity_values[start_location:finish_location]

			rp = resolution[start_location:finish_location]

			ms = self.get_mass_spectrum	
