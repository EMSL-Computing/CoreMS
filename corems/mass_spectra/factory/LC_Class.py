"""
Created on Jun 12, 2019
"""
from collections.abc import Mapping
from pathlib import Path
from corems.mass_spectra.calc.LC_Calc import LC_Calculations

__author__ = "Yuri E. Corilo"
__date__ = "Jun 25, 2019"

class LCMSBase(Mapping, LC_Calculations):
    """
    classdocs
    """

    def __init__(self, file_location, analyzer='Unknown', instrument_label='Unknown', sample_name=None):
        
        '''
         # Parameters
		----------
		file_location: text,  pathlib.Path(), or s3path.S3Path 
            Path object from pathlib containing the file location
        '''
		
        if  isinstance(file_location, str):
			# if obj is a string it defaults to create a Path obj, pass the S3Path if needed
            file_location = Path(file_location)

        if not file_location.exists():
        
            raise FileExistsError("File does not exist: " + file_location)
        
        self.file_location = file_location
        
        if sample_name: 
        
            self.sample_name = sample_name

        else: 
        
            self.sample_name = file_location.stem
        
        self.analyzer = analyzer
        self.instrument_label = instrument_label
        
        self._retention_time_list = []
        self._scans_number_list = []
        self._tic_list = []

        self._ms = {}
        """
        key is scan number; value is MassSpectrum Class
        """
    
    def __len__(self):
        
        return len(self._ms)
        
    def __getitem__(self, scan_number):
        
        return self._ms.get(scan_number)

    def __iter__(self):

         return iter(self._ms.values()) 

    def add_mass_spectrum(self, mass_spec):

        self._ms[mass_spec.scan_number] = mass_spec

    def set_tic_list_from_data(self):

        self.tic = [self._ms.get(i).tic for i in self.scans_number]
        
        # self.set_tic_list([self._ms.get(i).get_sumed_signal_to_noise() for i in self.get_scans_number()])

    def set_retention_time_from_data(self):

        retention_time_list = []

        for key_ms in sorted(self._ms.keys()):

            retention_time_list.append(self._ms.get(key_ms).rt)

        self.retention_time = retention_time_list 

        # self.set_retention_time_list(sorted(self._ms.keys()))

    def set_scans_number_from_data(self):
        
        self.scans_number = sorted(self._ms.keys())

    @property
    def scans_number(self):

        return self._scans_number_list

    @property
    def retention_time(self):

        return self._retention_time_list
    
    @property
    def tic(self):

        return self._tic_list

    @retention_time.setter
    def retention_time(self, l):
        # self._retention_time_list = linspace(0, 80, num=len(self._scans_number_list))
        self._retention_time_list = l

    @scans_number.setter
    def scans_number(self, l):

        self._scans_number_list = l

    @tic.setter
    def tic(self, l):

        self._tic_list = l    

    