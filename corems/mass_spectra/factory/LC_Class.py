"""
Created on Jun 12, 2019
"""
from corems.mass_spectra.calc.LC_Calc import LC_Calculations

__author__ = "Yuri E. Corilo"
__date__ = "Jun 25, 2019"

class LCMSBase(LC_Calculations):
    """
    classdocs
    """

    def __init__(self, file_location, analyzer='Unknown', instrument_label='Unknown', sample_name=None):
        
        """
        Constructor
        file_location: pathlib.Path() 
            Path object from pathlib containing the file location
        """
        
        if not file_location.exists():
        
            raise FileExistsError("File does not exist: " + file_location)
        
        self.file_location = file_location
        
        if sample_name: 
        
            self.sample_name = sample_name

        else: 
        
            self.sample_name = file_location.stem
        
        self.analyzer = analyzer
        self.instrument_label = instrument_label
        
        self.retention_time_list = []
        self.scans_number_list = []
        self.tic_list = []

        self.ms = {}
        """
        key is scan number; value is MassSpectrum Class
        """

    def __getitem__(self, i):
        
        scan_number = self.scans_number_list[i] 
        return self.ms[scan_number]
    
    def __iter__(self):
        self.cur_scan = 0
        return self 
        
    def __next__(self):
        
        i = self.cur_scan
        if i >= len(self.scans_number_list):
            raise StopIteration
        self.cur_scan += 1
        next_scan_number = self.scans_number_list[i] 
        return self.ms[next_scan_number]
        
    def add_mass_spectrum_for_scan(self, mass_spec):

        self.ms[mass_spec.scan_number] = mass_spec

    def get_mass_spec_by_scan_number(self, scan):

        return self.ms.get(scan)

    def set_tic_list_from_data(self):

        self.set_tic_list(
            [self.ms.get(i).tic for i in self.get_scans_number()]
        )

        # self.set_tic_list([self.ms.get(i).get_sumed_signal_to_noise() for i in self.get_scans_number()])

    def set_retention_time_from_data(self):

        retention_time_list = []

        for key_ms in sorted(self.ms.keys()):

            retention_time_list.append(self.ms.get(key_ms).rt)

        self.set_retention_time_list(retention_time_list)

        # self.set_retention_time_list(sorted(self.ms.keys()))

    def set_scans_number_from_data(self):
        self.set_scans_number_list(sorted(self.ms.keys()))

    def get_scans_number(self):

        return self.scans_number_list

    def get_retention_time(self):

        return self.retention_time_list

    def get_tic(self):

        return self.tic_list

    def set_retention_time_list(self, lista):
        # self.retention_time_list = linspace(0, 80, num=len(self.scans_number_list))
        self.retention_time_list = lista

    def set_scans_number_list(self, lista):

        self.scans_number_list = lista

    def set_tic_list(self, lista):

        self.tic_list = lista

    