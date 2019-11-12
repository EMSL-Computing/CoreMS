__author__ = "Yuri E. Corilo"
__date__ = "Oct 29, 2019"

from pandas import DataFrame
import h5py
from threading import Thread

from corems.encapsulation.constant import Labels
from corems.mass_spectra.factory.LC_Class import LCMSBase
from corems.encapsulation.settings.input import InputSetting
from corems.transient.input.BrukerSolarix import ReadBrukerSolarix

class ReadBruker_SolarixTransientMassSpectra(Thread):
    
    '''class docs'''
    
    def __init__(self, d_directory_location, analyzer='ICR', instrument_label='15T', 
                       auto_process=True, auto_noise=True, keep_profile=False):

        Thread.__init__(self)

        if not d_directory_location.exists():
            raise FileNotFoundError("File does not exist: " + str(d_directory_location))
        
        self.scan_attr = d_directory_location / "scan.xml"
        if not self.scan_attr.exists():
            raise FileExistsError("%s does not seem to be a valid Solarix Mass Spectra Experiment,\
                                maybe an Imaging experiment?\
                                please ReadBruker_SolarixTransientImage class for Imaging dataset " % d_directory_location)

        self.lcms = LCMSBase(d_directory_location, analyzer, instrument_label)

        self.auto_process = auto_process
        self.auto_noise = auto_noise
        self.keep_profile = keep_profile
    
    def get_scan_attr(self):
    
        from bs4 import BeautifulSoup
        
        soup = BeautifulSoup(self.scan_attr.open(),'xml')

        list_rt = [float(rt.text) for rt in soup.find_all('minutes')]
        list_tic = [float(tic.text) for tic in soup.find_all('tic')]
        list_scan = [int(scan.text) for scan in soup.find_all('count')]
    
        return list_scan, list_rt, list_tic
    
    def import_mass_spectra(self):
        
        list_scan, list_rt, list_tic = self.get_scan_attr()
        
        for scan_index, _ in enumerate(list_scan):
            
            mass_spec = self.get_mass_spectrum(scan_index)

            self.lcms.add_mass_spectrum_for_scan(mass_spec)

        self.lcms.set_retention_time_list(list_rt)
        self.lcms.set_tic_list(list_tic)
        self.lcms.set_scans_number_list(list_scan)
        
    def get_mass_spectrum(self, scan_index):
        
        bruker_reader = ReadBrukerSolarix(self.lcms.file_location)

        bruker_transient = bruker_reader.get_transient(scan_index=scan_index)

        mass_spec = bruker_transient.get_mass_spectrum(plot_result=False, 
                                                       auto_process=self.auto_process,
                                                       keep_profile=self.keep_profile, 
                                                       auto_noise=self.auto_noise)

        return mass_spec

    def run(self):
        '''creates the lcms obj'''
        self.import_mass_spectra()
            
    def get_lcms_obj(self):
        
        if self.lcms:
            return self.lcms
        else:
            raise Exception("returning a empty lcms class")
    