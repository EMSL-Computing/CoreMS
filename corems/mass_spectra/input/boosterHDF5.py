__author__ = "Yuri E. Corilo"
__date__ = "Oct 29, 2019"

from pandas import DataFrame
import h5py
from threading import Thread

from corems.encapsulation.constant import Labels
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecProfile
from corems.mass_spectra.factory.LC_Class import LCMSBase
from corems.encapsulation.settings.input import InputSetting

class ReadHDF_BoosterMassSpectra(Thread):
    
    '''class docs'''
    
    def __init__(self, file_location, auto_process=True):

        Thread.__init__(self)

        self.lcms = LCMSBase(file_location)

        self.hdf_obj =  h5py.File(file_location, 'r')

        self.list_scans =  sorted([int(i) for i in list(self.hdf_obj.keys())])

        self.initial_scan_number = self.list_scans [0]

        self.final_scan_number = self.list_scans [-1]

        self.file_location = file_location

        self.auto_process = True
    
    def get_polarity(self, file_location, scan):

        self.h5pydata = h5py.File(file_location, 'r')

        self.scans = list(self.h5pydata.keys())
        
        polarity = self.get_attr_data(scan,'r_h_polarity')
        
        if polarity == 'negative scan':
            
            return -1
        else:
            return +1    
    
    def get_attr_data(self, scan, attr_srt):

        return self.hdf_obj[str(scan)].attrs[attr_srt]

    def import_mass_spectra(self, d_parms):
        
        list_rt, list_tic = list(), list()
        
        for scan_number in self.list_scans:
            
            d_parms["rt"] =  list_rt.append(self.get_attr_data(scan_number, 'r_h_start_time'))

            d_parms["scan_number"] = scan_number

            d_parms['label'] = Labels.booster_profile
    
            d_parms["polarity"] = self.get_polarity(self.file_location, scan_number)

            d_parms["Aterm"] = self.get_attr_data(scan_number, 'r_cparams')[0]

            d_parms["Bterm"] = self.get_attr_data(scan_number, 'r_cparams')[1]

            list_rt.append(d_parms["rt"])

            list_tic.append(self.get_attr_data(scan_number, 'r_h_tic'))
            
            mass_spec = self.get_mass_spectrum(scan_number, d_parms)

            self.lcms.add_mass_spectrum_for_scan(mass_spec)

        self.lcms.set_retention_time_list(list_rt)
        self.lcms.set_tic_list(list_tic)
        self.lcms.set_scans_number_list(self.list_scans)
        
    def get_mass_spectrum(self, scan, d_parms):
        
        booster_data = self.hdf_obj[str(scan)]
        
        if booster_data.shape[0] is not 2:
            
            raise NotImplementedError('opening transient, needs read raw file here, get bandwidth, create transient class and then the mass spectrum')
       
        else:
            
            data_dict = {
                "m/z": booster_data[0],
                "Abundance": booster_data[1],
                "Resolving Power": None,
                "S/N": None,
            }
            
            data = DataFrame(data_dict)
            mass_spec = MassSpecProfile(data, d_parms, auto_process=self.auto_process)

        return mass_spec

    def run(self):
        '''creates the lcms obj'''

        d_parameters = InputSetting.d_parms(self.file_location)
        self.import_mass_spectra(d_parameters)
            
    def get_lcms_obj(self):
        
        if self.lcms.get_mass_spec_by_scan_number(self.initial_scan_number):
            return self.lcms
        else:
            raise Exception("returning a empty lcms class")
    