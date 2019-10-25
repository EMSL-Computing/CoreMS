

from corems.encapsulation.settings.input import InputParameters
from corems.encapsulation.constant import Labels

from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecProfile
from corems.mass_spectra.factory.LC_Class import LCMSBase

__author__ = "Yuri E. Corilo"
__date__ = "July 10, 2019"

class ImportBoosterMassSpectra(Thread):
    '''class docs'''
    
    def __init__(self, file_location, polarity, auto_process=True):

        Thread.__init__(self)

        self.lcms = LCMSBase(file_location)

        self.hdf_obj =  h5py.File(file_location, 'r')

        self._initial_scan_number = 0

        self._final_scan_number = self.get_scans_numbers()

        self.file_location = file_location

        self.polarity = polarity

        self.auto_process = True

     def get_scans_numbers(self):

        scan_numbers = len(list(self.hfd_obj.keys()))

        return scan_numbers
     
     def _import_mass_spectra(self, d_parms):

        spectra = self.Bruker_Library.MSSpectrumCollection

        list_scans, list_rt = list(), list()
        
        for scan_number in range(self.initial_scan_number, self.final_scan_number + 1):
            
                d_parms["rt"] = float(list(h5pydata.keys())[scan_number])

                d_parms["scan_number"] = scan_number

                d_parms['label'] = Labels.simulated_profile
        
                d_parms["polarity"] = self.polarity

                list_rt.append(d_parms["rt"])
                
                list_scans.append(scan_number)
                
                data_dict = self.get_data(scan_number)

                data = DataFrame(data_dict)
                mass_spec = MassSpecProfile(data, d_parms, auto_process=auto_process)
                mass_spec.process_mass_spec()
                self.LCMS.add_mass_spectrum_for_scan(mass_spec)

        self.lcms.set_retention_time_list(list_rt)
        self.lcms.set_tic_list_from_data()
        self.lcms.set_scans_number_list(list_scans)
        # return each_mass_spectrum
     
    def get_data(scan):
        
        scans = list(h5pydata.keys())

        booster_data = self.hdf_obj[scans[scan]]

        # index_to_cut = self.find_index_of_mass(1200, masslist[0])

        data_dict = {
            "m/z": booster_data[0],
            "Abundance": booster_data[1],
            "Resolving Power": None,
            "S/N": None,
        }

        return data_dict

    def run(self):
        '''creates the lcms obj'''

        d_parms = InputParameters.d_parms(self.file_location)
        self._import_mass_spectra(d_parameters)
            
    def get_lcms_obj(self):
        
        if self.lcms.get_mass_spec_by_scan_number(self._initial_scan_number):
            return self.lcms
        else:
            raise Exception("returning a empty lcms class")
    