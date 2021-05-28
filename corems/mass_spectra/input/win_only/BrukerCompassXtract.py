from ctypes import c_long
from threading import Thread

from comtypes.automation import BSTR
from comtypes.client import CreateObject
from numpy import array


from corems.encapsulation.factory.parameters import default_parameters
from corems.encapsulation.constant import Labels
from corems.mass_spectra.factory.LC_Class import LCMSBase
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecProfile

__author__ = "Yuri E. Corilo"
__date__ = "July 10, 2019"

class ImportLCMSBrukerCompassXtract(Thread):
    '''class docs'''
    
    def __init__(self, file_location, auto_process=True):

        Thread.__init__(self)

        self.lcms = LCMSBase(file_location)

        """Set up the COM object interface"""
        self.Bruker_Library = CreateObject("EDAL.MSAnalysis")

        self.res = self.Bruker_Library.Open(file_location)

        self.check_load_sucess()

        self._initial_scan_number = 1

        self._final_scan_number = self.get_scans_numbers()

        self.file_location = file_location

        self.auto_process = auto_process

    @property
    def initial_scan_number(self):
        return self._initial_scan_number

    @property
    def final_scan_number(self):
        return self._final_scan_number

    def check_scan(self, scan):

        scan_numbers = self.get_scans_numbers()
        return scan <= scan_numbers

    @initial_scan_number.setter
    def initial_scan_number(self, initial_scan_number):
        if self.check_scan(initial_scan_number):
            self._initial_scan_number = initial_scan_number
        else:
            raise Exception(
                "startscan and finalscan should be less than %s"
                % self.get_scans_numbers()
            )

    @final_scan_number.setter
    def final_scan_number(self, final_scan_number):

        if self.check_scan(final_scan_number):
            self._final_scan_number = final_scan_number
        else:
            raise Exception(
                "startscan and finalscan should be less than %s"
                % self.get_scans_numbers()
            )

    def get_scans_numbers(self):

        scan_numbers = self.Bruker_Library.MSSpectrumCollection.Count

        return scan_numbers

    def get_polarity_mode(self, spectrum):

        polarity_symbol = spectrum.Polarity

        if polarity_symbol == 0:

            return 1
            # return "POSITIVE_ION_MODE"

        elif polarity_symbol == 1:

            return -1
            # return "NEGATIVE_ION_MODE"

        else:

            raise IOError("Could not read mass spectrum polarity mode")

    def check_load_sucess(self):
        """ 0 if successful; otherwise, see Error Codes """

        if self.res == 0:

            self.break_it = False

        else:

            raise ImportError(str(self.res))

    def get_bruker_tics(self):

        strAnalysisData = BSTR("SumIntensity")

        if self.Bruker_Library.HasAnalysisData(strAnalysisData):

            tics_array = self.Bruker_Library.GetAnalysisData(strAnalysisData)

            tics_array = array(tics_array)

        return tics_array

    def get_bruker_retention_time(self):

        strAnalysisData = BSTR("RetentionTime")

        if self.Bruker_Library.HasAnalysisData(strAnalysisData):

            tics_array = self.Bruker_Library.GetAnalysisData(strAnalysisData)

            tics_array = array(tics_array)
        else:
            tics_array = [0]

        return tics_array

    @staticmethod
    def get_data(spectra, scan):
        """init_variable_from_get_spectrums
        # massList set up later
        #retention_time = spectrum.RetentionTime
        """

        spectrum = spectra[scan]

        is_profile = c_long(1)

        masslist = spectrum.GetMassIntensityValues(is_profile)

        # index_to_cut = self.find_index_of_mass(1200, masslist[0])

        data_dict = {
            Labels.mz: array(masslist[0]),
            Labels.abundance: array(masslist[1]),
            Labels.rp: None,
            Labels.s2n: None,
        }

        return data_dict

    def run(self):
        '''creates the lcms obj'''
        d_parameters = default_parameters(self.file_location)
        self._import_mass_spectra(d_parameters)

    def _import_mass_spectra(self, d_params):

        spectra = self.Bruker_Library.MSSpectrumCollection

        list_rt = self.get_bruker_retention_time()

        list_Tics = self.get_bruker_tics()

        list_scans = list()

        for scan_number in range(self.initial_scan_number, self.final_scan_number + 1):

            if spectra[scan_number].MSMSStage == 1:
                # this label needs to go inside a encapsulation class for consistence
                d_params["label"] = Labels.bruker_profile

                d_params["polarity"] = self.get_polarity_mode(spectra[scan_number])

                d_params["rt"] = list_rt[scan_number - 1]

                d_params["scan_number"] = scan_number

                list_scans.append(scan_number)

                data_dict = self.get_data(spectra, scan_number)

                mass_spec = MassSpecProfile(data_dict, d_params, auto_process=self.auto_process)
                
                mass_spec.process_mass_spec()
                
                self.lcms.add_mass_spectrum(mass_spec)

        self.lcms.retention_time = list_rt
        self.lcms.tic = list_Tics
        self.lcms.scans_number = list_scans 
        # return each_mass_spectrum

    def get_lcms_obj(self):
        """get_lc_ms_class method should only be used when using this class as a Thread, 
        otherwise use the run() method to return the lcms class"""

        if self.lcms.get(self._initial_scan_number):
            return self.lcms
        else:
            raise Exception("returning a empty lcms class")
