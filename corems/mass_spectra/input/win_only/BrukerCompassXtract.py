from ctypes import c_long
from threading import Thread

from comtypes.automation import BSTR
from comtypes.client import CreateObject
from numpy import array
from pandas import DataFrame

from corems.encapsulation.settings.input import InputParameters
from corems.mass_spectra.factory.LC_Class import LCMSBase
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecProfile

__author__ = "Yuri E. Corilo"
__date__ = "July 10, 2019"

class ImportLCMSBrukerCompassXtract(Thread):
    '''class docs'''
    
    def __init__(self, file_location):

        Thread.__init__(self)

        """Set up the COM object interface"""
        self.Bruker_Library = CreateObject("EDAL.MSAnalysis")

        self.res = self.Bruker_Library.Open(file_location)

        self.check_load_sucess()

        self.LCMS = LCMSBase(file_location)

        self._initial_scan_number = 1

        self._final_scan_number = self.get_scans_numbers()

        self.file_location = file_location

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
            "m/z": array(masslist[0]),
            "Abundance": array(masslist[1]),
            "Resolving Power": None,
            "S/N": None,
        }

        return data_dict

    def run(self):
        '''thread will automatically process mass spectrum
        use the get_mass_spectra class to import without processing mass spectrum
        use get_lcms to acess resulting class
        '''

        d_parameters = InputParameters.d_parms(self.file_location)
        self._import_mass_spectra(d_parameters)

        # return LCMS
    
    def get_mass_spectra(self,auto_process=True):

        d_parameters = InputParameters.d_parms(self.file_location)
        self._import_mass_spectra(d_parameters, auto_process=auto_process)
        return self.LCMS

    def _import_mass_spectra(self, d_parms, auto_process=True):

        spectra = self.Bruker_Library.MSSpectrumCollection

        list_RetentionTimeSeconds = self.get_bruker_retention_time()

        list_Tics = self.get_bruker_tics()

        list_scans = list()

        for scan_number in range(self.initial_scan_number, self.final_scan_number + 1):

            if spectra[scan_number].MSMSStage == 1:
                # this label needs to go inside a encapsulation class for consistence
                d_parms["label"] = "Bruker_Profile"

                d_parms["polarity"] = self.get_polarity_mode(spectra[scan_number])

                d_parms["rt"] = list_RetentionTimeSeconds[scan_number - 1]

                d_parms["scan_number"] = scan_number

                list_scans.append(scan_number)

                data_dict = self.get_data(spectra, scan_number)

                data = DataFrame(data_dict)
                mass_spec = MassSpecProfile(data, d_parms, auto_process=auto_process)
                mass_spec.process_mass_spec()
                self.LCMS.add_mass_spectrum_for_scan(mass_spec)

        self.LCMS.set_retention_time_list(list_RetentionTimeSeconds)
        self.LCMS.set_tic_list(list_Tics)
        self.LCMS.set_scans_number_list(list_scans)
        # return each_mass_spectrum

    def get_lcms(self):
        """get_lc_ms_class method should only be used when using this class as a Thread, 
        otherwise use the run() method to return the LCMS class"""

        if self.LCMS.get_mass_spec_by_scan_number(self._initial_scan_number):
            return self.LCMS
        else:
            raise Exception("returning a empty LCMS class")
