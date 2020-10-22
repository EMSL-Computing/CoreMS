from comtypes import byref
from comtypes.automation import BSTR, VARIANT
from comtypes.client import CreateObject
from ctypes import c_double, c_long

from corems.encapsulation.factory.parameters import default_parameters
from corems.mass_spectra.factory.LC_Class import LCMSBase
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecProfile, MassSpecCentroid
from corems.encapsulation.constant import Labels

from threading import Thread
import multiprocessing
import numpy

__author__ = "Yuri E. Corilo"
__date__ = "July 9, 2019"

class ImportMassSpectraThermoMSFileReader(Thread):
    
    """     Read FULL and PROFILE (it ignores all other scans) raw file data and store it return a LCMS class
    *  Default behavior is to load all scans numbers

    *  set start_scan_number  and final_scan_number to change it before calling start(), or run()

    *  Noise threshold will break the mass_spec.process_mass_spec() if the method in the
    MassSpecSetting class is set to something other than Relative Abundance
    (it needs to be fixed to work with all methods)
    """

    def __init__(self, file_location):

        Thread.__init__(self)

        self.thermo_Library = CreateObject("MSFileReader.XRawfile")

        self.thermo_Library.open(file_location)

        self.res = self.thermo_Library.SetCurrentController(0, 1)

        self.check_load_success()

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
    def initial_scan_number(self, start_scan_number):
        if self.check_scan(start_scan_number):
            self._initial_scan_number = start_scan_number
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

    
    def run(self):
        '''thread will automatically process mass spectrum
        use the get_mass_spectra class to import without processing mass spectrum'''

        d_parameters = default_parameters(self.file_location)
        self._import_mass_spectra(d_parameters)

        # return self.LCMS

    def get_mass_spectra(self,auto_process=True):

        d_parameters = default_parameters(self.file_location)
        self._import_mass_spectra(d_parameters, auto_process=auto_process)
        return self.LCMS

    def check_load_success(self):
        """ 0 if successful; otherwise, see Error Codes on MSFileReader Manual """
        if self.res == 0:

            self.break_it = False
            return True
        else:

            raise ImportError(str(self.res))

    def get_filter_for_scan_num(self, scan_number):
        """Returns the closest matching run time that corresponds to scan_number for the current
        controller. This function is only supported for MS device controllers.
        e.g.  ['FTMS', '-', 'p', 'NSI', 'Full', 'ms', '[200.00-1000.00]']
        """
        str_filter = BSTR(None)
        error = self.thermo_Library.GetFilterForScanNum(scan_number, byref(str_filter))
        if error:
            raise IOError(
                "scan %i GetFilterForScanNum error : %s" % (scan_number, str(error))
            )
        else:
            return str(str_filter.value).split()

    def check_full_scan(self, scan_number):

        scan_mode_symbol = self.get_filter_for_scan_num(scan_number)[4]

        return scan_mode_symbol == "Full"

    def get_polarity_mode(self, scan_number):

        polarity_symbol = self.get_filter_for_scan_num(scan_number)[1]

        if polarity_symbol == "+":

            return 1
            # return "POSITIVE_ION_MODE"

        elif polarity_symbol == "-":

            return -1

        else:

            raise Exception("Polarity Mode Unknown, please set it manually")

    def get_data(self, scan, d_parameter):

        scan = c_long(scan)
        pvarLabels = VARIANT()
        pvarFlags = VARIANT()

        self.thermo_Library.GetLabelData(pvarLabels, pvarFlags, scan)
        scans_labels = numpy.array(pvarLabels.value)

        mz = scans_labels[0]
        magnitude = scans_labels[1]
        rp = scans_labels[2]
        base_noise = scans_labels[3]
        noise = scans_labels[4]
        # charge = scans_labels[5]

        array_noise_std = (numpy.array(noise) - numpy.array(base_noise)) / 3
        l_signal_to_noise = numpy.array(magnitude) / array_noise_std

        d_parameter["baselise_noise"] = numpy.average(array_noise_std)

        d_parameter["baselise_noise_std"] = numpy.std(array_noise_std)

        data_dict = {
            Labels.mz: mz,
            Labels.abundance: magnitude,
            Labels.rp: rp,
            Labels.s2n: l_signal_to_noise,
        }

        return data_dict

    def get_scans_numbers(self):

        nScans = c_long()
        self.thermo_Library.GetNumSpectra(nScans)

        return int(nScans.value)

    def get_ScanHeaderInfoForScanNum(self, scan_number):

        nScanNumber = c_long(scan_number)  # get info for the twelfth scan
        nPackets = c_long(0)
        dRetantionTime = c_double(0.0)
        dLowMass = c_double(0.0)
        dHighMass = c_double(0.0)
        dTIC = c_double(0.0)
        dBasePeakMass = c_double(0.0)
        dBasePeakIntensity = c_double(0.0)
        nChannels = c_long(0)
        bUniformTime = c_long(False)
        dFrequency = c_double(0.0)
        self.thermo_Library.GetScanHeaderInfoForScanNum(
            nScanNumber,
            nPackets,
            dRetantionTime,
            dLowMass,
            dHighMass,
            dTIC,
            dBasePeakMass,
            dBasePeakIntensity,
            nChannels,
            bUniformTime,
            dFrequency,
        )

        return dRetantionTime.value, dTIC.value

    def is_profile_scan_for_scan_num(self, scan_number):

        IsProfileScan = c_long()
        error = self.thermo_Library.IsProfileScanForScanNum(
            c_long(scan_number), byref(IsProfileScan)
        )
        if error:
            raise IOError("IsProfileScanForScanNum error :", error)
        # print (IsProfileScan.value, bool(1))
        return bool(IsProfileScan.value)

    def _import_mass_spectra(self, d_params, auto_process=True):
        results = []
        # Each_Mass_Spectrum = namedtuple('each_mass_spectrum', ['mass_list', 'abundance_list', 'retention_time', 'scan_number', 'tic_number'])

        if self.check_load_success():

            """get number of scans"""

            list_Tics = list()

            list_RetentionTimeSeconds = list()

            list_scans = list()

            """key = scan_number or retention time"""
            # print(self.initial_scan_number, self.final_scan_number)
            for scan_number in range(
                self.initial_scan_number, self.final_scan_number + 1
            ):
                #print(scan_number)
                # scan_number = scan_number + 1

                "only import FULL scans and Profile Mode, it ignores all others"

                if self.check_full_scan(scan_number):

                    if self.is_profile_scan_for_scan_num(scan_number):

                        d_params["label"] = Labels.thermo_centroid

                        d_params["polarity"] = self.get_polarity_mode(scan_number)

                        d_params["rt"], TIC = self.get_ScanHeaderInfoForScanNum(
                            scan_number
                        )

                        d_params["scan_number"] = scan_number

                        list_RetentionTimeSeconds.append(d_params.get("rt"))

                        list_Tics.append(TIC)

                        list_scans.append(scan_number)

                        data_dict = self.get_data(scan_number, d_params)

                        #results.append((data, d_params))
                        
                        mass_spec = MassSpecCentroid(data_dict, d_params)
                        
                        self.LCMS.add_mass_spectrum(mass_spec)

            #pool = multiprocessing.Pool(5)
            #result = pool.starmap(MassSpecCentroid, results)
            #for ms in result:
            #self.LCMS.add_mass_spectrum(ms)
            
            self.LCMS.retention_time = list_RetentionTimeSeconds
            self.LCMS.set_tic_list(list_Tics)
            self.LCMS.set_scans_number_list(list_scans)

    def get_lcms(self):
        """get_lc_ms_class method should only be used when using this class as a Thread, 
        otherwise use the run() method to return the LCMS class"""

        if self.LCMS.get(self._initial_scan_number):
            return self.LCMS
        else:
            raise Exception("returning a empty LCMS class")
