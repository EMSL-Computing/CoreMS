
import sys
sys.path.append(".")

from enviroms.emsl.yec.structure.factory.mass_spec.MS_Profile import MassSpecProfile
from enviroms.emsl.yec.structure.factory.LC_Class import LCMSBase

from collections import namedtuple
from ctypes import c_double, c_long
from threading import Thread

import numpy
import pandas
from comtypes import byref
from comtypes.automation import BSTR, VARIANT
from comtypes.client import CreateObject
from pandas import DataFrame




__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"


class ReadAllScansThermoRawData():
    ''' Read Raw File'''

    def __init__(self, file_location):

        # Thread.__init__(self)

        self.file_location = file_location

        self.lCMSBase = LCMSBase(file_location)

        self.break_it = False

    def d_parms(self):

        output_parameters = dict()

        output_parameters["Aterm"] = None

        output_parameters["Bterm"] = None

        output_parameters["Cterm"] = None

        output_parameters["exc_high_freq"] = None

        output_parameters["exc_low_freq"] = None

        output_parameters["bandwidth"] = None

        output_parameters["number_data_points"] = 0

        output_parameters["polarity"] = None

        output_parameters["filename_path"] = self.file_location

        """scan_number and rt will be need to lc ms"""

        output_parameters["mobility_scan"] = 0

        output_parameters["mobility_rt"] = 0

        output_parameters["scan_number"] = 0

        output_parameters["rt"] = 0

        return output_parameters

    def run(self):

        self.thermo_Library = CreateObject('MSFileReader.XRawfile')

        self.thermo_Library.open(self.file_location)

        self.res = self.thermo_Library.SetCurrentController(0, 1)

        d_parms = self.d_parms()

        if self.check_load_sucess():

            self.get_mass_spectrums(self.file_location, d_parms)

        else:

            self.break_it = True

        del self.thermo_Library

        return self.lCMSBase

    def check_load_sucess(self):
        ''' 0 if successful; otherwise, see Error Codes '''
        if self.res == 0:

            return True

        else:

            return False

    def get_filter_for_scan_num(self, scan_number):
        """Returns the closest matching run time that corresponds to scan_number for the current
        controller. This function is only supported for MS device controllers.
        e.g.  ['FTMS', '-', 'p', 'NSI', 'Full', 'ms', '[200.00-1000.00]']
        """
        pbstrFilter = BSTR(None)
        error = self.thermo_Library.GetFilterForScanNum(
            scan_number, byref(pbstrFilter))
        if error:
            raise IOError("scan {}, GetFilterForScanNum error : {}".format(
                scan_number, str(error)))
        else:
            return str(pbstrFilter.value).split()

    def check_full_scan(self, scan_number):

        scan_mode_symbol = self.get_filter_for_scan_num(scan_number)[4]

        return scan_mode_symbol == 'Full'

    def get_polarity_mode(self, scan_number):

        polarity_symbol = self.get_filter_for_scan_num(scan_number)[1]

        if polarity_symbol == "+":

            return 1
            # return "POSITIVE_ION_MODE"

        elif polarity_symbol == "-":

            return -1

        else:

            raise Exception("Polarity Mode Unknown, please set it manually")

    def get_data(self, scan):

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
        #charge = scans_labels[5]

        array_noise_std = (numpy.array(noise) - numpy.array(base_noise)) / 3
        l_signal_to_noise = numpy.array(magnitude)/array_noise_std

        data_dict = {"m/z": mz, "Abundance": magnitude,
                     "Resolving Power": rp, "S/N": l_signal_to_noise}

        return data_dict

    def get_scans_numbers(self):

        nScans = c_long()
        self.thermo_Library.GetNumSpectra(nScans)

        return nScans.value

    def get_initial_rt(self):

        initial_rt = c_double()
        self.thermo_Library.GetStartTime(initial_rt)

        return initial_rt.value

    def get_final_rt(self):

        final_rt = c_double()
        self.thermo_Library.GetEndTime(final_rt)

        return final_rt.value

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
        self.thermo_Library.GetScanHeaderInfoForScanNum(nScanNumber,
                                                        nPackets, dRetantionTime,
                                                        dLowMass,
                                                        dHighMass, dTIC,
                                                        dBasePeakMass,
                                                        dBasePeakIntensity,
                                                        nChannels,
                                                        bUniformTime,
                                                        dFrequency)

        return dRetantionTime.value, dTIC.value

    def IsProfileScanForScanNum(self, scan_number):
        IsProfileScan = c_long()
        error = self.thermo_Library.IsProfileScanForScanNum(
            c_long(scan_number), byref(IsProfileScan))
        if error:
            raise IOError("IsProfileScanForScanNum error :", error)
        return bool(IsProfileScan.value)

    def get_mass_spectrums(self, file_name, d_parms):

        # Each_Mass_Spectrum = namedtuple('each_mass_spectrum', ['mass_list', 'abundance_list', 'retention_time', 'scan_number', 'tic_number'])

        if self.check_load_sucess():

            '''get number of scans'''

            numero_scans = int(self.get_scans_numbers())

            if numero_scans > 0:

                list_Tics = list()

                list_RetentionTimeSeconds = list()

                list_scans = list()

                '''key = scan_number or retention time'''

                for scan_number in range(numero_scans):

                    scan_number = scan_number + 1

                    "only import FULL scans"

                    if self.check_full_scan(scan_number):

                        d_parms["polarity"] = self.get_polarity_mode(
                            scan_number)

                        d_parms["rt"], TIC = self.get_ScanHeaderInfoForScanNum(
                            scan_number)

                        d_parms["scan_number"] = scan_number

                        list_RetentionTimeSeconds.append(d_parms.get("rt"))

                        list_Tics.append(TIC)

                        list_scans.append(scan_number)

                        # list_masses, list_Intensities = self.get_segment_mass_list_from_scan_num(scan)
                        data_dict = self.get_data(scan_number)
                        #data = DataFrame(data_dict)
                        #mass_spec = MassSpecProfile(data, d_parms)
                        # mass_spec.process_mass_spec()
                        '''find a way to get polarity from raw file'''

                        # self.lCMSBase.set_mass_spectrum_for_scan(mass_spec)

                # return each_mass_spectrum, list_Tics, array(list_RetentionTimeSeconds)
                self.lCMSBase.set_retention_time_list(
                    list_RetentionTimeSeconds)
                self.lCMSBase.set_tic_list(list_Tics)
                self.lCMSBase.set_scans_number_list(list_scans)

            else:

                self.break_it = True
                # return None, None, None

        else:
            pass
            # s elf.break_it = True
            # return None, None, None


a = ReadAllScansThermoRawData(
    "C:\\Users\\eber373\\Desktop\\WK_ps_lignin_190301112616.raw")
a.run()
