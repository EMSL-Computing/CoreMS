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

from threading import Thread
from typing import List, Dict, Any, Tuple
import numpy as np

class ImportMassSpectraThermoMSFileReader(Thread):
    """
    A class for importing mass spectra from Thermo MS file reader.

    Parameters:
    -----------
    file_location : str
        The file location of the Thermo MS file.

    Attributes:
    -----------
    thermo_Library : object
        The Thermo MS file reader library object.
    res : int
        The result of setting the current controller.
    LCMS : LCMSBase
        The LCMSBase object for storing the imported mass spectra.
    _initial_scan_number : int
        The initial scan number for importing mass spectra.
    _final_scan_number : int
        The final scan number for importing mass spectra.
    file_location : str
        The file location of the Thermo MS file.

    Properties:
    -----------
    initial_scan_number : int
        The initial scan number for importing mass spectra.
    final_scan_number : int
        The final scan number for importing mass spectra.

    Methods:
    --------
    check_scan(scan: int) -> bool:
        Check if the given scan number is valid.
    run() -> None:
        Automatically process mass spectrum in a separate thread.
    get_mass_spectra(auto_process: bool=True) -> LCMSBase:
        Get the imported mass spectra.
    check_load_success() -> bool:
        Check if the Thermo MS file was loaded successfully.
    get_filter_for_scan_num(scan_number: int) -> List[str]:
        Get the filter for the given scan number.
    check_full_scan(scan_number: int) -> bool:
        Check if the scan is a full scan.
    get_polarity_mode(scan_number: int) -> int:
        Get the polarity mode for the given scan number.
    get_data(scan: int, d_parameter: Dict[str, Any]) -> Dict[str, np.ndarray]:
        Get the data for the given scan number.
    get_scans_numbers() -> int:
        Get the total number of scans in the Thermo MS file.
    get_ScanHeaderInfoForScanNum(scan_number: int) -> Tuple[float, float]:
        Get the retention time and TIC for the given scan number.
    is_profile_scan_for_scan_num(scan_number: int) -> bool:
        Check if the scan is a profile scan.
    _import_mass_spectra(d_params: Dict[str, Any], auto_process: bool=True) -> None:
        Import the mass spectra from the Thermo MS file.
    get_lcms() -> LCMSBase:
        Get the LCMSBase object.

    """

    def __init__(self, file_location: str):
        Thread.__init__(self)
        self.thermo_Library = CreateObject("MSFileReader.XRawfile")
        self.thermo_Library.open(file_location)
        self.res: int = self.thermo_Library.SetCurrentController(0, 1)
        self.check_load_success()
        self.LCMS: LCMSBase = LCMSBase(file_location)
        self._initial_scan_number: int = 1
        self._final_scan_number: int = self.get_scans_numbers()
        self.file_location: str = file_location

    @property
    def initial_scan_number(self) -> int:
        return self._initial_scan_number

    @initial_scan_number.setter
    def initial_scan_number(self, start_scan_number: int) -> None:
        if self.check_scan(start_scan_number):
            self._initial_scan_number = start_scan_number
        else:
            raise Exception(
                "startscan and finalscan should be less than %s"
                % self.get_scans_numbers()
            )

    @property
    def final_scan_number(self) -> int:
        return self._final_scan_number

    @final_scan_number.setter
    def final_scan_number(self, final_scan_number: int) -> None:
        if self.check_scan(final_scan_number):
            self._final_scan_number = final_scan_number
        else:
            raise Exception(
                "startscan and finalscan should be less than %s"
                % self.get_scans_numbers()
            )

    def check_scan(self, scan: int) -> bool:
        """
        Check if the given scan number is valid.

        Parameters:
        -----------
        scan : int
            The scan number to check.

        Returns:
        --------
        bool
            True if the scan number is valid, False otherwise.
        """
        scan_numbers = self.get_scans_numbers()
        return scan <= scan_numbers

    def run(self) -> None:
        """
        Automatically process mass spectrum in a separate thread.
        """
        d_parameters = default_parameters(self.file_location)
        self._import_mass_spectra(d_parameters)

    def get_mass_spectra(self, auto_process: bool=True) -> LCMSBase:
        """
        Get the imported mass spectra.

        Parameters:
        -----------
        auto_process : bool, optional
            Whether to automatically process the mass spectra, by default True.

        Returns:
        --------
        LCMSBase
            The LCMSBase object containing the imported mass spectra.
        """
        d_parameters = default_parameters(self.file_location)
        self._import_mass_spectra(d_parameters, auto_process=auto_process)
        return self.LCMS

    def check_load_success(self) -> bool:
        """
        Check if the Thermo MS file was loaded successfully.

        Returns:
        --------
        bool
            True if the Thermo MS file was loaded successfully, False otherwise.
        """
        if self.res == 0:
            self.break_it = False
            return True
        else:
            raise ImportError(str(self.res))

    def get_filter_for_scan_num(self, scan_number: int) -> List[str]:
        """
        Get the filter for the given scan number.

        Parameters:
        -----------
        scan_number : int
            The scan number.

        Returns:
        --------
        List[str]
            The filter for the given scan number.
        """
        str_filter = BSTR(None)
        error = self.thermo_Library.GetFilterForScanNum(scan_number, byref(str_filter))
        if error:
            raise IOError(
                "scan %i GetFilterForScanNum error : %s" % (scan_number, str(error))
            )
        else:
            return str(str_filter.value).split()

    def check_full_scan(self, scan_number: int) -> bool:
        """
        Check if the scan is a full scan.

        Parameters:
        -----------
        scan_number : int
            The scan number.

        Returns:
        --------
        bool
            True if the scan is a full scan, False otherwise.
        """
        scan_mode_symbol = self.get_filter_for_scan_num(scan_number)[4]
        return scan_mode_symbol == "Full"

    def get_polarity_mode(self, scan_number: int) -> int:
        """
        Get the polarity mode for the given scan number.

        Parameters:
        -----------
        scan_number : int
            The scan number.

        Returns:
        --------
        int
            The polarity mode (-1 for negative, 1 for positive).
        """
        polarity_symbol = self.get_filter_for_scan_num(scan_number)[1]
        if polarity_symbol == "+":
            return 1
        elif polarity_symbol == "-":
            return -1
        else:
            raise Exception("Polarity Mode Unknown, please set it manually")

    def get_data(self, scan: int, d_parameter: Dict[str, Any]) -> Dict[str, np.ndarray]:
        """
        Get the data for the given scan number.

        Parameters:
        -----------
        scan : int
            The scan number.
        d_parameter : Dict[str, Any]
            The dictionary of parameters.

        Returns:
        --------
        Dict[str, np.ndarray]
            The data dictionary containing the mass list, abundance list, retention time, and signal-to-noise ratio.
        """
        scan = c_long(scan)
        pvarLabels = VARIANT()
        pvarFlags = VARIANT()
        self.thermo_Library.GetLabelData(pvarLabels, pvarFlags, scan)
        scans_labels = np.array(pvarLabels.value)
        mz = scans_labels[0]
        magnitude = scans_labels[1]
        rp = scans_labels[2]
        base_noise = scans_labels[3]
        noise = scans_labels[4]
        array_noise_std = (np.array(noise) - np.array(base_noise)) / 3
        l_signal_to_noise = np.array(magnitude) / array_noise_std
        d_parameter["baseline_noise"] = np.average(array_noise_std)
        d_parameter["baseline_noise_std"] = np.std(array_noise_std)
        data_dict = {
            Labels.mz: mz,
            Labels.abundance: magnitude,
            Labels.rp: rp,
            Labels.s2n: l_signal_to_noise,
        }
        return data_dict

    def get_scans_numbers(self) -> int:
        """
        Get the total number of scans in the Thermo MS file.

        Returns:
        --------
        int
            The total number of scans.
        """
        nScans = c_long()
        self.thermo_Library.GetNumSpectra(nScans)
        return int(nScans.value)

    def get_ScanHeaderInfoForScanNum(self, scan_number: int) -> Tuple[float, float]:
        """
        Get the retention time and TIC for the given scan number.

        Parameters:
        -----------
        scan_number : int
            The scan number.

        Returns:
        --------
        Tuple[float, float]
            The retention time and TIC.
        """
        nScanNumber = c_long(scan_number)
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

    def is_profile_scan_for_scan_num(self, scan_number: int) -> bool:
        """
        Check if the scan is a profile scan.

        Parameters:
        -----------
        scan_number : int
            The scan number.

        Returns:
        --------
        bool
            True if the scan is a profile scan, False otherwise.
        """
        IsProfileScan = c_long()
        error = self.thermo_Library.IsProfileScanForScanNum(
            c_long(scan_number), byref(IsProfileScan)
        )
        if error:
            raise IOError("IsProfileScanForScanNum error :", error)
        # print (IsProfileScan.value, bool(1))
        return bool(IsProfileScan.value)

    def _import_mass_spectra(self, d_params: Dict[str, Any], auto_process: bool=True) -> None:
        """
        Import the mass spectra from the Thermo MS file.

        Parameters:
        -----------
        d_params : Dict[str, Any]
            The dictionary of parameters.
        auto_process : bool, optional
            Whether to automatically process the mass spectra, by default True.
        """
        results = []
        if self.check_load_success():
            list_Tics = list()
            list_RetentionTimeSeconds = list()
            list_scans = list()
            for scan_number in range(
                self.initial_scan_number, self.final_scan_number + 1
            ):
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
                        mass_spec = MassSpecCentroid(data_dict, d_params)
                        self.LCMS.add_mass_spectrum(mass_spec)
            self.LCMS.retention_time = list_RetentionTimeSeconds
            self.LCMS.set_tic_list(list_Tics)
            self.LCMS.set_scans_number_list(list_scans)

    def get_lcms(self) -> LCMSBase:
        """
        Get the LCMSBase object.

        Returns:
        --------
        LCMSBase
            The LCMSBase object.
        """
        if self.LCMS.get(self._initial_scan_number):
            return self.LCMS
        else:
            raise Exception("returning an empty LCMS class")
