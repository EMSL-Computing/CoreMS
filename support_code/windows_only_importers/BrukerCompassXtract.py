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
    """Class for importing LC-MS data from Bruker Compass Xtract files.

    Parameters:
    -----------
    file_location : str
        The path to the Bruker Compass Xtract file.
    auto_process : bool, optional
        Flag indicating whether to automatically process the imported data.
        Default is True.

    Attributes:
    -----------
    initial_scan_number : int
        The initial scan number to import.
    final_scan_number : int
        The final scan number to import.

    Methods:
    --------
    check_scan(scan: int) -> bool:
        Check if the given scan number is valid.
    get_scans_numbers() -> int:
        Get the total number of scans in the file.
    get_polarity_mode(spectrum) -> int:
        Get the polarity mode of a given spectrum.
    check_load_success() -> None:
        Check if the file was loaded successfully.
    get_bruker_tics() -> array:
        Get the total ion current (TIC) array from the file.
    get_bruker_retention_time() -> array:
        Get the retention time array from the file.
    get_data(spectra, scan: int) -> dict:
        Get the mass spectrum data for a given scan number.
    run() -> None:
        Run the import process.
    _import_mass_spectra(d_params: dict) -> None:
        Import the mass spectra from the file.
    get_lcms_obj() -> LCMSBase:
        Get the LCMSBase object.

    """

    def __init__(self, file_location: str, auto_process: bool = True):
        Thread.__init__(self)
        self.lcms: LCMSBase = LCMSBase(file_location)
        self.Bruker_Library = CreateObject("EDAL.MSAnalysis")
        self.res = self.Bruker_Library.Open(file_location)
        self.check_load_success()
        self._initial_scan_number: int = 1
        self._final_scan_number: int = self.get_scans_numbers()
        self.file_location: str = file_location
        self.auto_process: bool = auto_process

    @property
    def initial_scan_number(self) -> int:
        return self._initial_scan_number

    @initial_scan_number.setter
    def initial_scan_number(self, initial_scan_number: int) -> None:
        if self.check_scan(initial_scan_number):
            self._initial_scan_number = initial_scan_number
        else:
            raise Exception(
                f"startscan and finalscan should be less than {self.get_scans_numbers()}"
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
                f"startscan and finalscan should be less than {self.get_scans_numbers()}"
            )

    def check_scan(self, scan: int) -> bool:
        scan_numbers = self.get_scans_numbers()
        return scan <= scan_numbers

    def get_scans_numbers(self) -> int:
        scan_numbers = self.Bruker_Library.MSSpectrumCollection.Count
        return scan_numbers

    def get_polarity_mode(self, spectrum) -> int:
        polarity_symbol = spectrum.Polarity
        if polarity_symbol == 0:
            return 1
        elif polarity_symbol == 1:
            return -1
        else:
            raise IOError("Could not read mass spectrum polarity mode")

    def check_load_success(self) -> None:
        if self.res == 0:
            self.break_it = False
        else:
            raise ImportError(str(self.res))

    def get_bruker_tics(self) -> array:
        strAnalysisData = "SumIntensity"
        if self.Bruker_Library.HasAnalysisData(strAnalysisData):
            tics_array = self.Bruker_Library.GetAnalysisData(strAnalysisData)
            tics_array = array(tics_array)
        return tics_array

    def get_bruker_retention_time(self) -> array:
        strAnalysisData = "RetentionTime"
        if self.Bruker_Library.HasAnalysisData(strAnalysisData):
            tics_array = self.Bruker_Library.GetAnalysisData(strAnalysisData)
            tics_array = array(tics_array)
        else:
            tics_array = [0]
        return tics_array

    @staticmethod
    def get_data(spectra, scan: int) -> dict:
        spectrum = spectra[scan]
        is_profile = c_long(1)
        masslist = spectrum.GetMassIntensityValues(is_profile)
        data_dict = {
            Labels.mz: array(masslist[0]),
            Labels.abundance: array(masslist[1]),
            Labels.rp: None,
            Labels.s2n: None,
        }
        return data_dict

    def run(self) -> None:
        d_parameters = default_parameters(self.file_location)
        self._import_mass_spectra(d_parameters)

    def _import_mass_spectra(self, d_params: dict) -> None:
        spectra = self.Bruker_Library.MSSpectrumCollection
        list_rt = self.get_bruker_retention_time()
        list_Tics = self.get_bruker_tics()
        list_scans = list()
        for scan_number in range(self.initial_scan_number, self.final_scan_number + 1):
            if spectra[scan_number].MSMSStage == 1:
                d_params["label"] = Labels.bruker_profile
                d_params["polarity"] = self.get_polarity_mode(spectra[scan_number])
                d_params["rt"] = list_rt[scan_number - 1]
                d_params["scan_number"] = scan_number
                list_scans.append(scan_number)
                data_dict = self.get_data(spectra, scan_number)
                mass_spec = MassSpecProfile(
                    data_dict, d_params, auto_process=self.auto_process
                )
                mass_spec.process_mass_spec()
                self.lcms.add_mass_spectrum(mass_spec)
        self.lcms.retention_time = list_rt
        self.lcms.tic = list_Tics
        self.lcms.scans_number = list_scans

    def get_lcms_obj(self) -> LCMSBase:
        if self.lcms.get(self._initial_scan_number):
            return self.lcms
        else:
            raise Exception("Returning an empty LCMSBase object.")
