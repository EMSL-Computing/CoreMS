__author__ = "Katherine Heal"
__date__ = "November 8, 2023"

from abc import ABC, abstractmethod


class SpectraParserInterface(ABC):
    """
    Interface for parsing mass spectra data into MassSpectraBase objects. This is an abstract class and should not be instantiated directly.

    Methods
    -------
    * load().
        Load mass spectra data.
    * run().
        Parse mass spectra data.
    * get_mass_spectra_obj().
        Return MassSpectraBase object with several attributes populated
    * get_mass_spectrum_from_scan(scan_number).
        Return MassSpecBase data object from scan number.
    """

    @abstractmethod
    def load(self):
        """
        Load mass spectra data.
        """
        pass

    @abstractmethod
    def run(self):
        """
        Parse mass spectra data.
        """
        pass

    @abstractmethod
    def get_mass_spectra_obj(self):
        """
        Return mass spectra data object.
        """
        pass

    @abstractmethod
    def get_mass_spectrum_from_scan(
        self, scan_number, spectrum_mode, auto_process=True
    ):
        """
        Return mass spectrum data object from scan number.
        """
        pass
