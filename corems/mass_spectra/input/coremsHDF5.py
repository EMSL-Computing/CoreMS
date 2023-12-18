__author__ = "Yuri E. Corilo"
__date__ = "Oct 29, 2019"


from threading import Thread

import h5py

from corems.encapsulation.constant import Labels
from corems.mass_spectrum.input.coremsHDF5 import ReadCoreMSHDF_MassSpectrum
from corems.mass_spectra.factory.LC_Class import LCMSBase


class ReadCoreMSHDF_MassSpectra(ReadCoreMSHDF_MassSpectrum, Thread):
    """
    Class for reading CoreMS HDF5 mass spectra.

    Parameters
    ----------
    file_location : str
        The file location of the CoreMS HDF5 file.

    Attributes
    ----------
    lcms : LCMSBase
        The LCMSBase object for storing mass spectra.
    list_scans : list
        The list of scan numbers in the HDF5 file.

    Methods
    -------
    * import_mass_spectra()
        Imports the mass spectra from the HDF5 file.
    * run()
        Creates the LCMS object by importing the mass spectra.
    * get_lcms_obj()
        Returns the LCMS object.

    Raises
    ------
    Exception
        If the LCMS object is empty.

    """

    def __init__(self, file_location: str):
        Thread.__init__(self)
        ReadCoreMSHDF_MassSpectrum.__init__(self, file_location)
        self.lcms = LCMSBase(self.file_location)
        self.list_scans = sorted([int(i) for i in list(self.h5pydata.keys())])

    def import_mass_spectra(self) -> None:
        """
        Imports the mass spectra from the HDF5 file.
        """
        list_rt, list_tic = list(), list()
        for scan_number in self.list_scans:
            mass_spec = self.get_mass_spectrum(scan_number)
            list_rt.append(mass_spec.retention_time)
            list_tic.append(mass_spec.tic)
            self.lcms.add_mass_spectrum(mass_spec)
        self.lcms.retention_time = list_rt
        self.lcms.tic = list_tic
        self.lcms.scans_number = self.list_scans

    def run(self) -> None:
        """
        Creates the LCMS object by importing the mass spectra.
        """
        self.import_mass_spectra()

    def get_lcms_obj(self) -> LCMSBase:
        """
        Returns the LCMS object.

        Returns
        -------
        LCMSBase
            The LCMS object.

        Raises
        ------
        Exception
            If the LCMS object is empty.
        """
        if self.lcms:
            return self.lcms
        else:
            raise Exception("Returning an empty LCMS class")
