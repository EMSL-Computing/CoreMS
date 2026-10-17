__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"

import sys

from pathlib import Path
from threading import Thread
import warnings
from corems.mass_spectrum.input.massList import ReadCoremsMasslist
from corems.mass_spectra.factory.lc_class import LCMSBase


class ReadCoremsMassSpectraText(ReadCoremsMasslist, Thread):
    """
    Class for reading CoreMS mass spectra from a text file.

    Parameters
    ----------
    file_location : str, pathlib.Path, or s3path.S3Path
        Path object from pathlib containing the file location
    analyzer : str, optional
        Name of the analyzer, by default 'Unknown'
    instrument_label : str, optional
        Label of the instrument, by default 'Unknown'

    Attributes
    ----------
    lcms : LCMSBase
        LCMSBase object for storing the mass spectra data.

    Methods
    -------
    * get_scans_filename(). Get the filenames of all the scan files associated with the CoreMS file.
    * set_filepath_datatype_and_delimiter(file_path_obj). Set the file path, data type, and delimiter based on the file path object.
    * import_mass_spectra(). Import the mass spectra from the scan files and add them to the LCMSBase object.
    * run(). Run the import_mass_spectra method to create the LCMSBase object.
    * get_lcms_obj(). Get the LCMSBase object.
    """

    def __init__(self, file_location, analyzer="Unknown", instrument_label="Unknown"):
        if isinstance(file_location, str):
            # if obj is a string it defaults to create a Path obj, pass the S3Path if needed
            file_location = Path(file_location)

        if not file_location.exists():
            raise FileNotFoundError("%s not found" % file_location)

        if not file_location.suffix == ".corems":
            raise TypeError("%s is not a valid CoreMS file" % file_location)

        Thread.__init__(self)

        ReadCoremsMasslist.__init__(self, file_location)

        self.lcms = LCMSBase(
            self.file_location, analyzer=analyzer, instrument_label=instrument_label
        )

    def get_scans_filename(self) -> list:
        all_other = self.file_location.glob("*_scan*[!.json]")

        scans_filepath = [
            (file_path_obj.stem.split("scan")[1], file_path_obj)
            for file_path_obj in all_other
        ]

        scans_filepath.sort(key=lambda m: int(m[0]))

        return scans_filepath

    def set_filepath_datatype_and_delimiter(self, file_path_obj) -> None:
        self.file_location = file_path_obj

        if file_path_obj.suffix == ".pkl":
            self.data_type == "dataframe"

        else:
            if file_path_obj.suffix == ".csv":
                self.data_type == "txt"
                self.delimiter = ","

            elif file_path_obj.suffix == ".xlsx":
                self.data_type == "excel"
                self.delimiter = ","

            elif file_path_obj.suffix == ".txt":
                self.data_type == "txt"
                self.delimiter = "\t"
                warnings.warn("using tab as delimiter")
            else:
                raise NotImplementedError(
                    "%s data not yet supported " % file_path_obj.suffix
                )

    def import_mass_spectra(self) -> None:
        list_rt, list_tic, list_scan = list(), list(), list()

        for scan_number, file_path_obj in self.get_scans_filename():
            self.set_filepath_datatype_and_delimiter(file_path_obj)

            mass_spec = self.get_mass_spectrum(int(scan_number))

            list_scan.append(int(scan_number))

            list_rt.append(mass_spec.retention_time)

            list_tic.append(mass_spec.tic)

            self.lcms.add_mass_spectrum(mass_spec)

        self.lcms.retention_time = list_rt
        self.lcms.tic_list = list_tic  # TODO: check if this is correct
        self.lcms.scans_number = list_scan

    def run(self) -> None:
        """Creates the LCMS object and imports mass spectra."""

        self.import_mass_spectra()

    def get_lcms_obj(self) -> LCMSBase:
        """
        Returns the LCMSBase object associated with the massList.

        If the LCMSBase object is already initialized, it is returned.
        Otherwise, an exception is raised.

        Raises:
            Exception: If the LCMSBase object is not initialized.
        """
        if self.lcms:
            return self.lcms
        else:
            raise Exception("returning an empty lcms class")
