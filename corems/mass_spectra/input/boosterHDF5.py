__author__ = "Yuri E. Corilo"
__date__ = "Oct 29, 2019"

from threading import Thread
from pathlib import Path
from io import BytesIO

import h5py
from s3path import S3Path

from corems.encapsulation.constant import Labels
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecProfile
from corems.mass_spectra.factory.lc_class import LCMSBase
from corems.encapsulation.factory.parameters import default_parameters


class ReadHDF_BoosterMassSpectra(Thread):
    """
    Class for reading HDF5 files containing booster mass spectra.

    Parameters
    ----------
    file_location : Path or S3Path
        The full path to the HDF5 file.
    analyzer : str, optional
        The type of analyzer used for the mass spectra. Defaults to "ICR".
    instrument_label : str, optional
        The label of the instrument. Defaults to "21T".
    auto_process : bool, optional
        Whether to automatically process the mass spectra. Defaults to True.
    """

    def __init__(
        self,
        file_location: Path | S3Path,
        analyzer="ICR",
        instrument_label="21T",
        auto_process=True,
    ):
        """
        Initialize the ReadHDF_BoosterMassSpectra class.

        Parameters
        ----------
        file_location : Path or S3Path
            The full path to the HDF5 file.
        analyzer : str, optional
            The type of analyzer used for the mass spectra. Defaults to "ICR".
        instrument_label : str, optional
            The label of the instrument. Defaults to "21T".
        auto_process : bool, optional
            Whether to automatically process the mass spectra. Defaults to True.
        """
        Thread.__init__(self)

        if isinstance(file_location, str):
            # if obj is a string it defaults to create a Path obj, pass the S3Path if needed
            self.file_location = Path(file_location)

        self.lcms = LCMSBase(
            file_location, analyzer=analyzer, instrument_label=instrument_label
        )

        if isinstance(file_location, S3Path):
            data = BytesIO(file_location.open("rb").read())
        else:
            data = file_location

        self.hdf_obj = h5py.File(data, "r")

        self.list_scans = sorted([int(i) for i in list(self.hdf_obj.keys())])

        self.initial_scan_number = self.list_scans[0]

        self.final_scan_number = self.list_scans[-1]

        self.file_location = file_location

        self.auto_process = True

        self.analyzer = analyzer

        self.instrument_label = instrument_label

    def get_polarity(self, file_location: Path | S3Path, scan: int):
        """
        Get the polarity of a scan.

        Parameters
        ----------
        file_location : Path or S3Path
            The full path to the HDF5 file.
        scan : int
            The scan number.

        """
        if isinstance(file_location, S3Path):
            data = BytesIO(file_location.open("rb").read())
        else:
            data = file_location

        self.h5pydata = h5py.File(data, "r")

        self.scans = list(self.h5pydata.keys())

        polarity = self.get_attr_data(scan, "r_h_polarity")

        if polarity == "negative scan":
            return -1
        else:
            return +1

    def get_attr_data(self, scan, attr_srt):
        """
        Get the attribute data of a scan.

        Parameters
        ----------
        scan : int
            The scan number.
        attr_srt : str
            The attribute name.

        """
        return self.hdf_obj[str(scan)].attrs[attr_srt]

    def import_mass_spectra(self, d_params: dict):
        """
        Import the mass spectra from the HDF5 file.

        Parameters
        ----------
        d_params : dict
            The parameters for importing the mass spectra.
        """
        list_rt, list_tic = list(), list()

        for scan_number in self.list_scans:
            d_params["rt"] = list_rt.append(
                self.get_attr_data(scan_number, "r_h_start_time")
            )

            d_params["scan_number"] = scan_number

            d_params["label"] = Labels.booster_profile

            d_params["polarity"] = self.get_polarity(self.file_location, scan_number)

            d_params["Aterm"] = self.get_attr_data(scan_number, "r_cparams")[0]

            d_params["Bterm"] = self.get_attr_data(scan_number, "r_cparams")[1]

            d_params["analyzer"] = self.analyzer

            d_params["instrument_label"] = self.instrument_label

            list_rt.append(d_params["rt"])

            list_tic.append(self.get_attr_data(scan_number, "r_h_tic"))

            mass_spec = self.get_mass_spectrum(scan_number, d_params)

            self.lcms.add_mass_spectrum(mass_spec)

        self.lcms.retention_time = list_rt
        self.lcms.tic = list_tic
        self.lcms.scans_number = self.list_scans

    def get_mass_spectrum(self, scan: int, d_params: dict):
        """
        Get the mass spectrum for a scan.

        Parameters
        ----------
        scan : int
            The scan number.
        d_params : dict
            The parameters for creating the mass spectrum.

        """
        booster_data = self.hdf_obj[str(scan)]

        if booster_data.shape[0] != 2:
            raise NotImplementedError(
                "opening transient, needs read raw file here, get bandwidth, create transient class and then the mass spectrum"
            )
        else:
            data_dict = {
                Labels.mz: booster_data[0],
                Labels.abundance: booster_data[1],
                Labels.rp: None,
                Labels.s2n: None,
            }

            mass_spec = MassSpecProfile(
                data_dict, d_params, auto_process=self.auto_process
            )

        return mass_spec

    def run(self):
        """
        Run the thread to create the LCMS object.
        """
        d_parameters = default_parameters(self.file_location)
        self.import_mass_spectra(d_parameters)

    def get_lcms_obj(self):
        """
        Get the LCMS object.

        """
        if len(self.lcms) > 0:
            return self.lcms
        else:
            raise Exception("Returning an empty LCMS class")
