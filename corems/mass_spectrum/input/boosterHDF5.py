from io import BytesIO

import h5py
from s3path import S3Path

from corems.encapsulation.constant import Labels
from corems.encapsulation.factory.parameters import default_parameters
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecProfile
from corems.mass_spectrum.input.baseClass import MassListBaseClass


class ReadHDF_BoosterMassSpectrum(MassListBaseClass):
    """The ReadHDF_BoosterMassSpectrum class parses the mass spectrum data from an HDF file and generate a mass spectrum object.

    Parameters
    ----------
    file_location : str
        The path to the HDF file.
    isCentroid : bool, optional
        Specifies whether the mass spectrum is centroided or not. Default is False.

    Attributes
    ----------
    polarity : int
        The polarity of the mass spectrum.
    h5pydata : h5py.File
        The HDF file object.
    scans : list
        The list of scan names in the HDF file.

    Methods
    -------
    * get_data_profile(mz, abundance, auto_process). Returns a MassSpecProfile object from the given m/z and abundance arrays.
    * get_attr_data(scan, attr_srt). Returns the attribute value for the given scan and attribute name.
    * get_polarity(file_location). Returns the polarity of the mass spectrum.
    * get_mass_spectrum(auto_process). Returns the mass spectrum as a MassSpecProfile object.
    * get_output_parameters(). Returns the default output parameters for the mass spectrum.
    """

    def __init__(self, file_location, isCentroid=False):
        self.polarity = self.get_polarity(file_location)
        super().__init__(file_location, isCentroid=False)

    def get_data_profile(self, mz, abundance, auto_process) -> MassSpecProfile:
        """
        Returns a MassSpecProfile object from the given m/z and abundance arrays.

        Parameters
        ----------
        mz : array_like
            The m/z values.
        abundance : array_like
            The abundance values.
        auto_process : bool
            Specifies whether to automatically process the mass spectrum.

        Returns
        -------
        MassSpecProfile
            The MassSpecProfile object.

        """
        data_dict = {Labels.mz: mz, Labels.abundance: abundance}
        output_parameters = self.get_output_parameters()
        return MassSpecProfile(data_dict, output_parameters, auto_process=auto_process)

    def get_attr_data(self, scan, attr_srt):
        """
        Returns the attribute value for the given scan and attribute name.

        Parameters
        ----------
        scan : int
            The scan index.
        attr_srt : str
            The attribute name.

        Returns
        -------
        object
            The attribute value.

        """
        return self.h5pydata[self.scans[scan]].attrs[attr_srt]

    def get_polarity(self, file_location: str | S3Path) -> int:
        """
        Returns the polarity of the mass spectrum.

        Parameters
        ----------
        file_location : str
            The path to the HDF file.

        Returns
        -------
        int
            The polarity of the mass spectrum.

        """
        if isinstance(file_location, S3Path):
            data = BytesIO(file_location.open("rb").read())
        else:
            data = file_location

        self.h5pydata = h5py.File(data, "r")
        self.scans = list(self.h5pydata.keys())

        polarity = self.get_attr_data(0, "r_h_polarity")

        if polarity == "negative scan":
            return -1
        else:
            return +1

    def get_mass_spectrum(self, auto_process=True) -> MassSpecProfile:
        """
        Returns the mass spectrum as a MassSpecProfile object.

        Parameters
        ----------
        auto_process : bool, optional
            Specifies whether to automatically process the mass spectrum. Default is True.

        Returns
        -------
        MassSpecProfile
            The MassSpecProfile object.

        """
        if len(self.scans) == 1:
            booster_data = self.h5pydata[self.scans[0]]

            if self.isCentroid:
                raise NotImplementedError
            else:
                mz = booster_data[0]
                abun = booster_data[1]
                return self.get_data_profile(mz, abun, auto_process)

    def get_output_parameters(self) -> dict:
        """
        Returns the default output parameters for the mass spectrum.

        Returns
        -------
        dict
            The default output parameters.

        """
        d_params = default_parameters(self.file_location)
        d_params["polarity"] = self.polarity
        d_params["filename_path"] = self.file_location
        d_params["mobility_scan"] = 0
        d_params["mobility_rt"] = 0
        d_params["scan_number"] = 0
        d_params["rt"] = self.get_attr_data(0, "r_h_start_time")
        d_params["label"] = Labels.booster_profile
        d_params["Aterm"] = self.get_attr_data(0, "r_cparams")[0]
        d_params["Bterm"] = self.get_attr_data(0, "r_cparams")[1]
        return d_params
