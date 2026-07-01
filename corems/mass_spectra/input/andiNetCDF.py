__author__ = "Yuri E. Corilo"
__date__ = "Feb 12, 2020"

from pathlib import Path
from threading import Thread
# from io import BytesIO

from netCDF4 import Dataset
from s3path import S3Path

from corems.encapsulation.constant import Labels
from corems.encapsulation.factory.parameters import default_parameters
from corems.mass_spectra.factory.GC_Class import GCMSBase
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecCentroidLowRes


class ReadAndiNetCDF(Thread):
    """
    A class for reading AndiNetCDF files and extracting mass spectra data.

    Parameters
    -----------
    file_location : str or Path
            The location of the AndiNetCDF file.
    analyzer : str, optional
            The type of analyzer used (default is 'Quadruple').
    instrument_label : str, optional
            The label of the instrument (default is 'GCMS-Agilent').
    auto_process : bool, optional
            Whether to automatically process the data (default is True).

    Attributes
    -----------
    file_location : Path
            The path to the AndiNetCDF file.
    net_cdf_obj : Dataset
            The NetCDF dataset object.
    ionization_type : str
            The ionization type used in the experiment.
    experiment_type : str
            The type of experiment.
    list_scans : range
            The range of scan numbers in the dataset.
    initial_scan_number : int
            The number of the initial scan.
    final_scan_number : int
            The number of the final scan.
    analyzer : str
            The type of analyzer used.
    instrument_label : str
            The label of the instrument.
    gcms : GCMSBase
            The GCMSBase object for storing mass spectra data.

    Methods
    --------
    * polarity().
            Get the polarity of the ionization.
    * get_mass_spectrum(mz, abun, rp, d_params).
            Add a mass spectrum to the GCMSBase object.
    * run().
            Populate the GCMSBase object with mass spectra data.
    * import_mass_spectra(d_params).
            Import mass spectra data from the AndiNetCDF file.
    * get_gcms_obj().
            Get the GCMSBase object.

    """

    def __init__(
        self,
        file_location: str | Path,
        analyzer="Quadruple",
        instrument_label="GCMS-Agilent",
        auto_process=True,
    ):
        Thread.__init__(self)

        if isinstance(file_location, str):
            self.file_location = Path(file_location)
        else:
            self.file_location = file_location

        if not self.file_location.exists():
            raise FileNotFoundError("File does not exist at %s", file_location)

        if isinstance(file_location, S3Path):
            bytes_io = self.file_location.open("rb").read()
            self.net_cdf_obj = Dataset(
                self.file_location.name,
                "r",
                diskless=True,
                memory=bytes_io,
                format="NETCDF3_CLASSIC",
            )
        else:
            self.net_cdf_obj = Dataset(
                self.file_location, "r", format="NETCDF3_CLASSIC"
            )

        self.ionization_type = self.net_cdf_obj.test_ionization_mode
        self.experiment_type = self.net_cdf_obj.experiment_type
        self.list_scans = range(
            len(self.net_cdf_obj.variables.get("actual_scan_number")[:])
        )
        self.initial_scan_number = self.list_scans[0]
        self.final_scan_number = self.list_scans[-1]
        self.analyzer = analyzer
        self.instrument_label = instrument_label
        self.gcms = GCMSBase(self.file_location, analyzer, instrument_label)

    @property
    def polarity(self):
        """
        Get the polarity of the ionization.

        """
        polarity = str(self.net_cdf_obj.test_ionization_polarity)
        if polarity == "Positive Polarity":
            return +1
        else:
            return -1

    def get_mass_spectrum(self, mz, abun, rp, d_params):
        """
        Add a mass spectrum to the GCMSBase object.

        Parameters
        -----------
        mz : array-like
                The m/z values of the mass spectrum.
        abun : array-like
                The abundance values of the mass spectrum.
        rp : array-like
                The resolution values of the mass spectrum.
        d_params : dict
                Additional parameters for the mass spectrum.

        """
        data_dict = {
            Labels.mz: mz,
            Labels.abundance: abun,
            Labels.rp: rp,
            Labels.s2n: None,
        }
        mass_spec = MassSpecCentroidLowRes(data_dict, d_params)
        self.gcms.add_mass_spectrum(mass_spec)

    def run(self):
        """
        Populate the GCMSBase object with mass spectra data.
        """
        d_parameters = default_parameters(self.file_location)
        self.import_mass_spectra(d_parameters)

    def import_mass_spectra(self, d_params):
        """
        Import mass spectra data from the AndiNetCDF file.

        Parameters
        -----------
        d_params : dict
                Additional parameters for the mass spectra.

        """
        ms_datapoints_per_scans = self.net_cdf_obj.variables.get("point_count")[:]
        list_tic = self.net_cdf_obj.variables.get("total_intensity")[:]
        list_rt = self.net_cdf_obj.variables.get("scan_acquisition_time")[:] / 60
        mass_values = self.net_cdf_obj.variables.get("mass_values")[:]
        intensity_values = self.net_cdf_obj.variables.get("intensity_values")[:]
        resolution = self.net_cdf_obj.variables.get("resolution")[:]
        individual_rp = len(mass_values) == len(resolution)
        finish_location = -1
        for scan_index in self.list_scans:
            datapoints = ms_datapoints_per_scans[scan_index]
            finish_location += datapoints
            start_location = finish_location - datapoints + 1
            d_params["rt"] = list_rt[scan_index]
            d_params["scan_number"] = scan_index
            d_params["label"] = Labels.gcms_centroid
            d_params["polarity"] = self.polarity
            d_params["analyzer"] = self.analyzer
            d_params["instrument_label"] = self.instrument_label
            mz = mass_values[start_location:finish_location]
            abun = intensity_values[start_location:finish_location]
            if individual_rp:
                rp = resolution[start_location:finish_location]
            else:
                rp = [resolution[scan_index]] * datapoints
            self.get_mass_spectrum(mz, abun, rp, d_params)
        self.gcms.retention_time = list_rt
        self.gcms.tic = list_tic
        self.gcms.scans_number = self.list_scans

    def get_gcms_obj(self):
        """
        Get the GCMSBase object.

        """
        return self.gcms
