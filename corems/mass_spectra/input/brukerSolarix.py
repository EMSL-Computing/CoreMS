__author__ = "Yuri E. Corilo"
__date__ = "Oct 29, 2019"

from threading import Thread
from pathlib import Path
from s3path import S3Path

# import h5py

# from corems.encapsulation.constant import Labels
from corems.mass_spectra.factory.LC_Class import LCMSBase

# from corems.encapsulation.factory.parameters import default_parameters
from corems.transient.input.brukerSolarix import ReadBrukerSolarix


class ReadBruker_SolarixTransientMassSpectra(Thread):
    """
    Class for reading Bruker Solarix Transient Mass Spectra.

    Parameters
    ----------
    d_directory_location : str, pathlib.Path, or s3path.S3Path
        Path object from pathlib containing the file location.
    analyzer : str, optional
        Type of analyzer used in the mass spectrometer. Defaults to "ICR".
    instrument_label : str, optional
        Label for the instrument. Defaults to "15T".
    auto_process : bool, optional
        Flag indicating whether to automatically process the mass spectra. Defaults to True.
    keep_profile : bool, optional
        Flag indicating whether to keep the profile data in the mass spectra. Defaults to False.
    """

    def __init__(
        self,
        d_directory_location: str | Path | S3Path,
        analyzer="ICR",
        instrument_label="15T",
        auto_process=True,
        keep_profile=False,
    ):
        
        Thread.__init__(self)

        if isinstance(d_directory_location, str):
            # if obj is a string it defaults to create a Path obj, pass the S3Path if needed
            d_directory_location = Path(d_directory_location)

        if not d_directory_location.exists():
            raise FileNotFoundError("File does not exist: " + str(d_directory_location))

        self.scan_attr = d_directory_location / "scan.xml"

        if not self.scan_attr.exists():
            raise FileExistsError(
                "%s does not seem to be a valid Solarix Mass Spectra Experiment,\
                                maybe an Imaging experiment?\
                                please ReadBruker_SolarixTransientImage class for Imaging dataset "
                % d_directory_location
            )

        self.lcms = LCMSBase(d_directory_location, analyzer, instrument_label)

        self.auto_process = auto_process
        self.keep_profile = keep_profile

    def get_scan_attr(self) -> dict:
        """
        Get the scan attributes from the scan.xml file.

        Returns
        -------
        dict
            Dictionary containing the scan number as key and a tuple of retention time and TIC as value.
        """
        from bs4 import BeautifulSoup

        soup = BeautifulSoup(self.scan_attr.open(), "xml")

        list_rt = [float(rt.text) for rt in soup.find_all("minutes")]
        list_tic = [float(tic.text) for tic in soup.find_all("tic")]
        list_scan = [int(scan.text) for scan in soup.find_all("count")]

        dict_scan_rt_tic = dict(zip(list_scan, zip(list_rt, list_tic)))

        return dict_scan_rt_tic

    def import_mass_spectra(self) -> None:
        """
        Import the mass spectra from the scan.xml file.
        """
        dict_scan_rt_tic = self.get_scan_attr()

        list_rt, list_tic = (
            list(),
            list(),
        )

        list_scans = sorted(list(dict_scan_rt_tic.keys()))

        for scan_number in list_scans:
            mass_spec = self.get_mass_spectrum(scan_number)

            self.lcms.add_mass_spectrum(mass_spec)

            list_rt.append(dict_scan_rt_tic.get(scan_number)[0])

            list_tic.append(dict_scan_rt_tic.get(scan_number)[1])

        self.lcms.retention_time = list_rt
        self.lcms.tic = list_tic
        self.lcms.scans_number = list_scans

    def get_mass_spectrum(self, scan_number: int):
        """
        Get the mass spectrum for a given scan number.

        Parameters
        ----------
        scan_number : int
            Scan number.

        """
        bruker_reader = ReadBrukerSolarix(self.lcms.file_location)

        bruker_transient = bruker_reader.get_transient(scan_number)

        mass_spec = bruker_transient.get_mass_spectrum(
            plot_result=False,
            auto_process=self.auto_process,
            keep_profile=self.keep_profile,
        )

        return mass_spec

    def run(self):
        """
        Run the import_mass_spectra method.
        """
        self.import_mass_spectra()

    def get_lcms_obj(self):
        """
        Get the LCMSBase object.

        Raises
        ------
        Exception
            If the LCMSBase object is empty.
        """
        if self.lcms:
            return self.lcms
        else:
            raise Exception("Returning an empty LCMSBase class.")
