__author__ = "Yuri E. Corilo"
__date__ = "Nov 11, 2019"

from copy import deepcopy
from io import BytesIO
from pathlib import Path

import chardet
from bs4 import BeautifulSoup
from pandas import read_csv, read_excel, read_pickle
from pandas.core.frame import DataFrame
from s3path import S3Path

from corems.encapsulation.constant import Labels
from corems.encapsulation.factory.parameters import default_parameters
from corems.encapsulation.factory.processingSetting import DataInputSetting
from corems.encapsulation.input.parameter_from_json import (
    load_and_set_parameters_class,
    load_and_set_parameters_ms,
    load_and_set_toml_parameters_class,
)


class MassListBaseClass:
    """The MassListBaseClass object reads mass list data types and returns the mass spectrum obj

    Parameters
    ----------
    file_location : Path or S3Path
        Full data path.
    isCentroid : bool, optional
        Determines the mass spectrum data structure. If set to True, it assumes centroid mode. If set to False, it assumes profile mode and attempts to peak pick. Default is True.
    analyzer : str, optional
        The analyzer used for the mass spectrum. Default is 'Unknown'.
    instrument_label : str, optional
        The label of the instrument used for the mass spectrum. Default is 'Unknown'.
    sample_name : str, optional
        The name of the sample. Default is None.
    header_lines : int, optional
        The number of lines to skip in the file, including the column labels line. Default is 0.
    isThermoProfile : bool, optional
        Determines the number of expected columns in the file. If set to True, only m/z and intensity columns are expected. Signal-to-noise ratio (S/N) and resolving power (RP) will be calculated based on the data. Default is False.
    headerless : bool, optional
        If True, assumes that there are no headers present in the file (e.g., a .xy file from Bruker) and assumes two columns: m/z and intensity. Default is False.

    Attributes
    ----------
    parameters : DataInputSetting
        The data input settings for the mass spectrum.
    data_type : str
        The type of data in the file.
    delimiter : str
        The delimiter used to read text-based files.

    Methods
    -------
    * set_parameter_from_toml(parameters_path). Sets the data input settings from a TOML file.
    * set_parameter_from_json(parameters_path). Sets the data input settings from a JSON file.
    * get_dataframe(). Reads the file and returns the data as a pandas DataFrame.
    * load_settings(mass_spec_obj, output_parameters). Loads the settings for the mass spectrum.
    * get_output_parameters(polarity, scan_index=0). Returns the output parameters for the mass spectrum.
    * clean_data_frame(dataframe). Cleans the data frame by removing columns that are not in the expected columns set.

    """

    def __init__(
        self,
        file_location: Path | S3Path,
        isCentroid: bool = True,
        analyzer: str = "Unknown",
        instrument_label: str = "Unknown",
        sample_name: str = None,
        header_lines: int = 0,
        isThermoProfile: bool = False,
        headerless: bool = False,
    ):
        self.file_location = (
            Path(file_location) if isinstance(file_location, str) else file_location
        )

        if not self.file_location.exists():
            raise FileExistsError("File does not exist: %s" % file_location)

        # (newline="\n")

        self.header_lines = header_lines

        if isThermoProfile:
            self._expected_columns = {Labels.mz, Labels.abundance}

        else:
            self._expected_columns = {
                Labels.mz,
                Labels.abundance,
                Labels.s2n,
                Labels.rp,
            }

        self._delimiter = None

        self.isCentroid = isCentroid

        self.isThermoProfile = isThermoProfile

        self.headerless = headerless

        self._data_type = None

        self.analyzer = analyzer

        self.instrument_label = instrument_label

        self.sample_name = sample_name

        self._parameters = deepcopy(DataInputSetting())

    @property
    def parameters(self):
        return self._parameters

    @parameters.setter
    def parameters(self, instance_DataInputSetting):
        self._parameters = instance_DataInputSetting

    def set_parameter_from_toml(self, parameters_path):
        self._parameters = load_and_set_toml_parameters_class(
            "DataInput", self.parameters, parameters_path=parameters_path
        )

    def set_parameter_from_json(self, parameters_path):
        self._parameters = load_and_set_parameters_class(
            "DataInput", self.parameters, parameters_path=parameters_path
        )

    @property
    def data_type(self):
        return self._data_type

    @data_type.setter
    def data_type(self, data_type):
        self._data_type = data_type

    @property
    def delimiter(self):
        return self._delimiter

    @delimiter.setter
    def delimiter(self, delimiter):
        self._delimiter = delimiter

    def encoding_detector(self, file_location) -> str:
        """
        Detects the encoding of a file.

        Parameters
        --------
        file_location : str
            The location of the file to be analyzed.

        Returns
        --------
        str
            The detected encoding of the file.
        """

        with file_location.open("rb") as rawdata:
            result = chardet.detect(rawdata.read(10000))
        return result["encoding"]

    def set_data_type(self):
        """
        Set the data type and delimiter based on the file extension.

        Raises
        ------
        TypeError
            If the data type could not be automatically recognized.
        """
        if self.file_location.suffix == ".csv":
            self.data_type = "txt"
            self.delimiter = ","
        elif self.file_location.suffix == ".txt":
            self.data_type = "txt"
            self.delimiter = "\t"
        elif self.file_location.suffix == ".tsv":
            self.data_type = "txt"
            self.delimiter = "\t"
        elif self.file_location.suffix == ".xlsx":
            self.data_type = "excel"
        elif self.file_location.suffix == ".ascii":
            self.data_type = "txt"
            self.delimiter = "  "
        elif self.file_location.suffix == ".pkl":
            self.data_type = "dataframe"
        elif self.file_location.suffix == ".pks":
            self.data_type = "pks"
            self.delimiter = "          "
            self.header_lines = 9
        elif self.file_location.suffix == ".xml":
            self.data_type = "xml"
            # self.delimiter = None
            # self.header_lines = None
        elif self.file_location.suffix == ".xy":
            self.data_type = "txt"
            self.delimiter = " "
            self.header_lines = None
        else:
            raise TypeError(
                "Data type could not be automatically recognized for %s; please set data type and delimiter manually."
                % self.file_location.name
            )

    def get_dataframe(self) -> DataFrame:
        """
        Get the data as a pandas DataFrame.

        Returns
        -------
        pandas.DataFrame
            The data as a pandas DataFrame.

        Raises
        ------
        TypeError
            If the data type is not supported.
        """

        if not self.data_type or not self.delimiter:
            self.set_data_type()

        if isinstance(self.file_location, S3Path):
            data = BytesIO(self.file_location.open("rb").read())
        else:
            data = self.file_location

        if self.data_type == "txt":
            if self.headerless:
                dataframe = read_csv(
                    data,
                    skiprows=self.header_lines,
                    delimiter=self.delimiter,
                    header=None,
                    names=["m/z", "I"],
                    encoding=self.encoding_detector(self.file_location),
                    engine="python",
                )
            else:
                dataframe = read_csv(
                    data,
                    skiprows=self.header_lines,
                    delimiter=self.delimiter,
                    encoding=self.encoding_detector(self.file_location),
                    engine="python",
                )

        elif self.data_type == "pks":
            names = [
                "m/z",
                "I",
                "Scaled Peak Height",
                "Resolving Power",
                "Frequency",
                "S/N",
            ]
            clean_data = []
            with self.file_location.open() as maglabfile:
                for i in maglabfile.readlines()[8:-1]:
                    clean_data.append(i.split())
            dataframe = DataFrame(clean_data, columns=names)

        elif self.data_type == "dataframe":
            dataframe = read_pickle(data)

        elif self.data_type == "excel":
            dataframe = read_excel(data)

        elif self.data_type == "xml":
            dataframe = self.read_xml_peaks(data)

        else:
            raise TypeError("Data type %s is not supported" % self.data_type)

        return dataframe

    def load_settings(self, mass_spec_obj, output_parameters):
        """
        #TODO loading output parameters from json file is not functional
        Load settings from a JSON file and apply them to the given mass_spec_obj.

        Parameters
        ----------
        mass_spec_obj : MassSpec
            The mass spectrum object to apply the settings to.

        """
        import json
        import warnings

        settings_file_path = self.file_location.with_suffix(".json")

        if settings_file_path.exists():
            self._parameters = load_and_set_parameters_class(
                "DataInput", self._parameters, parameters_path=settings_file_path
            )

            load_and_set_parameters_ms(
                mass_spec_obj, parameters_path=settings_file_path
            )

        else:
            warnings.warn(
                "auto settings loading is enabled but could not locate the file:  %s. Please load the settings manually"
                % settings_file_path
            )

        # TODO this will load the setting from SettingCoreMS.json
        # coreMSHFD5 overrides this function to import the attrs stored in the h5 file
        # loaded_settings = {}
        # loaded_settings['MoleculaSearch'] = self.get_scan_group_attr_data(scan_index,  time_index, 'MoleculaSearchSetting')
        # loaded_settings['MassSpecPeak'] = self.get_scan_group_attr_data(scan_index,  time_index, 'MassSpecPeakSetting')

        # loaded_settings['MassSpectrum'] = self.get_scan_group_attr_data(scan_index, time_index, 'MassSpectrumSetting')
        # loaded_settings['Transient'] = self.get_scan_group_attr_data(scan_index, time_index, 'TransientSetting')

    def get_output_parameters(self, polarity: int, scan_index: int = 0) -> dict:
        """
        Get the output parameters for the mass spectrum.

        Parameters
        ----------
        polarity : int
            The polarity of the mass spectrum +1 or -1.
        scan_index : int, optional
            The index of the scan. Default is 0.

        Returns
        -------
        dict
            A dictionary containing the output parameters.

        """
        from copy import deepcopy

        output_parameters = default_parameters(self.file_location)

        if self.isCentroid:
            output_parameters["label"] = Labels.corems_centroid
        else:
            output_parameters["label"] = Labels.bruker_profile

        output_parameters["analyzer"] = self.analyzer

        output_parameters["instrument_label"] = self.instrument_label

        output_parameters["sample_name"] = self.sample_name

        output_parameters["Aterm"] = None

        output_parameters["Bterm"] = None

        output_parameters["Cterm"] = None

        output_parameters["polarity"] = polarity

        # scan_number and rt will be need to lc ms====

        output_parameters["mobility_scan"] = 0

        output_parameters["mobility_rt"] = 0

        output_parameters["scan_number"] = scan_index

        output_parameters["rt"] = 0

        return output_parameters

    def clean_data_frame(self, dataframe):
        """
        Clean the input dataframe by removing columns that are not expected.

        Parameters
        ----------
        pandas.DataFrame
            The input dataframe to be cleaned.

        """

        for column_name in dataframe.columns:
            expected_column_name = self.parameters.header_translate.get(column_name)
            if expected_column_name not in self._expected_columns:
                del dataframe[column_name]

    def check_columns(self, header_labels: list[str]):
        """
        Check if the given header labels match the expected columns.

        Parameters
        ----------
        header_labels : list
            The header labels to be checked.

        Raises
        ------
        Exception
            If any expected column is not found in the header labels.
        """
        found_label = set()

        for label in header_labels:
            if not label in self._expected_columns:
                user_column_name = self.parameters.header_translate.get(label)
                if user_column_name in self._expected_columns:
                    found_label.add(user_column_name)
            else:
                found_label.add(label)

        not_found = self._expected_columns - found_label

        if len(not_found) > 0:
            raise Exception(
                "Please make sure to include the columns %s" % ", ".join(not_found)
            )

    def read_xml_peaks(self, data: str) -> DataFrame:
        """
        Read peaks from a Bruker .xml file and return a pandas DataFrame.

        Parameters
        ----------
        data : str
            The path to the .xml file.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the peak data with columns: 'm/z', 'I', 'Resolving Power', 'Area', 'S/N', 'fwhm'.
        """
        from numpy import nan

        with open(data, "r") as file:
            content = file.readlines()
            content = "".join(content)
            bs_content = BeautifulSoup(content, features="xml")
        peaks_xml = bs_content.find_all("pk")

        # initialise lists of the peak variables
        areas = []
        fwhms = []
        intensities = []
        mzs = []
        res = []
        sn = []
        # iterate through the peaks appending to each list
        for peak in peaks_xml:
            areas.append(
                float(peak.get("a", nan))
            )  # Use a default value if key 'a' is missing
            fwhms.append(
                float(peak.get("fwhm", nan))
            )  # Use a default value if key 'fwhm' is missing
            intensities.append(
                float(peak.get("i", nan))
            )  # Use a default value if key 'i' is missing
            mzs.append(
                float(peak.get("mz", nan))
            )  # Use a default value if key 'mz' is missing
            res.append(
                float(peak.get("res", nan))
            )  # Use a default value if key 'res' is missing
            sn.append(
                float(peak.get("sn", nan))
            )  # Use a default value if key 'sn' is missing

        # Compile pandas dataframe of these values
        names = ["m/z", "I", "Resolving Power", "Area", "S/N", "fwhm"]
        df = DataFrame(columns=names, dtype=float)
        df["m/z"] = mzs
        df["I"] = intensities
        df["Resolving Power"] = res
        df["Area"] = areas
        df["S/N"] = sn
        df["fwhm"] = fwhms
        return df

    def get_xml_polarity(self):
        """
        Get the polarity from an XML peaklist.

        Returns
        -------
        int
            The polarity of the XML peaklist. Returns -1 for negative polarity, +1 for positive polarity.

        Raises
        ------
        Exception
            If the data type is not XML peaklist in Bruker format or if the polarity is unhandled.
        """

        # Check its an actual xml
        if not self.data_type or not self.delimiter:
            self.set_data_type()

        if isinstance(self.file_location, S3Path):
            # data = self.file_location.open('rb').read()
            data = BytesIO(self.file_location.open("rb").read())

        else:
            data = self.file_location

        if self.data_type != "xml":
            raise Exception("This function is only for XML peaklists (Bruker format)")

        with open(data, "r") as file:
            content = file.readlines()
            content = "".join(content)
            bs_content = BeautifulSoup(content, features="xml")
        polarity = bs_content.find_all("ms_spectrum")[0]["polarity"]
        if polarity == "-":
            return -1
        elif polarity == "+":
            return +1
        else:
            raise Exception("Polarity %s unhandled" % polarity)
