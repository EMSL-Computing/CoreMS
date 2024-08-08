import json

from pandas import DataFrame
import h5py
from io import BytesIO
from s3path import S3Path

from corems.encapsulation.input.parameter_from_json import _set_dict_data_ms
from corems.mass_spectrum.input.massList import ReadCoremsMasslist
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecCentroid
from corems.encapsulation.constant import Labels
from corems.encapsulation.factory.parameters import default_parameters


class ReadCoreMSHDF_MassSpectrum(ReadCoremsMasslist):
    """Class for reading mass spectrum data from a CoreMS HDF5 file.

    Attributes
    ----------
    h5pydata : h5py.File
        The HDF5 file object.
    scans : list
        List of scan labels in the HDF5 file.

    Parameters
    ----------
    file_location : str or S3Path
        The path to the CoreMS HDF5 file.

    Methods
    -------
    * load_raw_data(mass_spectrum, scan_index=0) Load raw data into the mass spectrum object.
    * get_mass_spectrum(scan_number=0, time_index=-1, auto_process=True, load_settings=True, load_raw=True).Get a mass spectrum object.
    * load_settings(mass_spectrum, scan_index=0, time_index=-1). Load settings into the mass spectrum object.
    * get_dataframe(scan_index=0, time_index=-1). Get a pandas DataFrame representing the mass spectrum.
    * get_time_index_to_pull(scan_label, time_index). Get the time index to pull from the HDF5 file.
    * get_high_level_attr_data(attr_str). Get high-level attribute data from the HDF5 file.
    * get_scan_group_attr_data(scan_index, time_index, attr_group, attr_srt=None). Get scan group attribute data from the HDF5 file.
    * get_raw_data_attr_data(scan_index, attr_group, attr_str). Get raw data attribute data from the HDF5 file.
    * get_output_parameters(polarity, scan_index=0). Get the output parameters for the mass spectrum.
    """

    def __init__(self, file_location):
        super().__init__(file_location)

        if isinstance(self.file_location, S3Path):
            data = BytesIO(self.file_location.open("rb").read())
        else:
            data = self.file_location

        self.h5pydata = h5py.File(data, "r")

        self.scans = list(self.h5pydata.keys())

        print(self.scans)

    def load_raw_data(self, mass_spectrum, scan_index=0):
        """
        Load raw data into the mass spectrum object.

        Parameters
        ----------
        mass_spectrum : MassSpecCentroid
            The mass spectrum object to load the raw data into.
        scan_index : int, optional
            The index of the scan to load the raw data from. Default is 0.
        """

        scan_label = self.scans[scan_index]

        mz_profile = self.h5pydata[scan_label]["raw_ms"][0]

        abundance_profile = self.h5pydata[scan_label]["raw_ms"][1]

        mass_spectrum.mz_exp_profile = mz_profile

        mass_spectrum.abundance_profile = abundance_profile

    def get_mass_spectrum(
        self,
        scan_number=0,
        time_index=-1,
        auto_process=True,
        load_settings=True,
        load_raw=True,
    ):
        """
        Get a mass spectrum object.

        Parameters
        ----------
        scan_number : int, optional
            The index of the scan to retrieve the mass spectrum from. Default is 0.
        time_index : int, optional
            The index of the time point to retrieve the mass spectrum from. Default is -1.
        auto_process : bool, optional
            Whether to automatically process the mass spectrum. Default is True.
        load_settings : bool, optional
            Whether to load the settings into the mass spectrum object. Default is True.
        load_raw : bool, optional
            Whether to load the raw data into the mass spectrum object. Default is True.

        Returns
        -------
        MassSpecCentroid
            The mass spectrum object.
        """

        dataframe = self.get_dataframe(scan_number, time_index=time_index)

        if not set(
            ["H/C", "O/C", "Heteroatom Class", "Ion Type", "Is Isotopologue"]
        ).issubset(dataframe.columns):
            raise ValueError(
                "%s it is not a valid CoreMS file" % str(self.file_location)
            )

        dataframe.rename(columns=self.parameters.header_translate, inplace=True)

        polarity = dataframe["Ion Charge"].values[0]

        output_parameters = self.get_output_parameters(polarity, scan_index=scan_number)

        mass_spec_obj = MassSpecCentroid(
            dataframe.to_dict(orient="list"), output_parameters
        )

        if load_settings:
            self.load_settings(
                mass_spec_obj, scan_index=scan_number, time_index=time_index
            )

        if load_raw:
            self.load_raw_data(mass_spec_obj, scan_index=scan_number)

        self.add_molecular_formula(mass_spec_obj, dataframe)

        return mass_spec_obj

    def load_settings(self, mass_spectrum, scan_index=0, time_index=-1):
        """
        Load settings into the mass spectrum object.

        Parameters
        ----------
        mass_spectrum : MassSpecCentroid
            The mass spectrum object to load the settings into.
        scan_index : int, optional
            The index of the scan to load the settings from. Default is 0.
        time_index : int, optional
            The index of the time point to load the settings from. Default is -1.
        """

        loaded_settings = {}
        loaded_settings["MoleculaSearch"] = self.get_scan_group_attr_data(
            scan_index, time_index, "MoleculaSearchSetting"
        )
        loaded_settings["MassSpecPeak"] = self.get_scan_group_attr_data(
            scan_index, time_index, "MassSpecPeakSetting"
        )
        loaded_settings["MassSpectrum"] = self.get_scan_group_attr_data(
            scan_index, time_index, "MassSpectrumSetting"
        )
        loaded_settings["Transient"] = self.get_scan_group_attr_data(
            scan_index, time_index, "TransientSetting"
        )

        _set_dict_data_ms(loaded_settings, mass_spectrum)

    def get_dataframe(self, scan_index=0, time_index=-1):
        """
        Get a pandas DataFrame representing the mass spectrum.

        Parameters
        ----------
        scan_index : int, optional
            The index of the scan to retrieve the DataFrame from. Default is 0.
        time_index : int, optional
            The index of the time point to retrieve the DataFrame from. Default is -1.

        Returns
        -------
        DataFrame
            The pandas DataFrame representing the mass spectrum.
        """

        columnsLabels = self.get_scan_group_attr_data(
            scan_index, time_index, "ColumnsLabels"
        )

        scan_label = self.scans[scan_index]

        index_to_pull = self.get_time_index_to_pull(scan_label, time_index)

        corems_table_data = self.h5pydata[scan_label][index_to_pull]

        list_dict = []
        for row in corems_table_data:
            data_dict = {}
            for data_index, data in enumerate(row):
                label = columnsLabels[data_index]
                # if data starts with a b' it is a byte string, so decode it
                if isinstance(data, bytes):
                    data = data.decode("utf-8")
                if data == "nan":
                    data = None
                data_dict[label] = data

            list_dict.append(data_dict)

        return DataFrame(list_dict)

    def get_time_index_to_pull(self, scan_label, time_index):
        """
        Get the time index to pull from the HDF5 file.

        Parameters
        ----------
        scan_label : str
            The label of the scan.
        time_index : int
            The index of the time point.

        Returns
        -------
        str
            The time index to pull.
        """

        time_data = sorted(
            [(i, int(i)) for i in self.h5pydata[scan_label].keys() if i != "raw_ms"],
            key=lambda m: m[1],
        )

        index_to_pull = time_data[time_index][0]

        return index_to_pull

    def get_high_level_attr_data(self, attr_str):
        """
        Get high-level attribute data from the HDF5 file.

        Parameters
        ----------
        attr_str : str
            The attribute string.

        Returns
        -------
        dict
            The attribute data.

        Raises
        ------
        KeyError
            If the attribute string is not found in the HDF5 file.
        """

        return self.h5pydata.attrs[attr_str]

    def get_scan_group_attr_data(
        self, scan_index, time_index, attr_group, attr_srt=None
    ):
        """
        Get scan group attribute data from the HDF5 file.

        Parameters
        ----------
        scan_index : int
            The index of the scan.
        time_index : int
            The index of the time point.
        attr_group : str
            The attribute group.
        attr_srt : str, optional
            The attribute string. Default is None.

        Returns
        -------
        dict
            The attribute data.

        Notes
        -----
        This method retrieves attribute data from the HDF5 file for a specific scan and time point.
        The attribute data is stored in the specified attribute group.
        If an attribute string is provided, only the corresponding attribute value is returned.
        If no attribute string is provided, all attribute data in the group is returned as a dictionary.
        """

        scan_label = self.scans[scan_index]

        index_to_pull = self.get_time_index_to_pull(scan_label, time_index)

        if attr_srt:
            return json.loads(
                self.h5pydata[scan_label][index_to_pull].attrs[attr_group]
            )[attr_srt]

        else:
            data = self.h5pydata[scan_label][index_to_pull].attrs.get(attr_group)
            if data:
                return json.loads(data)
            else:
                return {}

    def get_raw_data_attr_data(self, scan_index, attr_group, attr_str):
        """
        Get raw data attribute data from the HDF5 file.

        Parameters
        ----------
        scan_index : int
            The index of the scan.
        attr_group : str
            The attribute group.
        attr_str : str
            The attribute string.

        Returns
        -------
        dict
            The attribute data.

        Raises
        ------
        KeyError
            If the attribute string is not found in the attribute group.

        Notes
        -----
        This method retrieves the attribute data associated with a specific scan, attribute group, and attribute string
        from the HDF5 file. It returns the attribute data as a dictionary.

        Example usage:
        >>> data = get_raw_data_attr_data(0, "group1", "attribute1")
        >>> print(data)
        {'key1': 'value1', 'key2': 'value2'}
        """
        scan_label = self.scans[scan_index]
        try:
            json.loads(self.h5pydata[scan_label]["raw_ms"].attrs[attr_group])[attr_str]
        except KeyError:
            attr_str = attr_str.replace("baseline", "baselise")
        return json.loads(self.h5pydata[scan_label]["raw_ms"].attrs[attr_group])[attr_str]

    def get_output_parameters(self, polarity, scan_index=0):
        """
        Get the output parameters for the mass spectrum.

        Parameters
        ----------
        polarity : str
            The polarity of the mass spectrum.
        scan_index : int, optional
            The index of the scan. Default is 0.

        Returns
        -------
        dict
            The output parameters.
        """

        d_params = default_parameters(self.file_location)
        d_params["filename_path"] = self.file_location
        d_params["scan_number"] = int(self.scans[scan_index])
        d_params['polarity'] = self.get_raw_data_attr_data( scan_index, 'MassSpecAttrs', 'polarity')
        d_params['rt'] =     self.get_raw_data_attr_data( scan_index, 'MassSpecAttrs', 'rt')
        
        d_params['tic'] =  self.get_raw_data_attr_data( scan_index, 'MassSpecAttrs', 'tic')
        
        d_params['mobility_scan'] =    self.get_raw_data_attr_data( scan_index, 'MassSpecAttrs', 'mobility_scan')
        d_params['mobility_rt'] =     self.get_raw_data_attr_data( scan_index, 'MassSpecAttrs', 'mobility_rt')
        d_params['Aterm'] =  self.get_raw_data_attr_data( scan_index, 'MassSpecAttrs', 'Aterm')
        d_params['Bterm'] =  self.get_raw_data_attr_data( scan_index, 'MassSpecAttrs', 'Bterm')
        d_params['Cterm'] = self.get_raw_data_attr_data( scan_index, 'MassSpecAttrs', 'Cterm')
        d_params['baseline_noise'] = self.get_raw_data_attr_data( scan_index, 'MassSpecAttrs', 'baseline_noise')
        d_params['baseline_noise_std'] = self.get_raw_data_attr_data( scan_index, 'MassSpecAttrs', 'baseline_noise_std')
        
        d_params['analyzer'] = self.get_high_level_attr_data('analyzer')
        d_params['instrument_label'] = self.get_high_level_attr_data('instrument_label')
        d_params['sample_name'] = self.get_high_level_attr_data('sample_name')

        return d_params
