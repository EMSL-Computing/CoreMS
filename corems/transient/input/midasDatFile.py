from os import path

from numpy import dtype, fromfile


__author__ = "Yuri E. Corilo"
__date__ = "Jun 19, 2019"


class ReadMidasDatFile:
    """[Not Implemented] Reads MIDAS .dat files (binary transient data)

    This class will read .dat binary format transient data, e.g. midas format from Predator or Thermo datastations
    This code is not yet implemented and is not fully functional.

    Parameters
    ----------
    filename_path : str
        The path to the .dat file

    Attributes
    ----------
    filename_path : str
        The path to the .dat file
    d_params : dict
        A dictionary with the parameters of the .dat file
    transient_data : numpy.ndarray
        The transient data

    Methods
    -------
    * read_file().
        Reads the .dat file and returns the transient data and the parameters
    * get_transient_data(data_file, d_params).
        Reads the transient data from the .dat file
    * parse_parameter(f).
        Parses the parameters from the .dat file

    Raises
    ------
    NotImplementedError
        This class is not yet implemented.

    """

    def __init__(self, filename_path):
        raise NotImplementedError(
            "This class is not yet implemented, if you want to use it please contact the author at corilo@pnnl.gov or feel free to implement it"
        )
        if not path.isfile(filename_path):
            raise Exception("File does not exist: " + filename_path)

        self.filename_path = filename_path

    def read_file(self):
        """Reads the .dat file and returns the transient data and the parameters

        Returns
        -------
        transient_data : numpy.ndarray
            The transient data
        d_params : dict
            A dictionary with the parameters of the .dat file

        """
        data_file = open(self.filename_path, "rb")

        # modo_de_ions = "POSITIVE ION MODE"
        d_params = self.parse_parameters(self.parameter_filename_location)

        transient_data = self.get_transient_data(data_file, d_params, d_params)

        return transient_data, d_params

    def get_transient_data(self, data_file, d_params):
        """Reads the transient data from the .dat file

        Parameters
        ----------
        data_file : file
            The .dat file
        d_params : dict
            A dictionary with the parameters of the .dat file

        Returns
        -------
        myarray : numpy.ndarray
            The transient data
        """

        # dt = np.dtype('<f')
        if d_params.get("storage_type").split()[0] == "int":
            dt = dtype("i2")

        else:
            dt = dtype("<f")
        # dt = np.dtype(int)

        myarray = fromfile(data_file, dtype=dt)

        data_file.close()

        if d_params.get("storage_type").split()[0] == "int":
            return myarray * d_params.get("VoltageScale")

        else:
            return myarray

    def parse_parameter(self, f):
        """Parses the parameters from the .dat file

        Parameters
        ----------
        f : file
            The .dat file

        Returns
        -------
        output_parameters : dict
            A dictionary with the parameters of the .dat file
        """

        output_parameters = {}
        output_parameters["filename_path"] = self.d_directory_location

        line = f.readline()

        while line != "Data:\n":
            if line[0:8] == "highfreq":
                final_frequency = float(line.split(":")[1])
                output_parameters["exc_high_freq"] = final_frequency

            elif line[0:7] == "lowfreq":
                initial_frequency = float(line.split(":")[1])
                output_parameters["exc_low_freq"] = initial_frequency

            elif line[0:9] == "sweeprate":
                sweeprate = float(line.split(":")[1])

                output_parameters["sweeprate"] = sweeprate

            elif line[0:13] == "Source Coeff0":
                Acoef = float(line.split(":")[1])
                output_parameters["Aterm"] = Acoef
                # print f.readline()
            elif line[0:13] == "Source Coeff1":
                output_parameters["Bterm"] = "Bcoef"

            elif line[0:13] == "Voltage Scale":
                voltage_scale = float(line.split(":")[1])
                output_parameters["VoltageScale"] = voltage_scale

            elif line[0:9] == "Bandwidth":
                bandwidth = float(line.split(":")[1])
                output_parameters["bandwidth"] = bandwidth

            elif line[0:11] == "Data Points":
                datapoints = float(line.split(":")[1])
                output_parameters["number_data_points"] = datapoints

            elif line[0:12] == "Storage Type":
                storage_type = line.split(":")[1]
                output_parameters["storage_type"] = storage_type

            elif line[0:12] == "Trap Voltage":
                trap_voltage = float(line.split(":")[1])
                # Bcoef = Bcoef*trap_voltage
                output_parameters["trap_voltage"] = trap_voltage

            line = f.readline()

        return output_parameters
