__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"
from copy import deepcopy
from datetime import datetime
from pathlib import Path

from numpy import genfromtxt, fromstring, dtype, fromfile, frombuffer, float64, float32
import pandas as pd
from s3path import S3Path
from xml.dom import minidom

from corems.transient.factory.TransientClasses import Transient
from corems.encapsulation.factory.parameters import default_parameters

class ReadBrukerSolarix(object):
    """    A class used to Read a single Transient from Bruker's FT-MS acquisition station (fid, or ser)
        
    Parameters
    ----------
    d_directory_location : str
        the full path of the .d folder
    
    Attributes
    --------
    d_directory_location : str
        the full path of the .d folder
    file_location : str
        the full path of the .d folder
    parameter_filename_location : str
        the full path of the apexAcquisition.method file
    transient_data_path : str
        the full path of the fid or ser file
    scan_attr : str
        the full path of the scan.xml file
    
    
    Methods
    -------
    * get_transient().
        Read the data and settings returning a Transient class  
    * get_scan_attr().
        Read the scan retention times, TIC values and scan indices.
    * locate_file(folder, type_file_name).
        Find the full path of a specific file within the acquisition .d folder or subfolders
    * parse_parameters(parameters_filename).
        Open the given file and retrieve all parameters from apexAcquisition.method
    * fix_freq_limits(d_parameters).
        Read and set the correct frequency limits for the spectrum
    * get_excite_sweep_range(filename).
        Determine excitation sweep range from ExciteSweep file
    
    """
    
    def __enter__(self ):
            
            return self.get_transient()
    
    def __exit__(self, exc_type, exc_val, exc_tb):
            
        return False

    def __init__(self, d_directory_location):
        
        if isinstance(d_directory_location, str):
            d_directory_location = Path(d_directory_location)
        
        if not d_directory_location.exists():
            raise FileNotFoundError("File does not exist: " + str(d_directory_location))

        self.d_directory_location = d_directory_location
        
        self.file_location = d_directory_location
        
        try:

            self.parameter_filename_location = self.locate_file(
                d_directory_location, "apexAcquisition.method"
            )
            self.transient_data_path = d_directory_location / "fid"
            
            if not self.transient_data_path.exists():

                self.transient_data_path = d_directory_location / "ser"

                if not self.transient_data_path.exists():
                    
                    raise FileNotFoundError("Could not locate transient data")

                else:
                    # get scan attributes
                    self.scan_attr = d_directory_location / "scan.xml"

        except:
            
            raise FileExistsError(
                "%s does not seem to be a valid Solarix Mass Spectrum"
                % (d_directory_location)
            )

    def get_scan_attr(self):
        """ Function to get the scan retention times, TIC values and scan indices. 

        Gets information from scan.xml file in the bruker .d folder.
        Note this file is only present in some .d format - e.g. for imaging mode data, it is not present.
        
        Returns
        -------
        dict_scan_rt_tic : dict
            a dictionary with scan number as key and rt and tic as values
        """
    
        from bs4 import BeautifulSoup
        
        try: 
            soup = BeautifulSoup(self.scan_attr.open(),'xml')
        except:
            raise FileNotFoundError("Dataset does not appear to contain a 'scan.xml' file or it is misformated")

        list_rt = [float(rt.text) for rt in soup.find_all('minutes')]
        list_tic = [float(tic.text) for tic in soup.find_all('tic')]
        list_scan = [int(scan.text) for scan in soup.find_all('count')]

        dict_scan_rt_tic = dict(zip(list_scan, zip(list_rt, list_tic)))
        
        return dict_scan_rt_tic
       
        
    def get_transient(self, scan_number=1):
        """ Function to get the transient data and parameters from a Bruker Solarix .d folder.
        
        Parameters
        ----------
        scan_number : int
            the scan number to be read. Default is 1.
        
        Returns
        -------
        Transient
            a transient object
        """

        file_d_params = self.parse_parameters(self.parameter_filename_location)

        self.fix_freq_limits(file_d_params)

        from sys import platform
        
        if platform == "win32":
            # Windows...
            dt = dtype("l")
        else:
            dt = dtype("i")

        # get rt, scan, and tic from scan.xml file, otherwise  using 0 defaults values 
        
        output_parameters = deepcopy(default_parameters(self.d_directory_location))

        if self.transient_data_path.name == 'ser':
            
            if self.scan_attr.exists():
                
                dict_scan_rt_tic = self.get_scan_attr()

                output_parameters["scan_number"] = scan_number

                output_parameters["rt"] = dict_scan_rt_tic.get(scan_number)[0]

                output_parameters["tic"] = dict_scan_rt_tic.get(scan_number)[1]
        
        output_parameters["analyzer"] = "ICR"

        output_parameters["label"] = "Bruker_Frequency"

        output_parameters["Aterm"] = float(file_d_params.get("ML1"))

        output_parameters["Bterm"] = float(file_d_params.get("ML2"))

        output_parameters["Cterm"] = float(file_d_params.get("ML3"))

        output_parameters["exc_high_freq"] = float(file_d_params.get("EXC_Freq_High"))

        output_parameters["exc_low_freq"] = float(file_d_params.get("EXC_Freq_Low"))
        try:
            output_parameters["qpd_enabled"] = float(file_d_params.get("QPD_Enabled"))
        except TypeError: # for older datasets which dont have this variable
            output_parameters["qpd_enabled"] = 0

        output_parameters["mw_low"] = float(file_d_params.get("MW_low"))

        output_parameters["mw_high"] = float(file_d_params.get("MW_high"))

        output_parameters["bandwidth"] = float(file_d_params.get("SW_h"))

        output_parameters["number_data_points"] = int(file_d_params.get("TD"))

        output_parameters["polarity"] = str(file_d_params.get("Polarity"))

        output_parameters["acquisition_time"] = file_d_params.get("acquisition_time")

        data_points = int(file_d_params.get("TD"))

        scan = output_parameters["scan_number"]
        from io import BytesIO
        if self.transient_data_path.name == 'ser':
            
            if isinstance(self.transient_data_path, S3Path):
                databin = BytesIO(self.transient_data_path.open('rb').read())
            
            else:
                databin = self.transient_data_path.open('rb')
               
            databin.seek((scan-1)*4*data_points)
            #read scan data and parse to 32int struct
            data = frombuffer(databin.read(4*data_points), dtype=dt)
        
        else:
            
            if isinstance(self.transient_data_path, S3Path):
                data = frombuffer(self.transient_data_path.open('rb').read(), dtype=dt)
            else:
                data = fromfile(self.transient_data_path, dtype=dt)
        
        return Transient(data, output_parameters)

    #    for key, values in default_parameters.items():
    #        print(key, values)
    def fix_freq_limits(self, d_parameters):
        """ Function to read and set the correct frequency limits for the spectrum
        
        Notes
        --------
        This is using the excitation limits from the apexAcquisition.method file,
        which may not match the intended detection limits in edge cases. 
        In default acquisitions, excitation and detection are the same. 
        But, they may not be in some cases with selective excitation, custom excite waveforms, or in 2DMS applications.
        
        Parameters
        ----------
        d_parameters : dict
            a dictionary with the parameters from the apexAcquisition.method file
        """

        highfreq = float(d_parameters.get("EXC_Freq_High"))

        lowfreq = float(d_parameters.get("EXC_Freq_Low"))

        # CR for compatibility with Apex format as there is no EXciteSweep file
        if not highfreq and lowfreq:

            excitation_sweep_filelocation = self.locate_file(
                self.d_directory_location, "ExciteSweep"
            )
            lowfreq, highfreq = self.get_excite_sweep_range(
                excitation_sweep_filelocation
            )
            d_parameters["EXC_Freq_High"] = highfreq
            d_parameters["EXC_Freq_Low"] = lowfreq

    @staticmethod
    def get_excite_sweep_range(filename):
        """ Function to determine excitation sweep range from ExciteSweep file

        This looks at the first and last rows of the ExciteSweep file to determine the excitation frequency range.
        Note that this assumes the excitation sweep was linear and the first and last rows are the lowest and highest frequencies.
        This is presumably always true, but again may be incorrect for edge cases with custom excitation waveforms.

        Parameters
        ----------
        filename : str
            the full path to the ExciteSweep file
        
        """
        ExciteSweep_lines = genfromtxt(filename, comments="*", delimiter="\n")
        # CR ready if we need the full array
        highfreq = fromstring(ExciteSweep_lines[0])
        lowfreq = fromstring(ExciteSweep_lines[-1])

        return lowfreq[0], highfreq[0]

    @staticmethod
    def locate_file(folder, type_file_name='apexAcquisition.method'):
        """ Function to locate a file in a folder

        Find the full path of a specific file within the acquisition .d folder or subfolders

        Parameters
        ----------
        folder : str
            the full path to the folder
        type_file_name : str
            the name of the file to be located
            Expected options: ExciteSweep or apexAcquisition.method

        Returns
        -------
        str
            the full path to the file

        Notes
        -----
        adapted from code from SPIKE library, https://github.com/spike-project/spike
                
        """
        
        from pathlib import Path
               
        #directory_location = folder.glob( '**/*apexAcquisition.method')
        directory_location = folder.glob( '**/*' + type_file_name)
        result = list(directory_location)
        if len(result) > 1:

            raise Exception(
                "You have more than 1 %s file in the %s folder, using the first one"
                % (type_file_name, folder)
            )

        elif len(result) == 0:

            raise Exception(
                "You don't have any %s file in the  %s folder, please double check the path"
                % (type_file_name, folder)
            )

        return result[0]

    @staticmethod
    def parse_parameters(parameters_filename):
        """ Function to parse the parameters from apexAcquisition.method file

        Open the given file and retrieve all parameters from apexAcquisition.method
            None is written when no value for value is found
            
            structure : <param name = "AMS_ActiveExclusion"><value>0</value></param>
        
        Parameters
        ----------
        parameters_filename : str
            the full path to the apexAcquisition.method file
        
        Returns
        -------
        dict
            a dictionary with the parameters and values
        
        Notes
        -----
        Adapted from code from SPIKE library, https://github.com/spike-project/spike.
        Code may not handle all possible parameters, but should be sufficient for most common use cases
        """
        
        #TODO: change to beautiful soup xml parsing
        
        
        xmldoc = minidom.parse(parameters_filename.open())

        x = xmldoc.documentElement
        parameter_dict = {}
        children = x.childNodes
        for child in children:
            # print( child.node)
            if child.nodeName == 'methodmetadata':
                sections = child.childNodes
                for section in sections:
                    for element in section.childNodes:
                        if element.nodeName == "date":
                        #if element.nodeName == "primarykey":
                            
                            date_time_str = (element.childNodes[0].nodeValue)
                            #parameter_dict["acquisition_time"] = pd.to_datetime(date_time_str, infer_datetime_format=True).to_pydatetime()
                            parameter_dict["acquisition_time"] = datetime.strptime(date_time_str, "%b_%d_%Y %H:%M:%S.%f")
                            
            
            if child.nodeName == "reportinfo":
                sections = child.childNodes
                for section in sections:
                    if section.nodeName == "section":
                        if section.getAttribute("title") == "Main":
                            for element in section.childNodes:
                                if element.nodeName == "section":
                                    if element.getAttribute("title") == "Polarity":
                                        if (
                                            str(
                                                element.childNodes[1].getAttribute(
                                                    "value"
                                                )
                                            )
                                            == "Negative"
                                        ):
                                            parameter_dict["Polarity"] = -1
                                        else:
                                            parameter_dict["Polarity"] = 1

            if child.nodeName == "paramlist":
                params = child.childNodes
                for param in params:
                    # print( param.nodeName)
                    if param.nodeName == "param":
                        paramenter_label = str(param.getAttribute("name"))
                        for element in param.childNodes:
                            if element.nodeName == "value":
                                try:
                                    parameter_value = str(element.firstChild.toxml())
                                    # print v
                                except:
                                    parameter_value = None

                            parameter_dict[paramenter_label] = parameter_value

        return parameter_dict


    def parse_sqlite(self, sqlite_filename="chromatography-data.sqlite"):
        """
        
        """
        import sqlite3

        def read_sqlite_file(file_path, table_name):
            """
            Read data from a SQLite database file and return it as a list of tuples

            Parameters
            ----------
            file_path : str
                the full path to the SQLite database file
            table_name : str
                the name of the table to be read
            
            Returns
            -------
            list
                a list of tuples with the data from the table
            """
             # Connect to the SQLite database file
            conn = sqlite3.connect(file_path)
            cursor = conn.cursor()

            # Execute a query to select data from a table (replace 'table_name' with your table's name)
            query = f"SELECT * FROM {table_name}"
            cursor.execute(query)

            # Fetch all rows from the result set
            rows = cursor.fetchall()
            stream = []
            # Print or process the fetched rows
            for row in rows:
                stream.append(row)
                # print(row)  # Print each row, you can also process it differently

            # Close the cursor and the connection
            cursor.close()
            conn.close()
            return stream
        
        def parse_binary(binary, type):
            """
            Parse binary data from the sqlite data streams
            """
            if type == "double":
                data = frombuffer(binary, dtype=float64)
            elif type == "float":
                data = frombuffer(binary, dtype=float32)
            return data
        
        sqlite_filelocation = self.locate_file(
                self.d_directory_location, sqlite_filename
            )
        table_name = "TraceSources"
        trace_sources = read_sqlite_file(sqlite_filelocation, table_name)
        table_name = "TraceChunks"
        trace_chunks = read_sqlite_file(sqlite_filelocation, table_name)
        times = []
        values = []
        trace_type = {}


        for index, source in enumerate(trace_sources):
            trace_id = source[0]
            trace_type[source[1]] = {"times": [], "values": []}
            for index, chunk in enumerate(trace_chunks):
                id = chunk[0]
                times = parse_binary(chunk[1], "double")
                values = parse_binary(chunk[2], "float")
                for time, value in zip(times, values):
                    if source[0] == id:
                        trace_type[source[1]]["times"].append(time)
                        trace_type[source[1]]["values"].append(value)

        return trace_type