__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"
from copy import deepcopy
from pathlib import Path

from numpy import genfromtxt, fromstring, dtype, fromfile, frombuffer
from s3path import S3Path
from xml.dom import minidom


from corems.transient.factory.TransientClasses import Transient
from corems.encapsulation.factory.parameters import default_parameters



class ReadBrukerSolarix(object):
    
    """
    A class used to Read a single Transient from Bruker's FT-MS acquisition station (fid, or ser)
        
    Parameters
    ----------
    d_directory_location : str
        the full path of the .d folder
    
    Methods
    -------
    get_transient()
        Read the data and settings returning a Transient class  
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
    
        from bs4 import BeautifulSoup
        
        soup = BeautifulSoup(self.scan_attr.open(),'xml')

        list_rt = [float(rt.text) for rt in soup.find_all('minutes')]
        list_tic = [float(tic.text) for tic in soup.find_all('tic')]
        list_scan = [int(scan.text) for scan in soup.find_all('count')]

        dict_scan_rt_tic = dict(zip(list_scan, zip(list_rt, list_tic)))
        
        return dict_scan_rt_tic
       
        
    def get_transient(self, scan_number=1):

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

        output_parameters["bandwidth"] = float(file_d_params.get("SW_h"))

        output_parameters["number_data_points"] = int(file_d_params.get("TD"))

        output_parameters["polarity"] = str(file_d_params.get("Polarity"))

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

    """
        for key, values in default_parameters.items():
            print(key, values)
    """
    def fix_freq_limits(self, d_parameters):

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
        """
        Function that returns the lower and higher frequency of the pulse generator
        """
        ExciteSweep_lines = genfromtxt(filename, comments="*", delimiter="\n")
        # CR ready if we need the full array
        highfreq = fromstring(ExciteSweep_lines[0])
        lowfreq = fromstring(ExciteSweep_lines[-1])

        return lowfreq[0], highfreq[0]

    @staticmethod
    def locate_file(folder, type_file_name):
        
        from pathlib import Path
        """
            type_of_file = ExciteSweep or apexAcquisition.method
            From the given folder this function return the absolute path to the ExciteSweep file, or the apexAcquisition.method file
            It should always be in a subfolder 
        
            ''' this code is a adaptation/automation of spike library
                https://bitbucket.org/delsuc/spike/
                from Marc Delsuc
            '''
        """
        
        directory_location = folder.glob( '**/*apexAcquisition.method')
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
        """ 
            TODO: change to beautiful soup xml parsing
            Open the given file and retrieve all parameters from apexAcquisition.method
            None is written when no value for value is found

            structure : <param name = "AMS_ActiveExclusion"><value>0</value></param>

            read_param returns  values in a dictionnary
            xml should be extinct, just a random option

            ''' this code is a adaptation/automation of spike library
                https://bitbucket.org/delsuc/spike/
                from Marc Delsuc
            '''
        """
        xmldoc = minidom.parse(parameters_filename.open())

        x = xmldoc.documentElement
        parameter_dict = {}
        children = x.childNodes
        for child in children:
            # print( child.node)
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

