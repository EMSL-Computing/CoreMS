
from glob import glob
from os import path
from xml.dom import minidom
from numpy import genfromtxt, fromstring, dtype, fromfile, linspace


class ReadBrukerSolarix():
    
    def __init__(self, d_directory_location):
        
        self.d_directory_location = d_directory_location
        
        if not path.exists(d_directory_location):
            raise Exception("File does not exist: "+ d_directory_location)
        
        try:
            
            self.parameter_filename_location = self.locate_file(d_directory_location, 'apexAcquisition.method')
            self.transient_data_filename = path.join(d_directory_location, "fid")
            
            if not path.isfile(self.transient_data_filename):
                
                raise Exception("Could not locate transient data")
            
        except:
            
            raise Exception('%s does not seem to be a valid Solarix spectrum'%(d_directory_location,))
            
    def read_file(self):
        
        d_parameters = self.parse_parameters(self.parameter_filename_location)
        
        number_data_points = int(d_parameters.get("TD"))
        
        bandwidth = float(d_parameters.get("SW_h"))
        
        highfreq = float(d_parameters.get("EXC_Freq_High"))
        
        lowfreq = float(d_parameters.get("EXC_Freq_Low"))
        
        #CR for compatibility with Apex format as there is no EXciteSweep file
        if not highfreq and lowfreq:
            
            excitation_sweep_filelocation = self.locate_file(self.d_directory_location, "ExciteSweep")
            lowfreq, highfreq = self.get_excite_sweep_range(excitation_sweep_filelocation)
            print (highfreq, lowfreq) 
        
        Atherm = float(d_parameters.get("ML1"))
        BTherm = float(d_parameters.get("ML2"))
        CTherm = float(d_parameters.get("ML3"))
        #buffer = zeros((number_data_points))
        
        dt = dtype('l')    
        data = fromfile(self.transient_data_filename, dtype=dt)
        print (len(data))
        #print(number_data_points)
        
        transient_time = (1 / bandwidth) * ((number_data_points) / 2)
        
        
        import matplotlib.pyplot as plt
        time_axis = linspace(0, transient_time, num=len(data))
        plt.plot(time_axis, data)
        plt.show()
        return data

    @staticmethod    
    def get_excite_sweep_range(filename):
        """
        Function that returns the lower and higher frequency of the pulse generator
        """
        ExciteSweep_lines = genfromtxt(filename, comments = "*", delimiter="\n") 
        #CR ready if we need the full array
        highfreq = fromstring(ExciteSweep_lines[0])
        lowfreq = fromstring(ExciteSweep_lines[-1])
        
        return lowfreq[0], highfreq[0]
    
    @staticmethod    
    def locate_file(folder, type_file_name):
        
        """
            type_of_file = ExciteSweep or apexAcquisition.method
            From the given folder this function return the absolute path to the ExciteSweep file, or the apexAcquisition.method file
            It should always be in a subfolder 
        """
        directory_location = glob(path.join(folder,"*",type_file_name))
        
        if len(directory_location)>1:
            
            raise Exception( "You have more than 1 %s file in the %s folder, using the first one" % (type_file_name, folder) )
        
        elif len(directory_location) == 0: 
        
            raise Exception( "You don't have any %s file in the  %s folder, please double check the path" % (type_file_name, folder) )
        
        return directory_location[0]
    
    @staticmethod  
    def parse_parameters(parameters_filename):
        """
            Open the given file and retrieve all parameters from apexAcquisition.method
            NC is written when no value for value is found
            
            structure : <param name = "AMS_ActiveExclusion"><value>0</value></param>
           
            read_param returns  values in a dictionnary
        """
        xmldoc = minidom.parse(parameters_filename)
        
        x = xmldoc.documentElement
        parameter_dict = {}
        children = x.childNodes
        for child in children:
            if (child.nodeName == 'paramlist'):
                params = child.childNodes
                for param in params:
                    if (param.nodeName == 'param'):
                        paramenter_label = str(param.getAttribute('name'))
                        for element in param.childNodes:
                            if element.nodeName == "value":
                                try:
                                    parameter_value = str(element.firstChild.toxml())
                                    #print v
                                except: 
                                    parameter_value = None
                        
                            parameter_dict[paramenter_label] = parameter_value
        
        return parameter_dict
    
    
    