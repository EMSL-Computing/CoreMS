'''
Created on Jun 19, 2019

@author: eber373
'''
from __future__ import print_function

from xml.dom import minidom

from emsl.yec.structure.input.parsers.BrukerSolarix import ReadBrukerSolarix


def read_param( parameters_filename):
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
        #print(children)
        for child in children:
            if (child.nodeName == 'paramlist'):
                params = child.childNodes
                for param in params:
                    if (param.nodeName == 'param'):
                        paramenter_label = str(param.getAttribute('name'))
                        for element in param.childNodes:
                            if element.nodeName == "value":
                                try:
                                    parameter_value = element.firstChild.toxml()
                                    #print v
                                except: 
                                    parameter_value = None
                        
                        parameter_dict[paramenter_label] = parameter_value
        
        return parameter_dict
        xmldoc.close()
        
filelocation = "C:\\Users\\eber373\\Desktop\\20190616_WK_ESFA_0pt2mgml_ESI_Neg_0pt8sFID_000001.d\\"    
filelocation = "C:\\Users\\eber373\\Desktop\\20190205_WK_SRFA_opt_000001.d\\"
#x = read_param(filelocation)
#print('XXXXXX')
#print(x.get("SW_h"))

read_instance = ReadBrukerSolarix(filelocation)
read_instance.read_file()