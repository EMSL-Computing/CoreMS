'''
Created on Jun 19, 2019

@author: eber373
'''
from __future__ import print_function

from emsl.yec.structure.input.BrukerSolarix import ReadBrukerSolarix
from emsl.yec.structure.farm.TD_Class import Transient
        
filelocation = "C:\\Users\\eber373\\Desktop\\20190616_WK_ESFA_0pt2mgml_ESI_Neg_0pt8sFID_000001.d\\"    
filelocation = "C:\\Users\\eber373\\Desktop\\20190205_WK_SRFA_opt_000001.d\\"


data, d_params = ReadBrukerSolarix(filelocation).read_file()
bruker_transient = Transient(data, d_params)
#bruker_transient.plot_transient(None)
bruker_transient.set_frequency_domain()
bruker_transient.plot_frequency_domain()
bruker_transient.plot_mz_domain()
