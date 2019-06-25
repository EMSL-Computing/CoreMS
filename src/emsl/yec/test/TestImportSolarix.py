'''
Created on Jun 19, 2019

@author: eber373
'''


import pickle

from emsl.yec.structure.factory.TransientClasses import Transient
from emsl.yec.structure.input.BrukerSolarix import ReadBrukerSolarix


#from emsl.yec.structure.input.MidasDatFile import ReadMidasDatFile
filelocation = "C:\\Users\\eber373\\Desktop\\20190616_WK_ESFA_0pt2mgml_ESI_Neg_0pt8sFID_000001.d\\"    
filelocation = "C:\\Users\\eber373\\Desktop\\20190205_WK_SRFA_opt_000001.d\\"

apodization_method = 'Hanning'
number_of_truncations = 0
number_of_zero_fills = 1


data, d_params = ReadBrukerSolarix(filelocation).read_file()

bruker_transient = Transient(data, d_params)
bruker_transient.set_processing_parameter(apodization_method, number_of_truncations, number_of_zero_fills)

mass_spec = bruker_transient.generate_mass_spec()
mass_spec.set_mspeaks()
mass_spec.plot_mz_domain_profile()
#mass_spec.set_processing_parameter()

with open('test.pkl', 'wb') as file:
    pickle.dump(bruker_transient, file, protocol=pickle.HIGHEST_PROTOCOL)


#transient = pickle.load( open( 'test.pkl', "rb" ) )
#do_something

'''

filelocation = "C:\\Users\\eber373\\Desktop\\20190616_WK_ESFA_0pt2mgml_ESI_Neg_0pt8sFID_000001.d\\"

data, d_params = ReadMidasDatFile(filelocation).read_file()

dat_transient = Transient(data, d_params)

dat_transient.set_frequency_domain()

dat_transient.plot_frequency_domain()

dat_transient.plot_mz_domain()
'''

