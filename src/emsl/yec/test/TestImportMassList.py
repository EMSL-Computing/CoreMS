'''
Created on Jun 28, 2019

@author: eber373
'''
from emsl.yec.structure.input.TextMassList import Read_MassList
from emsl.yec.structure.factory.mass_spec.MS_Centroid import MassSpecCentroid

directory = "C:\\Users\\eber373\\Documents\\Desenvolvimento\\Software Projects\\EnviroMS\\res\\"

#filename_ESFA = "20190616_WK_ESFA_0pt2mgml_ESI_Neg_0pt8sFID_000001.d"    
#filename_SRFA = "20190205_WK_SRFA_opt_000001.d"

file_name_ESFA_txt = "20190616_WK_ESFA_0pt2mgml_ESI_Neg_1pt4sFID_000001.ascii"

full_file_path = directory + file_name_ESFA_txt

polariy = -1
 
data, d_params = Read_MassList(full_file_path, polariy).read_file()

mass_spec = MassSpecCentroid(data, d_params)

mass_spec.assign_molecular_formulas()

mass_spec.plot_mz_domain_profile()