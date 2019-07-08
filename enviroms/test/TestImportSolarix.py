import pickle, sys

__author__ = "Yuri E. Corilo"
__date__ = "Jun 19, 2019"
sys.path.append(".")
from enviroms.emsl.yec.input.BrukerSolarix import ReadBrukerSolarix
if __name__ == '__main__':
    #from enviroms.emsl.yec.structure.input.MidasDatFile import ReadMidasDatFile
    filelocation = "C:\\Users\\eber373\\Desktop\\data\\20190616_WK_ESFA_0pt2mgml_ESI_Neg_0pt8sFID_000001.d\\"    
    filelocation = "C:\\Users\\eber373\\Desktop\\data\\20190205_WK_SRFA_opt_000001.d\\"
    
    apodization_method = 'Hanning'
    number_of_truncations = 0
    number_of_zero_fills = 1
    
    bruker_transient = ReadBrukerSolarix(filelocation)
    bruker_transient.set_processing_parameter(apodization_method, number_of_truncations, number_of_zero_fills)
    mass_spec = bruker_transient.generate_mass_spec(plot_result=False)
    
    mass_spec.cal_noise_treshould()
    mass_spec.find_peaks()
    
    print(mass_spec.mspeaks[0].exp_mz, mass_spec.mspeaks[-1].exp_mz)
    
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

