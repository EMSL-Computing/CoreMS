import pickle, sys, os

sys.path.append(".")
from enviroms.emsl.yec.transient.input.BrukerSolarix import ReadBrukerSolarix

__author__ = "Yuri E. Corilo"
__date__ = "Jun 19, 2019"

if __name__ == "__main__":
    
    # from enviroms.emsl.yec.structure.input.MidasDatFile import ReadMidasDatFile
    directory = os.path.join(os.getcwd(), "data/")

    #file_name = os.path.normcase("20190616_WK_ESFA_0pt2mgml_ESI_Neg_0pt8sFID_000001.d/")
    file_name = os.path.normcase("20190205_WK_SRFA_opt_000001.d/")

    file_location = directory + file_name
    
    #setting for signal processing
    apodization_method = "Hanning"
    number_of_truncations = 0
    number_of_zero_fills = 1

    bruker_reader = ReadBrukerSolarix(file_location)

    bruker_transient = bruker_reader.get_transient()

    bruker_transient.set_processing_parameter(
        apodization_method, number_of_truncations, number_of_zero_fills
    )

    mass_spec = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=True)

    #mass_spec.plot_mz_domain_profile_and_noise_threshold()

    #print(mass_spec.mspeaks[0].mz_exp, mass_spec.mspeaks[-1].mz_exp)

    #with open("test.pkl", "wb") as file:
    #    pickle.dump(bruker_transient, file, protocol=pickle.HIGHEST_PROTOCOL)

    # transient = pickle.load( open( 'test.pkl', "rb" ) )
    # do_something
