
import sys, os, time
sys.path.append(".")
from enviroms.emsl.yec.transient.input.BrukerSolarix import ReadBrukerSolarix
from enviroms.emsl.yec.molecular_id.calc.FindOxigenPeaks import FindOxygenPeaks
from enviroms.emsl.yec.mass_spectrum.input.TextMassList import Read_MassList
from enviroms.emsl.yec.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaLookupTableSettings

def creat_mass_spectrum(file_location):

    bruker_reader = ReadBrukerSolarix(file_location)

    bruker_transient = bruker_reader.get_transient()

    mass_spectrum_obj = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=True)
    #mass_spectrum_obj.plot_mz_domain_profile_and_noise_threshold()
    
    # polariy need to be set if reading a text file
    #polariy = -1
    # load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    #mass_list_reader = Read_MassList(file_location, polariy, delimiter="  ")
    #mass_spectrum_obj  = mass_list_reader.get_mass_spectrum(auto_process=True)

    return mass_spectrum_obj

if __name__ == "__main__":

    # MoleculaLookupTableSettings and MoleculaSearchSettings at
    # enviroms\emsl\yec\encapsulation\settings\molecular_id\MolecularIDSettings.py
    # for changing settings of the lookup table and searching algorithms

    directory = os.path.join(os.getcwd(), "data/")

    file_name = os.path.normcase("20190205_WK_SRFA_opt_000001.d/")

    #file_name = "20190616_WK_ESFA_0pt2mgml_ESI_Neg_1pt4sFID_000001.ascii"

    file_location = directory + file_name

    mass_spectrum = creat_mass_spectrum(file_location)
    
    time1 = time.time()

    find_formula_thread = FindOxygenPeaks(mass_spectrum)
    find_formula_thread.run()
    
    #while find_formula_thread.is_alive():
    #    for i in range(100):
    #        printProgressBar(i, 100)
    #find_formula_thread.join()

    mspeaks_results = find_formula_thread.get_list_found_peaks()
    
    mspeaks_results = sorted(mspeaks_results, key=lambda mp: mp.mz_exp)
    
    error = list()
    mass = list()
    abundance = list()
    for mspeak in mspeaks_results:
        
        for molecular_formula in mspeak:
            
            mass.append(mspeak.mz_exp)
            error.append(molecular_formula._calc_assigment_mass_error(mspeak.mz_exp))
            abundance.append(mspeak.abundance)

    #mass_spectrum.mz_exp_centroide
    #mass_spectrum.abundance_centroid
    from matplotlib import pylab
    pylab.plot(mass_spectrum.mz_exp, mass_spectrum.abundance) 
    pylab.plot(mass, abundance, "o")  
    pylab.show()  
    print('searching molecular formulas took %i seconds' % (time.time() - time1))
    pylab.plot(mass, error, "o")  
    pylab.show()  
    