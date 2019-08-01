
import sys, os, time
sys.path.append(".")
from enviroms.emsl.yec.transient.input.BrukerSolarix import ReadBrukerSolarix
from enviroms.emsl.yec.molecular_id.calc.FindOxigenPeaks import FindOxygenPeaks
from enviroms.emsl.yec.mass_spectrum.input.TextMassList import Read_MassList
from enviroms.emsl.yec.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaLookupTableSettings

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()

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
    find_formula_thread.start()
    
    #while find_formula_thread.is_alive():
    #    for i in range(100):
    #        printProgressBar(i, 100)
    find_formula_thread.join()

    mspeaks_results = find_formula_thread.get_list_found_peaks()
    mspeaks_results = sorted(mspeaks_results, key=lambda mp: mp.mz_exp)
    for mspeak in mspeaks_results:
        
        print('found',mspeak.molecular_formula_lowest_error.to_string, mspeak.mz_exp)
    
    print('searching molecular formulas took %i seconds' % (time.time() - time1))