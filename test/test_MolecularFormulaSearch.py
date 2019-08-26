__author__ = "Yuri E. Corilo"
__date__ = "Jul 25, 2019"

import os
import sys
import time
sys.path.append(".")

from enviroms.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaLookupTableSettings
from enviroms.mass_spectrum.input.TextMassList import Read_MassList
from enviroms.molecular_id.calc.MolecularFormulaSearch import SearchMolecularFormulas
from enviroms.molecular_id.calc.ClusterFilter import ClusteringFilter
from enviroms.transient.input.BrukerSolarix import ReadBrukerSolarix

def creat_mass_spectrum(file_location):

    bruker_reader = ReadBrukerSolarix(file_location)

    bruker_transient = bruker_reader.get_transient()

    mass_spectrum_obj = bruker_transient.get_mass_spectrum(
        plot_result=False, auto_process=True)

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
    kendrick_base =  {'C':1,'H':2,'O':1}   
    
    mass_spectrum = creat_mass_spectrum(file_location)
    mass_spectrum.change_kendrick_base_all_mspeaks(kendrick_base)
    
    ClusteringFilter().filter_kendrick(mass_spectrum)
    
    time1 = time.time()
    
    SearchMolecularFormulas().run_worker_mass_spectrum(mass_spectrum)
    mass_spectrum.reset_indexes()

    print('searching molecular formulas took %i seconds' % (time.time() - time1))
    
    i = 0
    j = 0
    error = list()
    mass = list()
    abundance = list()
    
    for mspeak in mass_spectrum:

        if mspeak.is_assigned:
            i += 1
                
            for mformula in mspeak:
                mass.append(mspeak.mz_exp)
                error.append(mformula._calc_assigment_mass_error(mspeak.mz_exp))
                abundance.append(mspeak.abundance)
        else:
            j += 1
            pass
        
    
    from matplotlib import pylab
    pylab.plot(mass_spectrum.mz_exp, mass_spectrum.abundance, color='g') 
    pylab.plot(mass, abundance, "o")  
    pylab.show()  
    pylab.plot(mass, error, "o")  
    pylab.show()  
    print('%i peaks assigned and %i peaks not assigned' % (i, j))
    #ClusteringFilter().filter_mass_error(mass_spectrum)
    