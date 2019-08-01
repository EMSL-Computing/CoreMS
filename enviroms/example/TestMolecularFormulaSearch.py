__author__ = "Yuri E. Corilo"
__date__ = "Jul 25, 2019"

import sys, time, os 
sys.path.append(".")
from enviroms.emsl.yec.transient.input.BrukerSolarix import ReadBrukerSolarix
from enviroms.emsl.yec.molecular_id.calc.MolecularFormulaSearch import SearchMolecularFormulas
from enviroms.emsl.yec.mass_spectrum.input.TextMassList import Read_MassList
from enviroms.emsl.yec.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaLookupTableSettings


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

    mass_spectrum = creat_mass_spectrum(file_location)
    time1 = time.time()
    SearchMolecularFormulas().runworker_mass_spectrum(mass_spectrum)

    print('searching molecular formulas took %i seconds' % (time.time() - time1))

    i = 0
    j = 0
    f = open("20190205_WK_SRFA_opt_000001_resolving_power_1p4sec_12T.txt", "w")
    for mspeak in mass_spectrum:

        if mspeak.is_assigned:
            i += 1
            #print(mspeak.mz_exp, len(mspeak))
            for formula in mspeak:
                print(formula.to_string, formula._calc_assigment_mass_error(mspeak.mz_exp))
                # need to change the calculation of error inside the formula
                if formula.is_isotopologue:
                    pass
                    print(formula.to_string, formula._calc_assigment_mass_error(mspeak.mz_exp))

        else:
            j += 1
            pass
            #print("No Hit")

    print('%i peaks assigned and %i peaks not assigned' % (i, j))
    f.close()
