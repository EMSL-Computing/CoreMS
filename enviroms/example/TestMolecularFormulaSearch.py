__author__ = "Yuri E. Corilo"
__date__ = "Jul 25, 2019"

import time, sys, os
sys.path.append(".")
from enviroms.emsl.yec.encapsulation.constant.Constants import Atoms
from enviroms.emsl.yec.molecular_id.calc.MolecularLookupTable import  MolecularCombinations
from enviroms.emsl.yec.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaLookupTableSetting
from enviroms.emsl.yec.mass_spectrum.input.TextMassList import Read_MassList
from enviroms.example.TestMolecularTableLookup import create_lookup_table
from enviroms.emsl.yec.encapsulation.constant.Constants import Labels

def creat_mass_spectrum():

    directory = os.path.join(os.getcwd(), "data/")
    
    file_name = "20190616_WK_ESFA_0pt2mgml_ESI_Neg_1pt4sFID_000001.ascii"
    
    file_location = directory + file_name
   
    #polariy need to be set or read from the file
    polariy = -1

    #load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    mass_list_reader = Read_MassList(file_location, polariy, delimiter="  ")

    mass_spectrum_obj  = mass_list_reader.get_mass_spectrum(auto_process=True)
    
    return mass_spectrum_obj




if __name__ == "__main__":
    
    
    precision = 1 #ppm
    time0 = time.time()
    
    dict_molecular_lookup_table = create_lookup_table()
    print(dict_molecular_lookup_table.keys())
    time1 = time.time()
    
    print("create the molecular lookup table took %.2f seconds", time1-time0)
    
    time2 = time.time()
    
    mass_spectrum_obj = creat_mass_spectrum()

    print("Loading the data took %.2f seconds", time2-time1)
    
    #number of peaks
    print(len(mass_spectrum_obj))
    
    ion_type = Labels.de_protonated_ion
    
    classes = dict_molecular_lookup_table.keys()
    
    min_abundance = mass_spectrum_obj.min_abundance
    
    for mspeak in mass_spectrum_obj:
        
        nominal_mz = mspeak.nominal_mz_exp
        
        for classe in classes:

            possible_formulas = dict_molecular_lookup_table.get(classe).get(ion_type).get(nominal_mz)
            
            if possible_formulas:
                
                for possible_formula in possible_formulas:
                    
                    if possible_formula:
                        
                        error = possible_formula._calc_assigment_mass_error(mspeak.mz_exp)
                        
                        if  -1*precision <= error <= precision:
                            
                            isotopologues = possible_formula.isotopologues(min_abundance, mspeak.abundance)
                            
                            time4 = time.time()
                            i = 0  
                            for isotopologue_formula in isotopologues:
                                i += 1
                                #move this outside to impove preformace
                                indexes_to_look = mass_spectrum_obj.get_nominal_mass_indexes(isotopologue_formula.mz_nominal_theo)
                                if indexes_to_look:
                                    
                                    for mspeaks_iso in mass_spectrum_obj[indexes_to_look[0]:indexes_to_look[-1]]:
                                        
                                        error = isotopologue_formula._calc_assigment_mass_error(mspeaks_iso.mz_exp)    
                                        
                                        if  -1*precision <= error <= precision:
                                                pass
                                                #print('YES')         
                                
                            
                            time5 = time.time()
                                
                            #print("Loading iso possible formulas took %.2f seconds with %i possible", time5-time4, i)
                                
                            #print(possible_formula.mz_theor, error, possible_formula.to_string)
        
    time3 = time.time()

    print("Loading possible formulas took %.2f seconds", time3-time2)
    pass#print(mspeak)

