'''
Created on May 17, 2017

@author: zuiur
'''
import sys
sys.path.append(".")
import csv, time
import itertools
from multiprocessing import freeze_support
from operator import itemgetter
from scipy.stats import norm
from enviroms.emsl.yec.molecular_id.calc.MolecularLookupTable import MolecularCombinations

def create_lookup_table():
    min_dbe = 0

    max_dbe = 50

    use_pah_line_rule = True

    #margin_error needs to be optimized by the data rp and sn
    margin_error = 0.3

    #min_mz,max_mz  needs to be optimized by the data
    min_mz = 100
    max_mz = 1300

    isRadical = True

    isProtonated = True

    ionization_type = "ESI"

    ionCharge = -1

    #this needs to be optimized by the data
    usedAtoms = {'C': (1, 100), 'H': (4, 200), 'O': (0, 20), 'N': (0, 2), 'S': (0, 2)}
    dict_molecular_formulas = MolecularCombinations().runworker(usedAtoms, ionCharge, ionization_type, isProtonated, isRadical, use_pah_line_rule, min_dbe, max_dbe)

    print("dict_nominal_mass_listformulae", len(dict_molecular_formulas))
    print(dict_molecular_formulas.keys())
    print(dict_molecular_formulas.get("O2").keys())
    for molecular_formulas in dict_molecular_formulas.get("O2").get(201):
        print( molecular_formulas.class_label)
        print( molecular_formulas.to_string)
        print( molecular_formulas.mz_theor)
        print(molecular_formulas._cal_isotopologues())
        #print( molecular_formulas.atoms)
        #print( molecular_formulas.ion_type)
        #print( molecular_formulas.ion_charge)

        
            
if __name__ == "__main__":
    
    time0 = time.time()
    create_lookup_table()
    print(time.time()-time0 )