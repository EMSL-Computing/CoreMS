__author__ = "Yuri E. Corilo"
__date__ = "Jul 22, 2019"


import time
from multiprocessing import freeze_support
import sys
sys.path.append(".")
from enviroms.emsl.yec.encapsulation.constant.Constants import Atoms
from enviroms.emsl.yec.molecular_id.calc.MolecularLookupTable import  MolecularCombinations
from enviroms.emsl.yec.encapsulation.settings.molecular_id.MolecularIDSettings import MolecularSpaceTableSetting

def create_lookup_table():
    
    
    #margin_error needs to be optimized by the data rp and sn
    margin_error = 0.3

    #min_mz,max_mz  needs to be optimized by the data
    min_mz = 100
    max_mz = 1300

    # C, H, N, O, S and P atoms are ALWAYS needed in the dictionary
    #the defaults values are defined at the encapsulation MolecularSpaceTableSetting    
    MolecularSpaceTableSetting.usedAtoms['C'] = (1,90)
    
    #some atoms has more than one valence state and the most commun will be used
    # adduct atoms needs valence 0
    MolecularSpaceTableSetting.usedAtoms['Cl'] = (0,0)
    possible_valences = Atoms.atoms_valence.get('Cl')
    valence_one = possible_valences[0]
    # if you want to specify it in needs to be changed here
    MolecularSpaceTableSetting.used_atom_valences['Cl'] =  valence_one
    #otherwise it will use the lowest valence, PS needs insure propagation to isotopologues
    
    dict_molecular_formulas = MolecularCombinations().runworker( )

    for molecular_formulas in dict_molecular_formulas.get("O10").get(401):
        print( molecular_formulas.class_label)
        print( molecular_formulas.to_string)
        print( molecular_formulas.mz_theor)
        for isotopologue in molecular_formulas.isotopologues:
            print("isotopologue", isotopologue.to_string)
        #print( molecular_formulas.atoms)
        #print( molecular_formulas.ion_type)
        #print( molecular_formulas.ion_charge)
            
if __name__ == "__main__":
    
    time0 = time.time()
    create_lookup_table()
    print(time.time()-time0 )
