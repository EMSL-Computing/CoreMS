__author__ = "Yuri E. Corilo"
__date__ = "Jul 22, 2019"

import sys, pytest
sys.path.append(".")
from enviroms.encapsulation.constant.Constants import Atoms
from enviroms.molecular_id.factory.MolecularFormulaFactory import MolecularFormula    


from copy import deepcopy

def test_molecular_formula():
    
    formula_dict = {'C':10, 'H':0, 'O':10,'Cl':2, 'IonType': 'Radical'}
    
    ion_charge = 1 
    formula_obj = MolecularFormula(formula_dict, ion_charge)
    print("ion_type", formula_obj.ion_type)
    print("mz_theor", formula_obj.mz_theor)
    
    min_abudance, current_abundance = 1,1 
    
    for isotopologue_obj in formula_obj.isotopologues(0.01, current_abundance):
        
        print("formula:", isotopologue_obj.to_string, 
              "mz_theor:", isotopologue_obj.mz_theor,
              "prprop_ratio:", isotopologue_obj.prop_ratio)
    
   

    