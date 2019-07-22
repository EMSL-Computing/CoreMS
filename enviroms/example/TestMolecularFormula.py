__author__ = "Yuri E. Corilo"
__date__ = "Jul 22, 2019"

import sys
sys.path.append(".")
from collections import OrderedDict
from enviroms.emsl.yec.encapsulation.constant.Constants import Atoms
from enviroms.emsl.yec.molecular_id.factory.MolecularFormulaFactory import MolecularFormula
from scipy.stats import binom
from numpy import where, exp
from IsoSpecPy import IsoSpecPy

from copy import deepcopy

if __name__ == "__main__":
    
    formula_dict = {'C':10, 'H':0, 'O':10,'Cl':1, 'IonType': 'Radical'}
    
    ion_charge = 1 
    formula_obj = MolecularFormula(formula_dict, ion_charge)
    print("ion_type", formula_obj.ion_type)
    print("mz_theor", formula_obj.mz_theor)
    
    for isotopologue_obj in formula_obj.isotopologues:
        
        print("formula:", isotopologue_obj.to_string, 
              "mz_theor:", isotopologue_obj.mz_theor,
              "prprop_ratio:", isotopologue_obj.prop_ratio)
    
   

    