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
    
    formula_dict = {'C':10, 'H':0, 'O':10,'P':1, 'IonType': 'Radical'}
    #formula_dict = OrderedDict(formula_dict)
    ion_charge = 1 
    formula_obj = MolecularFormula(formula_dict, ion_charge)
    
    #print(formula_obj.ion_type)
    #p#rint(formula_obj.mz_theor)
    
    #print(formula_obj.isotopologues[1]._d_molecular_formula == formula_obj.isotopologues[4]._d_molecular_formula)
    #print(formula_obj.isotopologues[1]._d_molecular_formula, formula_obj.isotopologues[4]._d_molecular_formula)
    #for i in mass.isotopologues(formula='C100'):
    #    print(i)
    #'''
    for isotopologue_obj in formula_obj.isotopologues:
        pass
        print(isotopologue_obj.mz_theor, isotopologue_obj.prop_ratio, isotopologue_obj.to_string)
    
   
    