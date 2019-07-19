import sys
sys.path.append(".")

from enviroms.emsl.yec.molecular_id.factory.MolecularFormulaFactory import MolecularFormula 

if __name__ == "__main__":
    
    formula_dict = {'C':20, 'H':12, 'O':100, 'IonType': 'Radical'}
    ion_charge = 1 
    formula_obj = MolecularFormula(formula_dict, ion_charge)
    print(formula_obj.ion_type)
    print(formula_obj.mz_theor)
    print(formula_obj._cal_isotopologues(formula_dict))