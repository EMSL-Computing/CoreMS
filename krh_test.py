from corems.molecular_formula.factory.MolecularFormulaFactory import MolecularFormulaBase

mol_form = MolecularFormulaBase({'C': 6, 'H': 10, 'O': 6, 'IonType': 'ADDUCT', 'Cl': 1}, ion_charge = -1)
print(mol_form._d_molecular_formula)
print(mol_form.string)
