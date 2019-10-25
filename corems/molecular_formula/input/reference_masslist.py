__author__ = "Yuri E. Corilo"
__date__ = "Oct 24, 2019"

from threading import Thread
from pathlib import Path
import sys

sys.path.append('.')
from corems.molecular_formula.factory.MolecularFormulaFactory import MolecularFormula 
from corems.encapsulation.constant import Labels

class ReadMassListReference(Thread):

    def __init__(self, mass_spectrum_obj, ref_file_location) :
            
            Thread.__init__(self)
            
            self.ref_file_location = Path(ref_file_location)
            
            if not self.ref_file_location.exists():
                tb = sys.exc_info()[2]
                raise FileNotFoundError(ref_file_location).with_traceback(tb)

            self.mass_spectrum_obj = mass_spectrum_obj
            
    def bruker_ref_file(self):

        import csv
        ref_f = open('samples.csv')
        rdr = csv.DictReader(filter(lambda row: row[0]!='#', ref_f))
        for row in rdr:
            print(row)
        ref_f.close()
        
if __name__ == "__main__":
    
    import csv
    
    with open('res/RefMassLists/PEG.ref') as ref_f:

        labels = ref_f.readline().strip('\n').split(';')
        
        for line in ref_f.readlines():
            
            if line != '\n':
    
                list_ref = (line.strip('\n').split(' '))
                
                if list_ref[2][-1] == '+': 
                    
                    ion_charge =  int(list_ref[2][:-1])
                
                else:
                    
                    ion_charge =  -1* int(list_ref[2][:-1])
                
                ion_mol_formula =   list_ref[3]

                print(ion_charge, ion_mol_formula)
                
                formula_dict = {'C':10, 'H':0, 'O':10,'Cl':2, Labels.ion_type: 'Unkown'}
                
                formula_obj = MolecularFormula(formula_dict, ion_charge)
        
        