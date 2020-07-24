__author__ = "Yuri E. Corilo"
__date__ = "Oct 24, 2019"

from threading import Thread
from pathlib import Path
import sys, re, json

sys.path.append('.')

from corems.molecular_formula.factory.MolecularFormulaFactory import MolecularFormula 
from corems.encapsulation.constant import Labels
from corems.encapsulation.constant import Atoms
from corems.molecular_id.factory.molecularSQL import CarbonHydrogen, MolecularFormulaLink, HeteroAtoms

class MolecularFormulaLinkProxy():
        
        def __init__(self, molecular_formula, mz):

            self.C = molecular_formula.get('C')
            self.H = molecular_formula.get('H')
            self.H_C = molecular_formula.get('H')/molecular_formula.get('C')
            self.classe = json.dumps(molecular_formula.class_dict)
            self.mass =  float(mz)                       
            self.dbe = molecular_formula.dbe
            self.formula_dict = molecular_formula.to_dict()

class ImportMassListRef():#Thread

    def __init__(self, ref_file_location) :
            
            #Thread.__init__(self)
            
            self.ref_file_location = Path(ref_file_location)
            
            if not self.ref_file_location.exists():
                tb = sys.exc_info()[2]
                raise FileNotFoundError(ref_file_location).with_traceback(tb)
    
    def molecular_formula_ref( self, mz, molecular_formula):
        
        return MolecularFormulaLinkProxy(molecular_formula, mz)
    
    def from_bruker_ref_file(self):

        import csv
        
        list_mf_obj = []

        with open(self.ref_file_location) as ref_f:

            labels = ref_f.readline().strip('\n').split(';')
            
            for line in ref_f.readlines():
                
                if line != '\n':
        
                    list_ref = (line.strip('\n').split(' '))
                    
                    if list_ref[2][-1] == '+': 
                        
                        ion_charge =  int(list_ref[2][:-1])
                    
                    else:
                        
                        ion_charge =  -1* int(list_ref[2][:-1])

                    
                    ion_mol_formula = list_ref[0]
                    mz = list_ref[1]
                    formula_dict = self.mformula_s_to_dict(ion_mol_formula)
                    
                    list_mf_obj.append(self.molecular_formula_ref(mz, MolecularFormula(formula_dict, ion_charge)))
        
        return  list_mf_obj           

    def from_corems_ref_file(self, delimiter="\t"): #pragma: no cover
        '''not being used'''
        import csv

        list_mf_obj = []

        with open('res/RefMassLists/Crude-Pos-ESI.ref') as ref_f:

            labels = ref_f.readline().strip('\n').split(delimiter)
            
            for line in ref_f.readlines():
                
                if line != '\n':
        
                    list_ref = (line.strip('\n').split(delimiter))
                    
                    formula_string = list_ref[0]
                    ion_charge = int(list_ref[1])
                    ion_type = list_ref[2]

                    molform = MolecularFormula(formula_string, ion_charge, ion_type=ion_type)
                    
                    list_mf_obj.append(self.molecular_formula_ref(molform))

        return  list_mf_obj           


    def split(self, delimiters, string, maxsplit=0): #pragma: no cover
    
        ''' does not work when formula has atoms with same caracaters:
            i.e - C10H21NNa
        '''
        regexPattern = '|'.join(map(re.escape, delimiters)) #pragma: no cover
        isotopes = re.findall(regexPattern, string) #pragma: no cover
        counts = re.split(regexPattern, string, maxsplit)  #pragma: no cover
        return isotopes, counts

    def mformula_s_to_dict(self, s_mformulatring):
        
        ''' 
            Converts a molecular formula string to a dict
            
            args:
            
            s_mformulatring: str
                'C10H21NNa'
            
            -   Does not work if the atomic mass number is passed i.e. 37Cl, 81Br
                The convention follow the light isotope labeling 35Cl is Cl, 12C is C, etc
            
           -    If you need to use heavy isotopes please use another reference file format that 
                separate the formula string by a blank space and parse it using the function corems_ref_file""   
        
        '''

        if s_mformulatring:
            
            #find the case C122
            all_atoms = re.findall(r'[A-Z]{1}[0-9]{1,10000}', s_mformulatring)
            
            #find the case Br2
            all_atoms2 = re.findall(r'[A-Z]{1}[a-z]{1}[0-9]{1,10000}', s_mformulatring)

            #find the case N
            single_digit_atoms_one = re.findall(r'[A-Z]{1}(?![0-9])(?![A-Z])', s_mformulatring)

            #find the case Na
            due_digit_atoms_one = re.findall(r'[A-Z]{1}[a-z]{1}(?![0-9])', s_mformulatring)

            all_atoms = all_atoms + all_atoms2 + single_digit_atoms_one + due_digit_atoms_one
            
            dict_res = {}
            
            for each_atom_count in all_atoms:
                
                
                count = re.findall(r'[0-9]{1,10000}', each_atom_count)
                atom = ''.join(re.findall(r'[A-z]', each_atom_count))
                
                if atom in Atoms.atoms_order:
                    
                    if count:
                        dict_res[atom] = int(count[0])
                    else:
                        dict_res[atom] = 1
                
                else:
                    
                    tb = sys.exc_info()[2]
                    raise TypeError("Atom %s does not exist in Atoms.atoms_order list" % atom).with_traceback(tb)
            
            dict_res[Labels.ion_type]  =  'unknown'

            return dict_res
        
        else: 
            
            tb = sys.exc_info()[2]
            raise Exception('Empty molecular formula').with_traceback(tb)

    
    
        
                