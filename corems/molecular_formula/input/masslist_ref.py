__author__ = "Yuri E. Corilo"
__date__ = "Oct 24, 2019"

from threading import Thread
from pathlib import Path
import sys, re, json
from typing import Dict, List

import pandas as pd

sys.path.append('.')

from corems.molecular_formula.factory.MolecularFormulaFactory import LCMSLibRefMolecularFormula, MolecularFormula 
from corems.encapsulation.constant import Labels
from corems.encapsulation.constant import Atoms
from corems.molecular_id.factory.molecularSQL import CarbonHydrogen, MolecularFormulaLink, HeteroAtoms

class MolecularFormulaLinkProxy:
        
        def __init__(self, molecular_formula, mz):

            self.C = molecular_formula.get('C')
            self.H = molecular_formula.get('H')
            self.H_C = molecular_formula.get('H')/molecular_formula.get('C')
            self.class_label = json.dumps(molecular_formula.class_dict)
            self.mz_calc =  float(mz)                       
            self.dbe = molecular_formula.dbe
            self.formula_dict = molecular_formula.to_dict()
        
        def to_dict(self):
            return self.formula_dict


class ImportMassListRef():#Thread

    def __init__(self, ref_file_location) :
            
            #Thread.__init__(self)
            
            self.ref_file_location = Path(ref_file_location)
            
            if not self.ref_file_location.exists():
                tb = sys.exc_info()[2]
                raise FileNotFoundError(ref_file_location).with_traceback(tb)
    
    def molecular_formula_ref( self, mz, molecular_formula):
        
        return MolecularFormulaLinkProxy(molecular_formula, mz)
    
    def from_lcms_lib_file(self, ion_charge: float, ion_types: List[str]) -> Dict [str, Dict[float, List[LCMSLibRefMolecularFormula] ] ]:
        '''
         return Dict[standard_name, Dict[m/z, List[MolecularFormula]]]:
            m/z: float:
                target m/z 
            standard_name: str
                name of the molecular standard mix 
            MolecularFormula: class
                corems molecular formula class    
        '''
        
        data = {}
        
        with open(self.ref_file_location) as ref_f:
            
            df = pd.read_csv(ref_f, header=0,  encoding= 'unicode_escape')
            
            print(df)

            for index, row in df.iterrows():
                
                formula_s = row["Neutral Formula"]
                formula_dict = self.mformula_s_to_dict(formula_s, Labels.neutral)
                name = row["Compound Name"]
                kegg_id = row["KEGG ID"]
                standard_name = row["NEW MIX"]
                cas = row["KEGG ID"]
                #print(row["Neutral Formula"], formula_dict)
                molf_formula = LCMSLibRefMolecularFormula(formula_dict, ion_charge, Labels.neutral, 
                                                name=name, kegg_id=kegg_id, cas=cas)
                #if round(molf_formula.mz_calc, 4) != round(row['Mass Adduct -H'],4):
                #    print(formula_s)
                #    print(round(molf_formula.mz_calc, 4) , round(row['Mass Adduct -H'],4))
        
                if standard_name in data.keys():

                    #TODO change it to target ion types and add ion type in the data structure   
                    mz_calc = molf_formula.protonated_mz
                    
                    if mz_calc in data.get(standard_name).keys():
                       
                       data.get(standard_name).get(mz_calc).append(molf_formula)
                    
                    else:   
                        data[standard_name][mz_calc] = [molf_formula]
                else:
                    
                    data[standard_name] = {molf_formula.mz_calc: [molf_formula]}
                #print(formula_s, formula_dict)
                #if molf_formula.ion_type != 'de-protonated':
                #    print( 'ha', molf_formula.ion_type )
                #print(formula_dict)
                #print(row['c1'], row['c2'])
        
        return data

    def from_bruker_ref_file(self) -> List[MolecularFormula]:

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
                    mz = float(list_ref[1])
                    formula_dict = self.mformula_s_to_dict(ion_mol_formula)
                    
                    list_mf_obj.append(MolecularFormula(formula_dict, ion_charge, external_mz=mz))
        
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

    def mformula_s_to_dict(self, s_mformulatring, iontype='unknown'):
        
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
            single_digit_atoms_one = re.findall(r'[A-Z]{1}(?![0-9])(?![a-z])', s_mformulatring)
            #print(single_digit_atoms_one)
            #find the case Na
            due_digit_atoms_one = re.findall(r'[A-Z]{1}[a-z]{1}(?![0-9])', s_mformulatring)
            
            all_atoms = all_atoms + all_atoms2 + due_digit_atoms_one +single_digit_atoms_one
            
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
            
            dict_res[Labels.ion_type]  = iontype

            return dict_res
        
        else: 
            
            tb = sys.exc_info()[2]
            raise Exception('Empty molecular formula').with_traceback(tb)

    
    
        
                