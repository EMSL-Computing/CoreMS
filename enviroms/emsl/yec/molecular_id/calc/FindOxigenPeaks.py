__author__ = "Yuri E. Corilo"
__date__ = "Jul 31, 2019"

from threading import Thread
from enviroms.emsl.yec.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaLookupTableSettings, MoleculaSearchSettings
from enviroms.emsl.yec.molecular_id.factory.MolecularFormulaFactory import MolecularFormula
from enviroms.emsl.yec.molecular_id.calc.MolecularLookupTable import  MolecularCombinations
from enviroms.emsl.yec.molecular_id.calc.MolecularFormulaSearch import SearchMolecularFormulas
from copy import deepcopy

class FindOxygenPeaks(Thread):
    '''class to walk 14Da units over oxygen space for negative ion mass spectrum of natural organic matter
        Returns a list of MSPeak class cotaining the possibles Molecular Formula class objects.  
    '''
    def __init__(self, mass_spectrum_obj):
        
        Thread.__init__(self)
        
        self.mass_spectrum_obj = mass_spectrum_obj
        
    def run(self):
        
        molecular_formula_obj_reference = self.find_most_abundant_peak(self.mass_spectrum_obj)
        
        possible_mol_formulas_objs = self.build_database(molecular_formula_obj_reference)
        
        self.list_found_mspeaks = self.look_for_it_now(molecular_formula_obj_reference, possible_mol_formulas_objs, self.mass_spectrum_obj)
        
        
    def find_most_abundant_peak(self, mass_spectrum_obj):
        
        #this function is intended for test only. 
        # Have to sort by Kendrick to be able to select the most abundant series 
        #then select the most abundant peak inside the series
        #or have the user select the reference mspeak somehow

        mspeak_most_abundant = mass_spectrum_obj.most_abundant_mspeak

        SearchMolecularFormulas().runworker_mspeak(mspeak_most_abundant, mass_spectrum_obj)

        if mspeak_most_abundant:

            return mspeak_most_abundant.molecular_formula_lowest_error 
        else:
            raise Exception("Could not find a possible molecular formula match for the most abudant peak of m/z %.5f"%mspeak_most_abundant.mz_exp )
        #return the first option
        #return mspeak_most_abundant[0]
        
    def get_list_found_peaks(self):
        
        return self.list_found_mspeaks
    
    def look_for_it_now(self, mol_formula_obj_reference, possible_mol_formulas_objs, mass_spectrum_obj):
       
        min_mz_error = mol_formula_obj_reference.mz_error - 2
        max_mz_error = mol_formula_obj_reference.mz_error + 2
        
        list_found_peaks = list()
        
        for possible_formula in possible_mol_formulas_objs:
            
            first_index, last_index = mass_spectrum_obj.get_nominal_mz_frist_last_indexes(possible_formula.mz_nominal_theo)
                        
            for mspeak in mass_spectrum_obj[first_index:last_index]:
                
                error = possible_formula._calc_assigment_mass_error(mspeak.mz_exp)    
                
                #need to define error distribution for abundance measurements
                if  min_mz_error <= error <= max_mz_error:
                    
                    mspeak.add_molecular_formula(possible_formula)

                    list_found_peaks.append(mspeak)

        return list_found_peaks     
            
            
    def build_database(self, initial_formula_obj):
        
        ion_charge = initial_formula_obj.ion_charge
        
        initial_mass = initial_formula_obj.mz_theor
            
        initial_formula = initial_formula_obj.to_dict
        
        min_mz = self.mass_spectrum_obj.min_mz_exp
        
        max_mz = self.mass_spectrum_obj.max_mz_exp
        
        list_formulas = []
            
        i = 0
        mass = initial_mass
        
        while mass > min_mz:
            #print "shit 1", mass, min_mz
            new_formula = deepcopy(initial_formula)
            new_formula['C'] = initial_formula.get('C') - i
            new_formula['O'] = initial_formula.get('O') - i
            molecular_formula = MolecularFormula(new_formula, ion_charge)
            list_formulas.append(molecular_formula)
            
            new_formula = deepcopy(initial_formula)
            new_formula['C'] = initial_formula.get('C') - i+1
            new_formula['H'] = initial_formula.get("H") + (2)
            new_formula['O'] = initial_formula.get("O") - i
            molecular_formula = MolecularFormula(new_formula, ion_charge)
            list_formulas.append(molecular_formula)
            mass = molecular_formula.mz_theor
            i = i + 1
        
        i = 1
        mass = initial_mass
        while mass < max_mz:
            
            new_formula = deepcopy(initial_formula)
            new_formula['C'] = initial_formula.get('C') + i
            new_formula['O'] = initial_formula.get('O') + i
            molecular_formula = MolecularFormula(new_formula, ion_charge)
            list_formulas.append(molecular_formula)
            
            
            new_formula = deepcopy(initial_formula)
            new_formula['C'] = initial_formula.get('C') + i+1
            new_formula['H'] = initial_formula.get("H") + (2)
            new_formula['O'] = initial_formula.get("O") + i
            molecular_formula = MolecularFormula(new_formula, ion_charge)
            list_formulas.append(molecular_formula)
            mass = molecular_formula.mz_theor
            i = i + 1
        
        list_formulas = sorted(list_formulas, key=lambda mf: mf.mz_theor)
        
        for mformula in list_formulas:
            print('created',mformula.to_string, mformula.mz_theor)
        
        return list_formulas
        