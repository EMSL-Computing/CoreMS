__author__ = "Yuri E. Corilo"
__date__ = "Jul 29, 2019"

import multiprocessing
import os
from os.path import join

from enviroms.emsl.yec.encapsulation.constant.Constants import Labels
from enviroms.emsl.yec.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaLookupTableSettings, MoleculaSearchSettings
from enviroms.emsl.yec.mass_spectrum.input.TextMassList import Read_MassList
from enviroms.emsl.yec.molecular_id.calc.MolecularLookupTable import  MolecularCombinations


class SearchMolecularFormulas:
     
    '''
    runworker()
    
    '''

    def runworker_mspeak(self, mspeak, mass_spectrum_obj):

        MoleculaLookupTableSettings.min_mz = mspeak.mz_exp-1
    
        MoleculaLookupTableSettings.max_mz = mspeak.mz_exp+1
        
        min_abundance = mass_spectrum_obj.min_abundance

        dict_molecular_lookup_table = MolecularCombinations().runworker()

        nominal_mz  = mspeak.nominal_mz_exp

        classes = list(dict_molecular_lookup_table.keys())
        
        '''
        waiting for python 3.8 release to set mass_spectrum_obj and dict_molecular_lookup_table on share memory
        pool = multiprocessing.Pool(number_of_process)
        args = [ (dict_molecular_lookup_table.get(classe).get(ion_type).get(nominal_mz), min_abundance, mass_spectrum_obj, mspeak_mz_exp, mspeak_abundance)  for classe in classes ]
        pool.map(SearchMolecularFormulaWorker(), args)

        pool.close()
        pool.join()
        '''
        for classe in classes:
            
            possible_formulas = list()    
            #we might need to increase the search space to -+1 m_z 
            if MoleculaSearchSettings.isRadical:
            
                ion_type = Labels.radical_ion
                
                formulas = dict_molecular_lookup_table.get(classe).get(ion_type).get(nominal_mz)
                
                if formulas:
                    
                    possible_formulas.extend(formulas)

            if MoleculaSearchSettings.isProtonated:
            
                ion_type = Labels.protonated_de_ion

                formulas = dict_molecular_lookup_table.get(classe).get(ion_type).get(nominal_mz)
                
                if formulas:
                    
                    possible_formulas.extend(formulas)
            
            
            if possible_formulas:
                
                SearchMolecularFormulaWorker().find_formulas(possible_formulas, min_abundance, mass_spectrum_obj, mspeak)
        


    def runworker_mass_spectrum(self, mass_spectrum_obj):

        
        #number_of_process = multiprocessing.cpu_count()

        '''loading this on a shared memory would be better than having to serialize it for every process
            waiting for python 3.8 release'''

        MoleculaLookupTableSettings.min_mz = mass_spectrum_obj.min_mz_exp
    
        MoleculaLookupTableSettings.max_mz = mass_spectrum_obj.max_mz_exp
        
        min_abundance = mass_spectrum_obj.min_abundance

        dict_molecular_lookup_table = MolecularCombinations().runworker()

        classes = list(dict_molecular_lookup_table.keys())

        print(len(mass_spectrum_obj))
        
        for mspeak in  mass_spectrum_obj:

            #print(mspeak) 
            nominal_mz  = mspeak.nominal_mz_exp
            
            '''
            waiting for python 3.8 release to set mass_spectrum_obj and dict_molecular_lookup_table on share memory
            pool = multiprocessing.Pool(number_of_process)
            args = [ (dict_molecular_lookup_table.get(classe).get(ion_type).get(nominal_mz), min_abundance, mass_spectrum_obj, mspeak_mz_exp, mspeak_abundance)  for classe in classes ]
            pool.map(SearchMolecularFormulaWorker(), args)

            pool.close()
            pool.join()
            '''
            for classe in classes:
                
                possible_formulas = list()    
                #we might need to increase the search space to -+1 m_z 
                if MoleculaSearchSettings.isRadical:
                
                    ion_type = Labels.radical_ion
                    
                    formulas = dict_molecular_lookup_table.get(classe).get(ion_type).get(nominal_mz)
                    
                    if formulas:
                        
                        possible_formulas.extend(formulas)

                if MoleculaSearchSettings.isProtonated:
                
                    ion_type = Labels.protonated_de_ion

                    formulas = dict_molecular_lookup_table.get(classe).get(ion_type).get(nominal_mz)
                    
                    if formulas:
                        
                        possible_formulas.extend(formulas)

                if possible_formulas:
                    
                    SearchMolecularFormulaWorker().find_formulas(possible_formulas, min_abundance, mass_spectrum_obj, mspeak)
            

class SearchMolecularFormulaWorker:

    # needs this wraper to pass the class to multiprocessing
    def __call__(self, args):

        return self.find_formulas(*args)  # ,args[1]

    def find_formulas(self, possible_formulas, min_abundance, mass_spectrum_obj, mspeak ):
        
        min_mz_error = MoleculaSearchSettings.min_mz_error
        max_mz_error = MoleculaSearchSettings.max_mz_error
        min_abun_error = MoleculaSearchSettings.min_abun_error
        max_abun_error = MoleculaSearchSettings.max_abun_error

        #f = open("abundance_error.txt", "a+")    
        mspeak_mz_exp, mspeak_abundance = mspeak.mz_exp, mspeak.abundance

        for possible_formula in possible_formulas:
            
            if possible_formula:
                
                error = possible_formula._calc_assigment_mass_error(mspeak_mz_exp)
                
                if  min_mz_error <= error <= max_mz_error:
                    
                    #add molecular formula match to mspeak
                    mspeak.add_molecular_formula(possible_formula)

                    #calculates and look for isotopologues
                    isotopologues = possible_formula.isotopologues(min_abundance, mspeak_abundance)
                    
                    for isotopologue_formula in isotopologues:
                        
                        #move this outside to impove preformace
                        #we need to increase the search space to -+1 m_z 
                        first_index, last_index = mass_spectrum_obj.get_nominal_mz_frist_last_indexes(isotopologue_formula.mz_nominal_theo)
                        
                        for mspeak_iso in mass_spectrum_obj[first_index:last_index]:
                            
                            error = isotopologue_formula._calc_assigment_mass_error(mspeak_iso.mz_exp)    
                            
                            #need to define error distribution for abundance measurements
                            if  min_mz_error <= error <= max_mz_error:
                                    
                                    abundance_error = isotopologue_formula._calc_abundance_error(mspeak_abundance,mspeak_iso.abundance )            
                                    
                                    # margin of error was set empirically/ needs statistical calculation
                                    #  of margin of error for the measurement of the abundances
                                    if min_abun_error <= abundance_error <= max_abun_error:
                                        
                                        mspeak_iso.add_molecular_formula(isotopologue_formula)
                                        
                                        #print(error, abundance_error)  
                                        
                                        #output_file.write("%.4f \t %.4f \t %.2f \n " % (mspeak_iso.mz_exp,  mspeak_iso.abundance,  abundance_error))     
        #f.close()                                