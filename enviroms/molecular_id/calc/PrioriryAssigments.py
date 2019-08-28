
from os import path
from threading import Thread
import os, sys
sys.path.append('.')
from enviroms.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaLookupTableSettings, MoleculaSearchSettings
from enviroms.molecular_id.calc.FindOxigenPeaks import FindOxygenPeaks

class OxigenPriorityAssignment(Thread):

    def __init__(self, mass_spectrum_obj, lookupTableSettings):
        
        Thread.__init__(self)
        self.mass_spectrum_obj = mass_spectrum_obj
        self.lookupTableSettings = lookupTableSettings

    def run(self):
        
        lowest_error = lambda msp:  msp.molecular_formula_lowest_error

        find_formula_thread = FindOxygenPeaks(self.mass_spectrum_obj, self.lookupTableSettings)
        find_formula_thread.run()
        find_formula_thread.set_mass_spec_indexes_by_found_peaks()
        
        mspeaks =[ mspeak.molecular_formula_lowest_error.class_label for mspeak  in sorted(self.mass_spectrum_obj, key=lambda msp: msp.mz_exp)]
        
        print(mspeaks)

        print (self.mass_spectrum_obj[0].molecular_formula_lowest_error['O'])
        
        min_o_mspeak = min(self.mass_spectrum_obj,key=lambda msp: msp.molecular_formula_lowest_error['O']) 
        max_o_mspeak = max(self.mass_spectrum_obj,key=lambda msp: msp.molecular_formula_lowest_error['O']) 

        dbe = [lowest_error(mspeak).dbe for mspeak in self.mass_spectrum_obj]
        oxigen = [lowest_error(mspeak)['O'] for mspeak in self.mass_spectrum_obj]
        carbons = [lowest_error(mspeak)['C'] for mspeak in self.mass_spectrum_obj]
        
        print(dbe)
        print(oxigen)
        print(carbons)
        
        self.lookupTableSettings.usedAtoms['O'] = (min_o-1, max_o+1)

        print (min_o.molecular_formula_lowest_error['O'], max_o.molecular_formula_lowest_error['O'])

    def run_worker_mass_spectrum(self):
        
        last_dif = 0
    
        last_error = 0
        
        closest_error = 0

        error_average = MoleculaSearchSettings.mz_error_average
        
        nbValues = 0

        settings.min_mz = mass_spectrum_obj.min_mz_exp
    
        settings.max_mz = mass_spectrum_obj.max_mz_exp
        
        min_abundance = mass_spectrum_obj.min_abundance

        dict_molecular_lookup_table = MolecularCombinations().runworker(settings)

        classes = list(dict_molecular_lookup_table.keys())

        print(len(mass_spectrum_obj))
        
        for ms_peak in sorted(mass_spectrum_obj, key=lambda m :m.mz_exp):

            #print(ms_peak) 
            nominal_mz  = ms_peak.nominal_mz_exp
            
            '''
            waiting for python 3.8 release to set mass_spectrum_obj and dict_molecular_lookup_table on share memory
            pool = multiprocessing.Pool(number_of_process)
            args = [ (dict_molecular_lookup_table.get(classe).get(ion_type).get(nominal_mz), min_abundance, mass_spectrum_obj, ms_peak_mz_exp, ms_peak_abundance)  for classe in classes ]
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
                    
                    SearchMolecularFormulaWorker().find_formulas(possible_formulas, min_abundance, mass_spectrum_obj, ms_peak, last_error, last_dif, closest_error, error_average, nbValues)
            

if __name__ == "__main__":
    
    
    from enviroms.transient.input.BrukerSolarix import ReadBrukerSolarix
    from enviroms.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaSearchSettings, MoleculaLookupTableSettings

    file_location = os.path.join(os.getcwd(), "tests/tests_data/") + os.path.normcase("ESI_NEG_SRFA.d/")
    
    bruker_reader = ReadBrukerSolarix(file_location)

    bruker_transient = bruker_reader.get_transient()

    mass_spectrum_obj = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=True)

    lookupTableSettings = MoleculaLookupTableSettings()
    
    lookupTableSettings.usedAtoms['O'] = (1, 30)

    assignOx = OxigenPriorityAssignment(mass_spectrum_obj, lookupTableSettings)

    assignOx.start()
