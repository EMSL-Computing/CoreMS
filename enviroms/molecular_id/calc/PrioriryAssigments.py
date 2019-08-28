
import os,  sys
sys.path.append('.')
from threading import Thread
from numpy.core.umath import invert
from enviroms.molecular_id.calc.FindOxigenPeaks import FindOxygenPeaks
from enviroms.molecular_id.calc.MolecularFormulaSearch import SearchMolecularFormulas

class OxigenPriorityAssignment(Thread):

    def __init__(self, mass_spectrum_obj, lookupTableSettings):
        
        Thread.__init__(self)
        self.mass_spectrum_obj = mass_spectrum_obj
        self.lookupTableSettings = lookupTableSettings

    def run(self):
        
        lowest_error = lambda msp:  msp.molecular_formula_lowest_error

        find_formula_thread = FindOxygenPeaks(self.mass_spectrum_obj, self.lookupTableSettings)
        find_formula_thread.run()
        #mass spec obj indexes are set to interate over only the peaks with a molecular formula candidate
        find_formula_thread.set_mass_spec_indexes_by_found_peaks()
        
        dict_class_and_ms_peak_index = {}
        
        for mspeak in sorted(self.mass_spectrum_obj, key=lambda msp: msp.abundance, reverse=True):
            
            classe = mspeak.molecular_formula_lowest_error.class_label
            
            if classe not in dict_class_and_ms_peak_index.keys():
                dict_class_and_ms_peak_index[classe] = mspeak
        
        # reset back the massspec obj indexes since we already collected the mspeak indexes
        self.mass_spectrum_obj.reset_indexes()
        
        for classe, mspeak in dict_class_and_ms_peak_index.items():
            
            oxigen_number = mspeak[0]['O']
            carbon_number = mspeak[0]['C']
            dbe = mspeak[0].dbe

            self.lookupTableSettings.use_pah_line_rule = False
            self.lookupTableSettings.min_dbe = dbe - 7
            self.lookupTableSettings.max_dbe = dbe + 7
            self.lookupTableSettings.usedAtoms['O'] =  (oxigen_number, oxigen_number)
            self.lookupTableSettings.usedAtoms['C'] =  (carbon_number-4, carbon_number+8) 
            
            SearchMolecularFormulas(first_hit=True).run_worker_mass_spectrum(self.mass_spectrum_obj, self.lookupTableSettings)
            
            print(classe) 
        
        print('before', len(self.mass_spectrum_obj))
        print('before', len(self.mass_spectrum_obj._mspeaks))
        
        self.mass_spectrum_obj.reset_indexes()
        print('after', len(self.mass_spectrum_obj))
        print('after', len(self.mass_spectrum_obj._mspeaks))

if __name__ == "__main__":
    
    from enviroms.transient.input.BrukerSolarix import ReadBrukerSolarix
    from enviroms.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaSearchSettings, MoleculaLookupTableSettings
    from enviroms.mass_spectrum.calc.CalibrationCalc import MZDomain_Calibration, FreqDomain_Calibration
    from matplotlib import pyplot, colors as mcolors

    def calibrate():
        
        MoleculaSearchSettings.error_method = 'average'
        MoleculaSearchSettings.min_mz_error = -5
        MoleculaSearchSettings.max_mz_error = 1
        MoleculaSearchSettings.mz_error_range = 1

        find_formula_thread = FindOxygenPeaks(mass_spectrum_obj, lookupTableSettings)
        find_formula_thread.run()
        mspeaks_results = find_formula_thread.get_list_found_peaks()
        
        calibrate = FreqDomain_Calibration(mass_spectrum_obj, mspeaks_results)
        calibrate.ledford_calibration()
        mass_spectrum_obj.clear_molecular_formulas()

        MoleculaSearchSettings.error_method = 'symmetrical'
        MoleculaSearchSettings.min_mz_error = -1
        MoleculaSearchSettings.max_mz_error = 1
        MoleculaSearchSettings.mz_error_range = 2
        MoleculaSearchSettings.mz_error_average = 0
        MoleculaSearchSettings.min_abun_error = -30 # percentage
        MoleculaSearchSettings.max_abun_error = 70 # percentage
        MoleculaSearchSettings.isProtonated = True
        MoleculaSearchSettings.isRadical= True
    

    def plot():
        colors = list(mcolors.XKCD_COLORS.keys())
        oxigens = range(6,21)
        for o in oxigens:
            o_c = list()
            for mspeak in mass_spectrum_obj:
                if mspeak:
                    #molecular_formula = mspeak.molecular_formula_lowest_error
                    for molecular_formula in mspeak:
                        if molecular_formula['O'] == o:
                            if  not molecular_formula.is_isotopologue:
                                
                                pyplot.plot(molecular_formula['C'], molecular_formula.dbe, "o",   color=colors[molecular_formula['O']])
                                pyplot.plot(molecular_formula['C'], molecular_formula.dbe, "o",   color=colors[molecular_formula['O']])
                                pyplot.annotate(molecular_formula['O'], (molecular_formula['C']+0.5, molecular_formula.dbe+0.5))

            pyplot.show()                
    
    file_location = os.path.join(os.getcwd(), "tests/tests_data/") + os.path.normcase("ESI_NEG_SRFA.d/")
    
    bruker_reader = ReadBrukerSolarix(file_location)

    bruker_transient = bruker_reader.get_transient()

    mass_spectrum_obj = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=True)

    lookupTableSettings = MoleculaLookupTableSettings()
    
    lookupTableSettings.usedAtoms['O'] = (1, 22)
    
    calibrate()

    assignOx = OxigenPriorityAssignment(mass_spectrum_obj, lookupTableSettings)

    assignOx.start()

    assignOx.join()

    plot()



    