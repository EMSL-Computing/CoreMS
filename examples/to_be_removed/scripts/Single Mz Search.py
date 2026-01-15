import pprint, sys
from pathlib import Path


sys.path.append("./")
from corems.encapsulation.factory.parameters import MSParameters
from corems.encapsulation.input import parameter_from_json
from corems.mass_spectrum.input.numpyArray import ms_from_array_centroid
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas


def run_molecular_formula_search(mz, parameters_filepath=None):
    
    mz = [mz]
    abundance = [1]
    rp, s2n = [[1],[1]]
    
    MSParameters.mass_spectrum.noise_threshold_method = 'relative_abundance'
    MSParameters.mass_spectrum.noise_threshold_absolute_abundance = 0 
    
    MSParameters.molecular_search.url_database = ''
    MSParameters.molecular_search.error_method = 'None'
    MSParameters.molecular_search.min_ppm_error  = -10
    MSParameters.molecular_search.max_ppm_error = 10
    MSParameters.molecular_search.mz_error_range = 1
    MSParameters.molecular_search.isProtonated = True 
    MSParameters.molecular_search.isRadical= False 
    MSParameters.molecular_search.isAdduct= False 

    usedatoms = {'C': (1,100) , 'H': (4,200), 'O': (0,10), 'N': (0,1), 'P': (0,1)}
    MSParameters.molecular_search.usedAtoms = usedatoms
    MSParameters.molecular_search.usedAtoms = usedatoms
    mass_spectrum_obj = ms_from_array_centroid(mz, abundance, rp, s2n, 'single mf search', polarity=1, auto_process=True)
    
    if parameters_filepath:
        
        parameter_from_json.load_and_set_parameters_ms(mass_spectrum_obj, parameters_path=parameters_filepath)

    mass_spectrum_obj.settings.noise_threshold_method = 'relative threshold'
    mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = False
    mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = 10
    mass_spectrum_obj.molecular_search_settings.use_isotopologue_filter = False

    print('Searching for molecular formulas within %.3f and %.3f ppm' % (mass_spectrum_obj.molecular_search_settings.min_ppm_error, mass_spectrum_obj.molecular_search_settings.max_ppm_error))

    SearchMolecularFormulas(mass_spectrum_obj, find_isotopologues=True).run_worker_ms_peaks([mass_spectrum_obj[0]])
    
    ms_peak = mass_spectrum_obj[0]
    
    if ms_peak:
        
        header = ['Molecular Formula',  'Calculated m/z', 'Mass Error', 'DBE', 'Ion Type']
        
        results = []
        
        for formula in ms_peak:
            
            results.append([formula.string, formula.mz_calc, formula.mz_error, formula.dbe, formula.ion_type])
              
        pprint.pprint(results)
        
        
    else:        
        
        print("Could not find a possible molecular formula match for the m/z %.5f" % mz[0])
        
if __name__ == "__main__":
    
    run_molecular_formula_search(760.58156938877)        