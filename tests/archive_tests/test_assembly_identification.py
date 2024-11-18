
import os
import sys
sys.path.append(".")
import pytest

from corems.encapsulation.factory.parameters import MSParameters
from corems.mass_spectrum.calc.CalibrationCalc import FreqDomain_Calibration
from corems.mass_spectrum.output.export import HighResMassSpecExport 
from corems.molecular_id.search.findOxygenPeaks import FindOxygenPeaks
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.molecular_id.calc.ClusterFilter import ClusteringFilter
from corems.transient.input.brukerSolarix import ReadBrukerSolarix

def calibrate(mass_spectrum_obj):
    
    find_formula_thread = FindOxygenPeaks(mass_spectrum_obj)
    find_formula_thread.run()
    mspeaks_results = find_formula_thread.get_list_found_peaks()
    
    #calibrate = FreqDomain_Calibration(mass_spectrum_obj, mspeaks_results)
    #calibrate.ledford_calibration()
    #mass_spectrum_obj.clear_molecular_formulas()

def filter_by_resolving_power():

    mass_spectrum_obj.filter_by_max_resolving_power(15, T)
    mass_spectrum_obj.filter_by_min_resolving_power(1, T)

def filter_by_kendrick():

    kendrick_base =  {'C':1,'H':2,'O':1}   
    mass_spectrum_obj.change_kendrick_base_all_mspeaks(kendrick_base)
        
        # needs to be wrapped inside the mass_spec class
    ClusteringFilter().filter_kendrick(mass_spectrum_obj)

def assign_mf_pox(mass_spectrum_obj):
    
    MSParameters.molecular_search.usedAtoms['O'] = (4, 20)
    MSParameters.molecular_search.usedAtoms['N'] = (0, 0)
    MSParameters.molecular_search.usedAtoms['S'] = (0, 0)
    MSParameters.molecular_search.usedAtoms['Cl'] = (0,0)
    MSParameters.molecular_search.usedAtoms['P'] = (1, 1)
      
    MSParameters.molecular_search.isProtonated = True
    MSParameters.molecular_search.isRadical = True
    MSParameters.molecular_search.isAdduct = True

    SearchMolecularFormulas(mass_spectrum_obj, first_hit=True).run_worker_mass_spectrum()

def assign_mf_nsox(mass_spectrum_obj):
    
    #print(len(mass_spectrum_obj), 'before kendrick filter')
    filter_by_resolving_power()
    #print(len(mass_spectrum_obj), 'after kendrick filter')
    #print(len(mass_spectrum_obj), 'after resolving power filter')

    MSParameters.molecular_search.usedAtoms['O'] = (4, 20)
    MSParameters.molecular_search.usedAtoms['N'] = (1, 3)
    MSParameters.molecular_search.usedAtoms['S'] = (1, 5)
    MSParameters.molecular_search.usedAtoms['Cl'] = (0, 0)
        
    
    MSParameters.molecular_search.min_dbe = 0
    MSParameters.molecular_search.max_dbe = 36
    
    MSParameters.molecular_search.isProtonated = True
    MSParameters.molecular_search.isRadical = True
    MSParameters.molecular_search.isAdduct = True

    SearchMolecularFormulas(mass_spectrum_obj, first_hit=True).run_worker_mass_spectrum(mass_spectrum_obj,)

def assign_mf_sox(mass_spectrum_obj):
    
    MSParameters.molecular_search.usedAtoms['O'] = (1, 10)
    MSParameters.molecular_search.usedAtoms['N'] = (0, 0)
    MSParameters.molecular_search.usedAtoms['S'] = (0, 3)
    MSParameters.molecular_search.usedAtoms['Cl'] = (0,0)
    MSParameters.molecular_search.usedAtoms['P'] = (0, 0)
   
    
    MSParameters.molecular_search.isProtonated = True
    MSParameters.molecular_search.isRadical = True
    MSParameters.molecular_search.isAdduct = True

    assignOx = OxygenPriorityAssignment(mass_spectrum_obj)
    assignOx.create_data_base()

    filter_by_kendrick()
    
    filter_by_resolving_power()
    
    assignOx.run()
    
    ClusteringFilter().remove_assignment_by_mass_error(mass_spectrum_obj)
    #assignOx.start()
    #assignOx.join()

def assign_mf_ox(mass_spectrum_obj):
    
    MSParameters.molecular_search.usedAtoms['O'] = (1, 10)
    MSParameters.molecular_search.usedAtoms['N'] = (0, 0)
    MSParameters.molecular_search.usedAtoms['S'] = (0, 0)
    MSParameters.molecular_search.usedAtoms['Cl'] = (0,0)
    MSParameters.molecular_search.usedAtoms['P'] = (0, 0)
    
   
    MSParameters.molecular_search.isProtonated = True
    MSParameters.molecular_search.isRadical = True
    MSParameters.molecular_search.isAdduct = True

    assignOx = OxygenPriorityAssignment(mass_spectrum_obj)
    assignOx.create_data_base()

    filter_by_kendrick()
    
    filter_by_resolving_power()
    
    assignOx.run()
    
    #ClusteringFilter().remove_assignment_by_mass_error(mass_spectrum_obj)
    #assignOx.start()
    #assignOx.join()

def assign_mf_nox(mass_spectrum_obj):
    
    MSParameters.molecular_search.usedAtoms['O'] = (1, 10)
    MSParameters.molecular_search.usedAtoms['N'] = (0, 4)
    MSParameters.molecular_search.usedAtoms['S'] = (0, 0)
    MSParameters.molecular_search.usedAtoms['Cl'] = (0, 0)
    MSParameters.molecular_search.usedAtoms['P'] = (0, 0)
  
    MSParameters.molecular_search.isProtonated = True
    MSParameters.molecular_search.isRadical = True
    MSParameters.molecular_search.isAdduct = True

    assignOx = OxygenPriorityAssignment(mass_spectrum_obj)
    assignOx.create_data_base()

    filter_by_kendrick()
    
    filter_by_resolving_power()
    
    assignOx.run()
    
    ClusteringFilter().remove_assignment_by_mass_error(mass_spectrum_obj)
    #assignOx.start()
    #assignOx.join()

def search_hc(mass_spectrum_obj):

    #print(len(mass_spectrum_obj), 'before kendrick filter')
    filter_by_resolving_power()
    #print(len(mass_spectrum_obj), 'after kendrick filter')
    #print(len(mass_spectrum_obj), 'after resolving power filter')

    MSParameters.molecular_search.usedAtoms['O'] = (0,0)
    MSParameters.molecular_search.usedAtoms['N'] = (0, 0)
    MSParameters.molecular_search.usedAtoms['S'] = (0, 0)
    MSParameters.molecular_search.usedAtoms['Cl'] = (0, 0)
    #MSParameters.molecular_search.usedAtoms['F'] = (0, 1)
    #MSParameters.molecular_search.usedAtoms['P'] = (0, 0)
    
    
    MSParameters.molecular_search.min_dbe = 0
    MSParameters.molecular_search.max_dbe = 50
    
    MSParameters.molecular_search.isProtonated = True
    MSParameters.molecular_search.isRadical = True
    MSParameters.molecular_search.isAdduct = True

    SearchMolecularFormulas(mass_spectrum_obj, first_hit=True).run_worker_mass_spectrum()

def search_nx(mass_spectrum_obj):

    #print(len(mass_spectrum_obj), 'before kendrick filter')
    filter_by_resolving_power()
    #print(len(mass_spectrum_obj), 'after kendrick filter')
    #print(len(mass_spectrum_obj), 'after resolving power filter')

    MSParameters.molecular_search.usedAtoms['O'] = (0,0)
    MSParameters.molecular_search.usedAtoms['N'] = (1, 3)
    MSParameters.molecular_search.usedAtoms['S'] = (0, 0)
    MSParameters.molecular_search.usedAtoms['Cl'] = (0, 0)
    #MSParameters.molecular_search.usedAtoms['F'] = (0, 1)
    #MSParameters.molecular_search.usedAtoms['P'] = (0, 0)
    
    
    MSParameters.molecular_search.min_dbe = 0
    MSParameters.molecular_search.max_dbe = 50
    
    MSParameters.molecular_search.isProtonated = True
    MSParameters.molecular_search.isRadical = True
    MSParameters.molecular_search.isAdduct = True

    SearchMolecularFormulas(mass_spectrum_obj,first_hit=True).run_worker_mass_spectrum()

def search_mf(mass_spectrum_obj):

    #print(len(mass_spectrum_obj), 'before kendrick filter')
    filter_by_resolving_power()
    #print(len(mass_spectrum_obj), 'after kendrick filter')
    #print(len(mass_spectrum_obj), 'after resolving power filter')

    MSParameters.molecular_search.usedAtoms['O'] = (0,0)
    MSParameters.molecular_search.usedAtoms['N'] = (0, 0)
    MSParameters.molecular_search.usedAtoms['S'] = (0, 0)
    MSParameters.molecular_search.usedAtoms['Cl'] = (0, 0)
    #MSParameters.molecular_search.usedAtoms['F'] = (0, 1)
    #MSParameters.molecular_search.usedAtoms['P'] = (0, 0)
    
    
    MSParameters.molecular_search.min_dbe = 0
    MSParameters.molecular_search.max_dbe = 50
    
    MSParameters.molecular_search.isProtonated = True
    MSParameters.molecular_search.isRadical = True
    MSParameters.molecular_search.isAdduct = True

    SearchMolecularFormulas(mass_spectrum_obj, first_hit=True).run_worker_mass_spectrum()

def search_ox(mass_spectrum_obj):

    filter_by_resolving_power()
   
    MSParameters.molecular_search.usedAtoms['O'] = (1,10)
    MSParameters.molecular_search.usedAtoms['N'] = (0, 0)
    MSParameters.molecular_search.usedAtoms['S'] = (0, 0)
    MSParameters.molecular_search.usedAtoms['Cl'] = (0, 1)
    #MSParameters.molecular_search.usedAtoms['F'] = (0, 1)
    #MSParameters.molecular_search.usedAtoms['P'] = (0, 0)
    
    SearchMolecularFormulas(mass_spectrum_obj, first_hit=True).run_worker_mass_spectrum()

def search_sx(mass_spectrum_obj):

    #print(len(mass_spectrum_obj), 'before kendrick filter')
    filter_by_resolving_power()
    #print(len(mass_spectrum_obj), 'after kendrick filter')
    #print(len(mass_spectrum_obj), 'after resolving power filter')

    MSParameters.molecular_search.usedAtoms['O'] = (0,0)
    MSParameters.molecular_search.usedAtoms['N'] = (0, 0)
    MSParameters.molecular_search.usedAtoms['S'] = (1, 3)
    MSParameters.molecular_search.usedAtoms['Cl'] = (0, 0)
    #MSParameters.molecular_search.usedAtoms['F'] = (0, 1)
    #MSParameters.molecular_search.usedAtoms['P'] = (0, 0)
    
    
    MSParameters.molecular_search.min_dbe = 0
    MSParameters.molecular_search.max_dbe = 50
  
    SearchMolecularFormulas(mass_spectrum_obj, first_hit=True).run_worker_mass_spectrum()

def search_sox(mass_spectrum_obj):

    filter_by_resolving_power()
    
    MSParameters.molecular_search.usedAtoms['O'] = (1, 10)
    MSParameters.molecular_search.usedAtoms['N'] = (0, 0)
    MSParameters.molecular_search.usedAtoms['S'] = (1, 3)
    MSParameters.molecular_search.usedAtoms['Cl'] = (0, 0)
    
    SearchMolecularFormulas(mass_spectrum_obj, first_hit=True).run_worker_mass_spectrum()

def search_nox(mass_spectrum_obj):

    
    filter_by_resolving_power()
    
    MSParameters.molecular_search.usedAtoms['O'] = (1, 10)
    MSParameters.molecular_search.usedAtoms['N'] = (1, 3)
    MSParameters.molecular_search.usedAtoms['S'] = (0, 0)
    MSParameters.molecular_search.usedAtoms['Cl'] = (0, 0)
   
    SearchMolecularFormulas(mass_spectrum_obj, first_hit=True).run_worker_mass_spectrum()
    
def search_nsox(mass_spectrum_obj):

    filter_by_resolving_power()
    
    MSParameters.molecular_search.usedAtoms['O'] = (1, 10)
    MSParameters.molecular_search.usedAtoms['N'] = (1, 3)
    MSParameters.molecular_search.usedAtoms['S'] = (1, 3)
    MSParameters.molecular_search.usedAtoms['Cl'] = (0, 0)
    
    MSParameters.molecular_search.min_dbe = 0
    MSParameters.molecular_search.max_dbe = 50
    
    MSParameters.molecular_search.isProtonated = True
    MSParameters.molecular_search.isRadical = True
    MSParameters.molecular_search.isAdduct = True

    SearchMolecularFormulas(mass_spectrum_obj, first_hit=True).run_worker_mass_spectrum()


def plot_error_distribution():

    for mspeak in mass_spectrum_obj:
        
        if mspeak:
            
            for molecular_formula in mspeak:
                if not molecular_formula.is_isotopologue:
                        
                    #off_set +=  0.1
                    molecular_formula.mz_error
                    pyplot.plot(mspeak.mz_exp, molecular_formula._assignment_mass_error, "o")

    pyplot.ylabel("m/z Error (ppm)")
    pyplot.xlabel('m/z')
    pyplot.savefig(mass_spectrum_obj.filename+'/'+"_error_dist"+'.png')        

def plot_resolving_power():

    pyplot.scatter(mass_spectrum_obj.mz_exp, mass_spectrum_obj.resolving_power/1000,
                                         s=mass_spectrum_obj.signal_to_noise, 
                                         cmap='seismic')

    pyplot.plot(mass_spectrum_obj.mz_exp, mass_spectrum_obj.resolving_power_calc(B, T, "low")/1000,  c='r')

    #pyplot.show()
    #pyplot.savefig(mass_spectrum_obj.filename+'/'+"_error_dist"+'.png')

def percent_assigned_peaks():
    
    total = len(mass_spectrum_obj)
    total_assigned = 0
    sum_assigned_abun = 0
    sum_abundance = sum(mass_spectrum_obj.abundance)
    
    for mspeak in mass_spectrum_obj:
        if mspeak: 
            total_assigned +=1
            sum_assigned_abun +=  mspeak.abundance
    
    print("%.2f of peaks assigned %.2f relative abundance"% (total_assigned/total*100, sum_assigned_abun/sum_abundance*100))

def plot_mass_spectrum():
    
    pyplot.plot(mass_spectrum_obj.mz_exp_profile, mass_spectrum_obj.abundance_profile)
    
    for mspeak in mass_spectrum_obj:
        
        if  mspeak:
            #off_set = 0 
            #for molecular_formula in mspeak:
            pyplot.plot(mspeak.mz_exp, mspeak.abundance, 'o', c='g')                       

                #if not molecular_formula.is_isotopologue:
                #    if molecular_formula.mspeak_indexes_isotopologues:
                        #pyplot.annotate(molecular_formula.string, (mspeak.mz_exp + off_set, mspeak.abundance ))
                        
                #        off_set +=  0.1
                    #if molecular_formula['O'] == o:
                    #    if  not molecular_formula.is_isotopologue:
                    #        pyplot.plot(molecular_formula['C'], molecular_formula.dbe, "o",   color=colors[molecular_formula['O']])
                    #        pyplot.plot(molecular_formula['C'], molecular_formula.dbe, "o",   color=colors[molecular_formula['O']])
                    #        pyplot.annotate(molecular_formula.class_label, (molecular_formula['C']+0.5, molecular_formula.dbe+0.5))

        else:
            
            continue
            #pyplot.plot(mspeak.mz_exp, mspeak.abundance, 'o', c='r')

    #pyplot.show()                

def plot_c_dbe_classes(df):
        
        print(df.columns)
        for h_class in df['Heteroatom Class'].unique():
            #print(df[df['Heteroatom Class']==h_class])
            carbon_number = df[df['Heteroatom Class']==h_class]['C']
            dbe = df[df['Heteroatom Class']==h_class]['DBE']
            abun  = df[df['Heteroatom Class']==h_class]['Measured Abundance']
            pyplot.scatter(carbon_number, dbe, c=abun)
            pyplot.xlabel("Carbon Number")
            pyplot.ylabel("DBE")
            pyplot.savefig(mass_spectrum_obj.filename+'/'+h_class+'.png')
            pyplot.close()

if __name__ == "__main__":
    pass
    
    #TODO include search for Na/H exchange for carboxylic acids  (needs to load from close shell options)
    # i.e [M-H + Na + Cl]-
    from matplotlib import colors as mcolors
    from matplotlib import pyplot
    import matplotlib
    from pathlib import Path
    #matplotlib.use('Agg')

    #file_name = "20190911_Kew_APPI_Elliot_DRY_IAT100ms_50ulmin_000001.d"
    #file_name = "20190912_Kew_LDI_Elliot_DRY_C4_0_C4_000001.d"
    #file_name= "20190912_kew_ldi_elliot_spe_f4_anchorchip_0_f4_000001.d"
    #file_name = "20190912_Kew_LDI_Elliot_DRY_C4_AnchorChip_0_C4_000001.d"
    #file_name = "20190911_Kew_ESI_Elliot_DRY_IAT100ms_000001.d"
    #file_name = "20190911_Kew_ESI_Elliot_SPE_IAT100ms_000002.d"
    #file_name = "20190912_Kew_LDI_Elliot_Whole_B15_1500um_0_B15_000001.d"
    
    file_location = Path.cwd() / "tests/tests_data/ftms/ESI_NEG_SRFA.d"
    
    bruker_reader = ReadBrukerSolarix(file_location)

    bruker_transient = bruker_reader.get_transient()
    
    T = bruker_transient.transient_time
    
    print(T, "T")

    B = 15

    mass_spectrum_obj = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=True)

    #mass_spectrum_obj.plot_profile_and_noise_threshold()
    
    MolecularFormulaSearchSettings = MolecularFormulaSearchSettings()
    
    #plot_resolving_power()

    #MSParameters.molecular_search.error_method = 'None'
    #MSParameters.molecular_search.min_ppm_error  = -10
    #MSParameters.molecular_search.max_ppm_error = 10

    calibrate(mass_spectrum_obj)

    plot_mass_spectrum()

    #MSParameters.molecular_search.error_method = 'None'
    #MSParameters.molecular_search.min_ppm_error  = -1
    #MSParameters.molecular_search.max_ppm_error = 1
    #MSParameters.molecular_search.mz_error_range = 1
    #MSParameters.molecular_search.mz_error_average = 0
    #MSParameters.molecular_search.min_abun_error = -30 # percentage
    #MSParameters.molecular_search.max_abun_error = 70 # percentage

    #assign_mf_ox(mass_spectrum_obj)

    #assign_mf_nox(mass_spectrum_obj)

    #assign_mf_sox(mass_spectrum_obj)

    #assign_mf_nsox(mass_spectrum_obj)

    #assign_mf_pox(mass_spectrum_obj)

    #plot_mass_spectrum()
    
    #print('what???', MSParameters.molecular_search.min_mz, MSParameters.molecular_search.max_mz)

    #search_hc(mass_spectrum_obj)

    #search_sx(mass_spectrum_obj)

    #search_nx(mass_spectrum_obj)

    #search_ox(mass_spectrum_obj)

    #search_nox(mass_spectrum_obj)

    #search_sox(mass_spectrum_obj)
    
    #search_nsox(mass_spectrum_obj)

    percent_assigned_peaks()

    #plot_mass_spectrum()
    
    #if not os.path.isdir(mass_spectrum_obj.filename):
    #    os.mkdir(mass_spectrum_obj.filename)
    
    #HighResMassSpecExport(mass_spectrum_obj.filename+'/'+mass_spectrum_obj.filename, mass_spectrum_obj, 'excel').start()
    
    #df = HighResMassSpecExport(mass_spectrum_obj.filename, mass_spectrum_obj, 'excel').get_pandas_df()

    #plot_c_dbe_classes(df)

    #plot_error_distribution()
