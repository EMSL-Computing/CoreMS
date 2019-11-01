[![pipeline status](https://gitlab.pnnl.gov/mass-spectrometry/corems/badges/master/pipeline.svg)](https://gitlab.pnnl.gov/corilo/corems/commits/master)
[![coverage report](https://gitlab.pnnl.gov/mass-spectrometry/corems/badges/master/coverage.svg)](https://gitlab.pnnl.gov/corilo/corems/commits/master)

# CoreMS 
Mass Spectrometry ToolBox for Small Molecules Analysis

## Currrent Version 
[4.2.0.alpha]

# Basic usage:

```python

file_path= 'neg_esi_srfa_1ppm_test.d'

#Bruker Solarix class reader
bruker_reader = ReadBrukerSolarix(file_path) --> str
 
#access the transient object     
bruker_transient_obj = bruker_reader.get_transient()

#calculates the transient duration time     
T =  bruker_transient_obj.transient_time

#access the mass spectrum object      
mass_spectrum_obj = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=True)

# - search molecular formulas for monoisotopic mass spectral peaks
# - calculate fine isotopic structure  based on monoisotopic found 
# - search molecular formulas for correspondent isotopologues
# - settings are stored at SearchConfig.json
SearchMolecularFormulas(first_hit=False).run_worker_mass_spectrum(mass_spectrum_obj)

# iterate over mass spectral peaks objs
for mspeak in mass_spectrum_obj.sort_by_abundance():
    
    # returns true if there is at least one molecular formula associated
    # with the mass spectral peak
    # same as mspeak.is_assigned -- > bool
	if  mspeak:   
        # get the molecular formula with the lowest highest mass accuracy      
		molecular_formula = mspeak.molecular_formula_lowest_error
            # plot mz and peak height, use mass_spectrum_obj.mz_exp to access all mz
            # and mass_spectrum_obj.mz_exp_profile to access mz with all available datapoints    
		    pyplot.plot(mspeak.mz_exp, mspeak.abundance, 'o', c='g')                       
            # iterate over all molecular formulae associated with the ms peaks obj 
            for molecular_formula in mspeak:
                # check if the molecular formula is a isotopologue
                if molecular_formula.is_isotopologue:
                    # access the molecular formula text representation
                    print (molecular_formula.to_string()) 
                    # get 13C atoms count
                    print (molecular_formula[’13C’])  
      else:
        # get mz and peak height  
        print(mspeak.mz_exp,mspeak.abundance)       

# create export class
exportMS = MassSpecExport('neg_esi_srfa_1ppm_test',  mass_spectrum_obj.filter_by_sn())
# get pandas dataframe
df_obj = exportMS.get_pandas_df()
# set file output type
exportMS.output_type = ’excel’ #csv, txt
# save the file 
exportMS.save()

```