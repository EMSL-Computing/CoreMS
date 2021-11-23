from corems.transient.input.brukerSolarix import ReadBrukerSolarix
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.mass_spectrum.output.export import HighResMassSpecExport
from matplotlib import pyplot

file_path= 'tests/tests_data/ftms/ESI_NEG_SRFA.d'

#Bruker Solarix class reader
bruker_reader = ReadBrukerSolarix(file_path)

#access the transient object
bruker_transient_obj = bruker_reader.get_transient()

#calculates the transient duration time
T =  bruker_transient_obj.transient_time

#access the mass spectrum object
mass_spectrum_obj = bruker_transient_obj.get_mass_spectrum(plot_result=False, auto_process=True)

# - search monoisotopic molecular formulas for all mass spectral peaks
# - calculate fine isotopic structure based on monoisotopic molecular formulas found and current dynamic range
# - search molecular formulas of correspondent calculated isotopologues,
# - settings are stored at SearchConfig.json and can be changed directly on the file or inside the framework class

SearchMolecularFormulas(mass_spectrum_obj, first_hit=False).run_worker_mass_spectrum()

# iterate over mass spectral peaks objs
for mspeak in mass_spectrum_obj.sort_by_abundance():

    # returns true if there is at least one molecular formula associated
    # with the mass spectral peak
    # same as mspeak.is_assigned -- > bool
    if  mspeak:

        # get the molecular formula with the highest mass accuracy
        molecular_formula = mspeak.molecular_formula_lowest_error

        # plot mz and peak height, use mass_spectrum_obj.mz_exp to access all mz
        # and mass_spectrum_obj.mz_exp_profile to access mz with all available datapoints
        pyplot.plot(mspeak.mz_exp, mspeak.abundance, 'o', c='g')

        # iterate over all molecular formulae associated with the ms peaks obj
        for molecular_formula in mspeak:

            #check if the molecular formula is a isotopologue
            if molecular_formula.is_isotopologue:

                #access the molecular formula text representation
                print (molecular_formula.string)

                #get 13C atoms count
                print (molecular_formula['13C'])
    else:
        #get mz and peak height
        print(mspeak.mz_exp,mspeak.abundance)


#exporting data
mass_spectrum_obj.to_csv("filename")

mass_spectrum_obj.to_hdf("filename")
# save pandas Datarame to pickle
mass_spectrum_obj.to_pandas("filename")
# get pandas Dataframe
df = mass_spectrum_obj.to_dataframe()