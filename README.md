
# CoreMS

**CoreMS** is a comprehensive mass spectrometry framework for small molecules analysis.

## Current Version

### `6.2.0.alpha`

[![pipeline status](https://gitlab.pnnl.gov/mass-spectrometry/corems/badges/master/pipeline.svg)](https://gitlab.pnnl.gov/corilo/corems/commits/master) [![coverage report](https://gitlab.pnnl.gov/mass-spectrometry/corems/badges/master/coverage.svg)](https://gitlab.pnnl.gov/corilo/corems/commits/master)

### Data input formats

- Bruker Solarix ComprassXtract
- Bruker Solarix transients, ser and fid (FT and magnitude mode)
- ThermoFisher Raw
- Spectroswiss Signal booster data-acquisition station
- Midas (.dat) from MagLab ICR data-acquisition station (FT and magnitude mode)
- Mass list in Profile and Centroid Mode (include all delimiters types and Excel)
- CoreMS exported processed mass list files(Excel, csv, txt, hdf5)
- Panda dataframe

### Data output formats

- Text Files (csv, tab separated txt, etc)
- Microsoft Excel (xlsx)
- Hierarchical Data Format (.hdf5) with processing setting and class attributes
- Pandas data frame (can be saved using pickle, h5, etc)

### Data structure type

- LC-MS
- IMS-MS
- LC-IMS-MS (`TODO`)
- Transient
- Mass Spectra
- Mass Spectrum

## Available features

### Signal Processing

- Apodization, Zerofilling, and Magnitude mode FT
- Manual and automatic noise threshold calculation
- Peak picking apex quadratic fitting
- Resolving Power calculation

### Calibration

- Frequency and m/z domain calibration functions:
- ledford equation [ref]
- linear equation
- quadratic equation
- Automatic search most abundant **Ox** homologue serie
- step fit ('walking calibration") based on the ledford equation [ref]

### Molecular formulae search and assigment

- Automatic local (SQLite) or external (MongoDB or Postgres) database check, generation, and search
- Automatic molecular formulae assignments algorithm for ESI(-) MS for natural organic matter analysis
- Automatic fine isotopic structure calculation and search
- Flexible Kendrick normalization base
- Kendrick filter using density-based clustering

### Mass spectrum simulations

- Peak shape (Lorentz and Gaussian)
- Mass error distribution(`TODO`)
- ICR Resolving Power based on magnetic field (B), and transient time(T)

## Basic example

```python
from corems.transient.input.BrukerSolarix import ReadBrukerSolarix
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.mass_spectrum.output.export import MassSpecExport

file_path= 'neg_esi_srfa_1ppm_test.d'

#Bruker Solarix class reader
bruker_reader = ReadBrukerSolarix(file_path)

#access the transient object
bruker_transient_obj = bruker_reader.get_transient()

#calculates the transient duration time
T =  bruker_transient_obj.transient_time

#access the mass spectrum object
mass_spectrum_obj = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=True)

# - search monoisotopic molecular formulas for all mass spectral peaks
# - calculate fine isotopic structure based on monoisotopic molecular formulas found and current dynamic range
# - search molecular formulas of correspondent calculated isotopologues,
# - settings are stored at SearchConfig.json and can be changed directly on the file or inside the framework class

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

            #check if the molecular formula is a isotopologue
            if molecular_formula.is_isotopologue:

                #access the molecular formula text representation
                print (molecular_formula.to_string())

                #get 13C atoms count
                print (molecular_formula['13C'])
    else:
        #get mz and peak height
        print(mspeak.mz_exp,mspeak.abundance)

#create export class
exportMS = MassSpecExport('neg_esi_srfa_1ppm_test',  mass_spectrum_obj.filter_by_sn(4))
#get pandas dataframe
df_obj = exportMS.get_pandas_df()
#set file output type
exportMS.output_type = ’excel’ #csv, txt
#save the file
exportMS.save()
```
