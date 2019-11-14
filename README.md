# CoreMS

**CoreMS** is a comprehensive mass spectrometry framework for software development and data analysis of small molecules analysis.

- reproducible pipeline
- logical mass spectrometric data structure
- self-containing data and metadata storage
- modern molecular formulae assigment algorithms
- dinamic molecular search space database search and generator

## Current Version

### `0.1.2.beta`

[![pipeline status](https://gitlab.pnnl.gov/mass-spectrometry/corems/badges/master/pipeline.svg)](https://gitlab.pnnl.gov/corilo/corems/commits/master) [![coverage report](https://gitlab.pnnl.gov/mass-spectrometry/corems/badges/master/coverage.svg)](https://gitlab.pnnl.gov/corilo/corems/commits/master)

### Data input formats

- Bruker Solarix (CompassXtract)
- Bruker Solarix transients, ser and fid (FT magnitude mode only)
- ThermoFisher (.raw)
- Spectroswiss signal booster data-acquisition station (.hdf5)
- MagLab ICR data-acquisition station (FT and magnitude mode) (.dat)
- Generic mass list in profile and centroid mde (include all delimiters types and Excel formats)
- CoreMS exported processed mass list files(excel, .csv, .txt, pandas dataframe as .pkl)
- CoreMS self-containing Hierarchical Data Format (.hdf5)
- Pandas dataframe

### Data output formats

- Pandas data frame (can be saved using pickle, h5, etc)
- Text Files (.csv, tab separated .txt, etc)
- Microsoft Excel (xlsx)
- Automatic JSON for metadata storage and reusage
- Self-containing Hierarchical Data Format (.hdf5) inclusing raw data and ime-series datapoint for processed datasets with all associated metadata stored as json attributes

### Data structure types

- LC-MS
- IMS-MS (`TODO`)
- LC-IMS-MS (`TODO`)
- Collections (`TODO`)
- Transient
- Mass Spectra
- Mass Spectrum
- Mass Spectral Peak
- Molecular Formula
- Molecular Structure (`TODO`)

## Available features

### Signal Processing

- Apodization, Zerofilling, and Magnitude mode FT
- Manual and automatic noise threshold calculation
- Peak picking using apex quadratic fitting
- Experimental resolving power calculation

### Calibration

- Frequency and m/z domain calibration functions:
- Ledford equation [ref]
- Linear equation
- Quadratic equation
- Automatic search most abundant **Ox** homologue serie
- Step fit ('walking calibration") based on the ledford equation [ref]

### Molecular formulae search and assignment

- Automatic local (SQLite) or external (MongoDB or Postgres) database check, generation, and search
- Automatic molecular formulae assignments algorithm for ESI(-) MS for natural organic matter analysis
- Automatic fine isotopic structure calculation and search for all isotopes
- Flexible Kendrick normalization base
- Kendrick filter using density-based clustering

### Mass spectrum simulations

- Peak shape (Lorentz and Gaussian)
- Mass error distribution(`TODO`)
- Calculated ICR Resolving Power based on magnetic field (B), and transient time(T)

## Jupyter-CoreMS

If you don't have docker installed, the easiest way is to [install docker for desktop](https://hub.docker.com/?overlay=onboarding)

- Open a terminal and run:

    ```bash
    docker run -v home/user/yourdir -p 8888:8888 gitlab.pnnl.gov:4567/mass-spectrometry/corems:latest
    ```

- In your browser, open the URL address provided in the terminal: `http://127.0.1.2:8888/?token=<token>.`

- Open the CoreMS-Tutorial.ipynb and follow the code

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
exportMS.output_type = 'hdf5' #csv, txt, #pandas
#save the file
exportMS.save()
```
