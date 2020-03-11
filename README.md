# CoreMS

**CoreMS** is a comprehensive mass spectrometry framework for software development and data analysis of small molecules analysis.

CoreMS is a Python package of high-level building blocks for mass spectrometry data processing and software development that provides fast, flexible, and expressive data structures. The package is designed to make working with mass spectrometry software and data both easy and intuitive by removing the steep learning curve needed to accurately implement signal processing, data curation and annotation. All modules are designed using an intuitive mass spectrometric hierarchical class structure, using object-oriented programming paradigm. Each module contains classes to store the data, expose pertinent functions, and calculate conventional mass spectrometric parameters.

- reproducible pipeline
- logical mass spectrometric data structure
- self-containing data and metadata storage
- modern molecular formulae assignment algorithms
- dynamic molecular search space database search and generator

## Current Version

### `8.3.0.beta`

[![pipeline status](https://gitlab.pnnl.gov/mass-spectrometry/corems/badges/master/pipeline.svg)](https://gitlab.pnnl.gov/corilo/corems/commits/master) [![coverage report](https://gitlab.pnnl.gov/mass-spectrometry/corems/badges/master/coverage.svg)](https://gitlab.pnnl.gov/corilo/corems/commits/master)

### Data input formats

- Bruker Solarix (CompassXtract)
- Bruker Solarix transients, ser and fid (FT magnitude mode only)
- ThermoFisher (.raw)
- Spectroswiss signal booster data-acquisition station (.hdf5)
- MagLab ICR data-acquisition station (FT and magnitude mode) (.dat)
- ANDI NetCDF for GC-MS (.cdf)
- Generic mass list in profile and centroid mde (include all delimiters types and Excel formats)
- CoreMS exported processed mass list files(excel, .csv, .txt, pandas dataframe as .pkl)
- CoreMS self-containing Hierarchical Data Format (.hdf5)
- Pandas Dataframe

### Data output formats

- Pandas data frame (can be saved using pickle, h5, etc)
- Text Files (.csv, tab separated .txt, etc)
- Microsoft Excel (xlsx)
- Automatic JSON for metadata storage and reuse
- Self-containing Hierarchical Data Format (.hdf5) including raw data and ime-series data-point for processed data-sets with all associated metadata stored as json attributes

### Data structure types

- LC-MS
- GC-MS
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
- LedFord equation [ref]
- Linear equation
- Quadratic equation
- Automatic search most abundant **Ox** homologue series
- Step fit ('walking calibration") based on the LedFord equation [ref]

### Molecular formulae search and assignment

- Automatic local (SQLite) or external (MongoDB or PostgreSQL) database check, generation, and search
- Automatic molecular formulae assignments algorithm for ESI(-) MS for natural organic matter analysis
- Automatic fine isotopic structure calculation and search for all isotopes
- Flexible Kendrick normalization base
- Kendrick filter using density-based clustering
- Kendrick classification
- Heteroatoms classification and visualization

### Mass spectrum simulations

- Peak shape (Lorentz,  Gaussian, Voigt, and pseudo-Voigt)
- Peak fitting for peak shape definition
- Peak position in function of datapoints, signal to noise and resolving power (Lorentz and Gaussian)
- Prediction of mass error distribution
- Calculated ICR Resolving Power based on magnetic field (B), and transient time(T)

## Jupyter-CoreMS

A docker image containing a custom python distribution with all dependencies and a Jupyter server with notebook examples

If you don't have docker installed, the easiest way is to [install docker for desktop](https://hub.docker.com/?overlay=onboarding)

- Open a terminal and run:

    ```bash
    docker run --rm -v host/dir:/home/CoreMS/data: -p 8888:8888 gitlab.pnnl.gov:4567/mass-spectrometry/corems:latest
    ```

- In your browser, open the URL address provided in the terminal: `http://127.3.2.1:8888/?token=<token>.`

- Open the CoreMS-Tutorial.ipynb

## Basic example

```python
from corems.transient.input.brukerSolarix import ReadBrukerSolarix
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.mass_spectrum.output.export import MassSpecExport

file_path= 'neg_esi_srfa_1ppm_test.d'

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

SearchMolecularFormulas(first_hit=False).run_worker_mass_spectrum(mass_spectrum_obj)

# iterate over mass spectral peaks objs
for mspeak in mass_spectrum_obj.sort_by_abundance():

    # returns true if there is at least one molecular formula associated
    # with the mass spectral peak
    # same as mspeak.is_assigned -- > bool
    if  mspeak:

        # get the molecular formula with the lowest highest mass accuracy
        molecular_formula = mspeak.molecular_formula_lowest_error

        # plot mz and peak height, use mass_spectrum_obj.mz_exp to access all mz
        # and mass_spectrum_obj.mz_exp_profile to access mz with all available datapoints
        pyplot.plot(mspeak.mz_exp, mspeak.abundance, 'o', c='g')

        # iterate over all molecular formulae associated with the ms peaks obj
        for molecular_formula in mspeak:

            #check if the molecular formula is a isotopologue
            if molecular_formula.is_isotopologue:

                #access the molecular formula text representation
                print (molecular_formula.to_string)

                #get 13C atoms count
                print (molecular_formula['13C'])
    else:
        #get mz and peak height
        print(mspeak.mz_exp,mspeak.abundance)


#exporting data
mass_spectrum_obj.to_csv("filename")

mass_spectrum_obj.to_hdf5("filename")

mass_spectrum_obj.to_pandas("filename")
```
