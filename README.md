# CoreMS

**CoreMS** is a comprehensive mass spectrometry framework for software development and data analysis of small molecules analysis.

Data handling and software development for modern mass spectrometry (MS) is an interdisciplinary endeavor requiring skills in computational science and a deep understanding of MS. To enable scientific software development to keep pace with fast improvements in MS technology, we have developed a Python software framework named CoreMS. The goal of the framework is to provide a fundamental, high-level basis for working with all mass spectrometry data types, allowing custom workflows for data signal processing, annotation, and curation. The data structures were designed with an intuitive, mass spectrometric hierarchical structure, thus allowing organized and easy access to the data and calculations. Moreover, CoreMS supports direct access for almost all vendors’ data formats, allowing for the centralization and automation of all data processing workflows from the raw signal to data annotation and curation.

- reproducible pipeline
- logical mass spectrometric data structure
- self-containing data and metadata storage
- modern molecular formulae assignment algorithms
- dynamic molecular search space database search and generator

## Current Version

### `16.0.0.beta`

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

- Automatic local (SQLite) or external (PostgreSQL) database check, generation, and search
- Automatic molecular formulae assignments algorithm for ESI(-) MS for natural organic matter analysis
- Automatic fine isotopic structure calculation and search for all isotopes
- Flexible Kendrick normalization base
- Kendrick filter using density-based clustering
- Kendrick classification
- Heteroatoms classification and visualization

### Mass spectrum simulations

- Peak shape (Lorentz,  Gaussian, Voigt, and pseudo-Voigt)
- Peak fitting for peak shape definition
- Peak position in function of data points, signal to noise and resolving power (Lorentz and Gaussian)
- Prediction of mass error distribution
- Calculated ICR Resolving Power based on magnetic field (B), and transient time(T)

## CoreMS Installation
    
```bash
pip install corems
```

A external database is needed to run the molecular formula assignments workflow:  
```bash
docker-compose up -d
```
-  Database url = "postgres://coremsdb:coremsmolform@localhost:5432/molformula"  
    The database url is already set the default on the CoreMS parameter model class

To be able to open thermo file a installation of pythonnet is needed:
- Windows: 
    ```bash
    pip install pythonnet
    ```

- Mac and Linux:
    ```bash
    brew install mono
    pip install pythonnet   
    ```
 Another option is to run the docker stack that will start the CoreMS containers(see next section)

## Molecular Database and Jupyter Notebook

A docker container containing:
- A custom python distribution will all dependencies installed
- A Jupyter notebook server with workflow examples
- A PostgreSQL database for the molecular formulae assignment

If you don't have docker installed, the easiest way is to [install docker for desktop](https://hub.docker.com/?overlay=onboarding)

- Start the containers from the latest built docker image (easiest way): 

    ```bash
    docker-compose -f docker-compose-jupyter.yml up
    ```

- Build a new image with current changes: 

    1) Build the corems image:
        ```bash
        docker build -t corems:local .
        ```
    2. Start the database container:
        ```bash
        docker-compose up -d   
        ```
    3. Start the Jupyter Notebook:
        ```bash
        docker run --rm -v ./data:/home/CoreMS/data corems:local
        ```

- Open your browser, copy and past the URL address provided in the terminal: `http://localhost:8888/?token=<token>.`

- Open the CoreMS-Tutorial.ipynb

## Examples

More examples can be found under the directory docs/example

- Basic functionality example

```python
from corems.transient.input.brukerSolarix import ReadBrukerSolarix
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.mass_spectrum.output.export import HighResMassSpecExport

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

mass_spectrum_obj.to_hdf5("filename")
# save pandas Datarame to pickle
mass_spectrum_obj.to_pandas("filename")
# get pandas Dataframe
df = mass_spectrum_obj.to_dataframe("filename")
```
