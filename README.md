![CoreMS Logo](https://github.com/EMSL-Computing/CoreMS/blob/master/docs/CoreMS.COLOR_small.png?raw=true)  

<div align="left">

<br>
<br>
<a href="https://doi.org/10.5281/zenodo.14009575"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.14009575.svg" alt="DOI"></a>
<a href="https://github.com/EMSL-Computing/CoreMS/actions/workflows/ci.yml"><img src="https://github.com/EMSL-Computing/CoreMS/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
<a href="https://pypi.org/project/CoreMS/"><img src="https://img.shields.io/pypi/v/CoreMS.svg" alt="PyPI"></a>
<a href="https://pypi.org/project/CoreMS/"><img src="https://img.shields.io/pypi/pyversions/CoreMS.svg" alt="Python versions"></a>
<br>
</div>

# Table of Contents  
- Introduction
  - [CoreMS](#CoreMS)  
  - [Current Version](#current-version)  
  - [Contact Information](#main-developers/contact )  
  - [Documentation](#documentation)
  - [Contribution Information](#contributing)
  - [Data Input](#data-input-formats)  
  - [Data Output](#data-output-formats)  
  - [Data Structure](#data-structure-types)  
  - [Features](#available-features)  
- Installation
  - [Installation](#corems-installation)  
  - [Thermo Raw File on Mac and Linux](#thermo-raw-file-access)  
- Execution     
  - [Building and Running the Docker Image](#docker-image)
  - [Example for FT-ICR Data Processing](#simple-script-example)  
  - [Jupyter Notebook Examples](examples/notebooks)  
- Sibling Projects    
    - [EnviroMS](https://github.com/EMSL-Computing/EnviroMS)  
    - [MetaMS](https://github.com/EMSL-Computing/MetaMS)  

***

# CoreMS

**CoreMS** is a comprehensive mass spectrometry framework for software development and data analysis of small molecules analysis.

Data handling and software development for modern mass spectrometry (MS) is an interdisciplinary endeavor requiring skills in computational science and a deep understanding of MS. To enable scientific software development to keep pace with fast improvements in MS technology, we have developed a Python software framework named CoreMS. The goal of the framework is to provide a fundamental, high-level basis for working with all mass spectrometry data types, allowing custom workflows for data signal processing, annotation, and curation. The data structures were designed with an intuitive, mass spectrometric hierarchical structure, thus allowing organized and easy access to the data and calculations. Moreover, CoreMS supports direct access for almost all vendors’ data formats, allowing for the centralization and automation of all data processing workflows from the raw signal to data annotation and curation.

CoreMS aims to provide 
- logical mass spectrometric data structure
- self-containing data and metadata storage
- modern molecular formulae assignment algorithms
- dynamic molecular search space database search and generator

***

## Current Version

 `3.11.0`

***

## Main Developers/Contact 
- [Yuri. E. Corilo](mailto:corilo@pnnl.gov)  
- [William Kew](mailto:william.kew@pnnl.gov)
- [Katherine Heal](mailto:katherine.heal@pnnl.gov)

***

## Documentation

API documentation can be found [here](https://emsl-computing.github.io/CoreMS/corems.html).

Overview slides can be found [here](https://github.com/EMSL-Computing/CoreMS/blob/master/examples/CoreMS-Overview.pdf).

***

## Contributing

As an open source project, CoreMS welcomes contributions of all forms. Before contributing, please see our [Dev Guide](./CONTRIBUTING.md)

***

## Data formats
### Data input formats

- Bruker Solarix (CompassXtract)
- Bruker Solarix transients, ser and fid (FT magnitude mode only)
- ThermoFisher (.raw)
- Spectroswiss signal booster data-acquisition station (.hdf5)
- MagLab ICR data-acquisition station (FT and magnitude mode) (.dat)
- ANDI NetCDF for GC-MS (.cdf)
- mzml for LC-MS (.mzml)
- Generic mass list in profile and centroid mde (include all delimiters types and Excel formats)
- CoreMS exported processed mass list files(excel, .csv, .txt, pandas dataframe as .pkl)
- CoreMS self-containing Hierarchical Data Format (.hdf5)
- Pandas Dataframe
- Support for cloud Storage using s3path.S3path

### Data output formats

- Pandas data frame (can be saved using pickle, h5, etc)
- Text Files (.csv, tab separated .txt, etc)
- Microsoft Excel (xlsx)
- Automatic JSON for metadata storage and reuse
- Self-containing Hierarchical Data Format (.hdf5) including raw data and time-series data-point for processed data-sets with all associated metadata stored as json attributes

### Data structure types

- LC-MS
- GC-MS
- Transient
- Mass Spectra
- Mass Spectrum
- Mass Spectral Peak
- Molecular Formula

***

## Available features

### FT-MS Signal Processing, Calibration, and Molecular Formula Search and Assignment

- Apodization, Zerofilling, and Magnitude mode FT
- Manual and automatic noise threshold calculation
- Peak picking using apex quadratic fitting
- Experimental resolving power calculation
- Frequency and m/z domain calibration functions:
- LedFord equation
- Linear equation
- Quadratic equation
- Automatic search most abundant **Ox** homologue series
- Automatic local (SQLite) or external (PostgreSQL) database check, generation, and search
- Automatic molecular formulae assignments algorithm for ESI(-) MS for natural organic matter analysis
- Automatic fine isotopic structure calculation and search for all isotopes
- Flexible Kendrick normalization base
- Kendrick filter using density-based clustering
- Kendrick classification
- Heteroatoms classification and visualization

### GC-MS Signal Processing, Calibration, and Compound Identification

- Baseline detection, subtraction, smoothing 
- m/z based Chromatogram Peak Deconvolution,
- Manual and automatic noise threshold calculation
- First and second derivatives peak picking methods
- Peak Area Calculation
- Retention Index Calibration
- Automatic local (SQLite) or external (MongoDB or PostgreSQL) database check, generation, and search
- Automatic molecular match algorithm with all spectral similarity methods 

### High Resolution Mass Spectrum Simulations

- Peak shape (Lorentz,  Gaussian, Voigt, and pseudo-Voigt)
- Peak fitting for peak shape definition
- Peak position in function of data points, signal to noise and resolving power (Lorentz and Gaussian)
- Prediction of mass error distribution
- Calculated ICR Resolving Power based on magnetic field (B), and transient time(T)

### LC-MS Signal Processing, Molecular Formula Search and Assignment, and Spectral Similarity Searches
See walkthrough in [this notebook](examples/notebooks/LCMS_Tutorial.ipynb)
- Two dimensional (m/z and retention time) peak picking using persistent homology
- Smoothing, cetroid detection, and integration of extracted ion chromatograms
- Peak shape metric calculations including half peak height, tailing factor, and dispersity index
- MS1 deconvolution of mass features
- Idenfitication of <sup>13</sup>C isotopes within the mass features
- Compatibility with molecular formula searching on MS1 or MS2 spectra
- Spectral search capability using entropy similarity

***

## Installation 
    
```bash
pip install corems
```

Corems requires **Python 3.9 or later** (including Python 3.13) and is compatible with **NumPy 2.x**, **pandas 2.x**, and **SQLAlchemy 2.x**.

To install with development and testing extras:

```bash
pip install "corems[dev]"
```

By default the molecular formula database will be generated using SQLite.

To use PostgreSQL the easiest way is to build a docker container:

```bash
docker-compose up -d
```

- Change the url_database on `MSParameters.molecular_search.url_database` to: `"postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp"`
- Set the env variable `COREMS_DATABASE_URL` to: `"postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp"`

### Thermo Raw File Access:

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

***

## Building and Running the CoreMS Docker Image <a name="docker-image"></a>

CoreMS provides a Dockerfile that packages the entire application (including .NET 8 runtime for Thermo .raw file support) into a self-contained image. This is useful for running CoreMS in a reproducible environment without installing dependencies on your host system.

### Prerequisites
- Docker installed and running on your system.
- The CoreMS repository cloned locally.
- Navigate to the root of the CoreMS repository before running any commands.

### Building the Docker Image

The Makefile provides convenience targets for building the image. The image is tagged with the current version from `.bumpversion.cfg`.

**On Linux/Windows (standard build):**
```bash
make build-image-local
```

**On macOS (cross-platform build for linux/amd64):**
```bash
make build-image-mac-local
```

This runs `docker build` with the `--platform linux/amd64` flag, which is necessary when building on Apple Silicon (M1/M2/M3) Macs to ensure compatibility.

Alternatively, you can build manually with:
```bash
docker build -t corems:<version> .
```
Replace `<version>` with your desired tag (e.g., `3.11.0`).

### What the Dockerfile Does

The Dockerfile performs the following steps:
1. Starts from a `python:3.13-slim` base image.
2. Installs the .NET 8 runtime (required for Thermo .raw file support via PythonNET).
3. Copies the CoreMS source code into the image.
4. Installs CoreMS and all its dependencies via `pip install .`.
5. Installs `pytest-xdist` and `pytest-cov` for running tests.
6. Removes build-time dependencies (gcc, python3-dev) to keep the image lean.

### Running the Docker Image

**On Linux/Windows:**
```bash
make image-run-local
```

**On macOS:**
```bash
make image-run-mac-local
```

This launches an interactive bash shell inside the container:
```bash
docker run -it corems:<version>
```

From within the container, you can import and use CoreMS directly:
```python
python3 -c "import corems; print(corems.__version__)"
```

### Mounting Data into the Container

To process your own data files, mount a local directory into the container:
```bash
docker run -it -v /path/to/your/data:/data corilo/corems:<version>
```
Your files will then be accessible at `/data` inside the container.

### Managing the PostgreSQL Database with Docker Compose

The `docker-compose.yml` file defines a PostgreSQL database service for CoreMS. The Makefile provides targets to manage it:

**Start the database:**
```bash
make db-up
```

**Stop the database:**
```bash
make db-down
```

**View database logs:**
```bash
make db-logs
```

These are equivalent to running `docker-compose up -d`, `docker-compose down`, and `docker-compose logs -f` respectively.

***

## Example for FT-ICR Data Processing

More examples can be found in the [examples/notebooks](examples/notebooks) directory

- Basic functionality example

```python
from corems.transient.input.brukerSolarix import ReadBrukerSolarix
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.mass_spectrum.output.export import HighResMassSpecExport
from matplotlib import pyplot

file_path= 'tests/tests_data/ftms/ESI_NEG_SRFA.d'

# Instatiate the Bruker Solarix reader with the filepath
bruker_reader = ReadBrukerSolarix(file_path)

# Use the reader to instatiate a transient object
bruker_transient_obj = bruker_reader.get_transient()

# Calculate the transient duration time
T =  bruker_transient_obj.transient_time

# Use the transient object to instatitate a mass spectrum object
mass_spectrum_obj = bruker_transient_obj.get_mass_spectrum(plot_result=False, auto_process=True)

# The following SearchMolecularFormulas function does the following
# - searches monoisotopic molecular formulas for all mass spectral peaks
# - calculates fine isotopic structure based on monoisotopic molecular formulas found and current dynamic range
# - searches molecular formulas of correspondent calculated isotopologues
# - settings are stored at SearchConfig.json and can be changed directly on the file or inside the framework class

SearchMolecularFormulas(mass_spectrum_obj, first_hit=False).run_worker_mass_spectrum()

# Iterate over mass spectral peaks objs within the mass_spectrum_obj
for mspeak in mass_spectrum_obj.sort_by_abundance():

    # If there is at least one molecular formula associated, mspeak returns True
    if  mspeak:

        # Get the molecular formula with the highest mass accuracy
        molecular_formula = mspeak.molecular_formula_lowest_error

        # Plot mz and peak height
        pyplot.plot(mspeak.mz_exp, mspeak.abundance, 'o', c='g')

        # Iterate over all molecular formulas associated with the ms peaks obj
        for molecular_formula in mspeak:

            # Check if the molecular formula is a isotopologue
            if molecular_formula.is_isotopologue:

                # Access the molecular formula text representation and print
                print (molecular_formula.string)

                # Get 13C atoms count
                print (molecular_formula['13C'])
    else:
        # Get mz and peak height
        print(mspeak.mz_exp,mspeak.abundance)

# Save data
## to a csv file
mass_spectrum_obj.to_csv("filename")
mass_spectrum_obj.to_hdf("filename")
# to pandas Datarame pickle
mass_spectrum_obj.to_pandas("filename")

# Extract data as a pandas Dataframe
df = mass_spectrum_obj.to_dataframe()
```

***

## UML Diagrams

UML (unified modeling language) diagrams for Direct Infusion FT-MS and GC-MS classes can be found [here](docs/uml).

***

## Citing CoreMS

If you use CoreMS in your work, please use the following citation:

Version [3.11.0 Release on GitHub](https://github.com/EMSL-Computing/CoreMS/releases/tag/v3.11.0), archived on Zenodo:  

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14009575.svg)](https://doi.org/10.5281/zenodo.14009575)

Yuri E. Corilo, William R. Kew, Lee Ann McCue, Katherine R . Heal, James C. Carr (2024, October 29). EMSL-Computing/CoreMS: CoreMS 3.0.0 (Version v3.0.0), as developed on Github. Zenodo. http://doi.org/10.5281/zenodo.14009575

```

***


This material was prepared as an account of work sponsored by an agency of the
United States Government.  Neither the United States Government nor the United
States Department of Energy, nor Battelle, nor any of their employees, nor any
jurisdiction or organization that has cooperated in the development of these
materials, makes any warranty, express or implied, or assumes any legal
liability or responsibility for the accuracy, completeness, or usefulness or
any information, apparatus, product, software, or process disclosed, or
represents that its use would not infringe privately owned rights.

Reference herein to any specific commercial product, process, or service by
trade name, trademark, manufacturer, or otherwise does not necessarily
constitute or imply its endorsement, recommendation, or favoring by the United
States Government or any agency thereof, or Battelle Memorial Institute. The
views and opinions of authors expressed herein do not necessarily state or
reflect those of the United States Government or any agency thereof.

                 PACIFIC NORTHWEST NATIONAL LABORATORY
                              operated by
                                BATTELLE
                                for the
                   UNITED STATES DEPARTMENT OF ENERGY
                    under Contract DE-AC05-76RL01830

