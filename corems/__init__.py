__version__ = '4.2.1.beta'
__doc__ = '''

CoreMS - a powerful framework for mass spectrometry data processing and analysis of small molecules
=====================================================================

**CoreMS** CoreMS is a Python package is a high-level building block for mass spectrometry data processing and
 software development that provides fast, flexible, and expressive data structures. The package is designed to 
 make working with mass spectrometry software and data both easy and intuitive by removing the steep learning curve
 needed to accurately implement signal processing, data curation and annotation. All modules are designed using an 
 intuitive mass spectrometric hierarchical class structure, using object-oriented programming paradigm. Each module 
 contains classes to store the data, expose pertinent functions, and calculate conventional mass spectrometric parameters. 


Main Features
-------------
Here are just a few of the things that CoreMS does well:

    Data input formats

    - Bruker Solarix CompassXtract
    - Bruker Solarix transients, ser and fid (FT and magnitude mode)
    - ThermoFisher raw file
    - Spectroswiss Signal booster data-acquisition station HFD 5 
    - Midas (.dat) from MagLab ICR data-acquisition station (FT-magnitude mode)
    - Mass list in Profile and Centroid Mode (include all delimiters types and Excel)
    - CoreMS exported processed mass list files (Excel, csv, txt, etc)
    - Panda dataframe, processed or unprocessed

    Data output formats

    - Text Files (csv, tab separated txt, etc)
    - Microsoft Excel (.xlsx)
    - Hierarchical Data Format (.h5) (`TODO`)
    - Pandas data frame (can be saved using pickle, h5, etc)

    Data structure type

    - GC-MS     (`TODO`)
    - LC-MS
    - IMS-MS    (`TODO`)
    - LC-IMS-MS (`TODO`)
    - Transient
    - Mass Spectra
    - Mass Spectrum

    Available features

        Signal Processing

        - Apodization, Zerofilling, and Magnitude mode FT
        - Manual and automatic noise threshold calculation
        - Peak picking apex quadratic fitting
        - Resolving Power calculation

        Calibration

        - Frequency and m/z domain calibration functions:
        - ledford equation [ref]
        - linear equation
        - quadratic equation
        - Automatic search most abundant **Ox** homologue series
        - step fit ('walking calibration") based on the ledford equation [ref]

        Molecular formulae search and assignment

        - Automatic local or external database generation
        - Automatic molecular formulae assignments for ESI(-) MS for natural organic matter analysis
        - Automatic fine isotopic structure calculation and search
        - Flexible Kendrick normalization base
        - Kendrick filter using density-based clustering

        Mass spectrum simulations

        - Peak shape (Lorentz and Gaussian)
        - Mass error distribution
        - ICR Resolving Power based on magnetic field (B), and transient time(T)
'''
