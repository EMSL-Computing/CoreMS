# CoreMS Example Notebooks

This folder contains Jupyter notebooks demonstrating key workflows and capabilities of CoreMS for mass spectrometry data processing.

## Overview

CoreMS is a comprehensive framework for processing mass spectrometry data, including:
- Direct Infusion FT-ICR MS
- LC-MS/MS workflows with tandem MS
- GC-MS with metabolomics databases
- Advanced calibration and molecular formula assignment

## Notebooks

### 1. DirectInfusion_FTICR_Tutorial.ipynb
**Complete workflow for Direct Infusion FTICR-MS data**

Learn how to process high-resolution Fourier Transform Ion Cyclotron Resonance Mass Spectrometry data:
- Import and process Bruker transient files
- Apply noise thresholding and peak picking
- Perform mass calibration
- Assign molecular formulas
- Export and visualize results

**Data Format:** Bruker Solarix `.d` files  
**Recommended For:** High-resolution MS users, environmental/petroleum omics

---

### 2. GCMS_Tutorial.ipynb
**Automated compound identification using spectral matching**

Process low resolution GC-MS data with automated compound identification:
- Load ANDI NetCDF format GC-MS data
- Process chromatograms
- Perform retention index calibration using FAMES standards
- Match spectra against MetabRef library for compound identification
- Export results in multiple formats

**Data Format:** ANDI NetCDF `.cdf` files  
**Recommended For:** GC-MS metabolomics users  
**Note:** Uses MetabRef database interface

---

### 3. LCMS_Tutorial.ipynb
**Tandem mass spectrometry workflow for compound identification**

Process liquid chromatography mass spectrometry data with MS2 analysis:
- Load LC-MS data from Thermo RAW files
- Configure MS parameters for MS1 and MS2 data
- Process mass spectra and extract chromatographic features
- Perform molecular formula searches on MS1 data
- Extract and associate MS2 fragmentation data
- Annotate compounds using MS2 spectral matching
- Export and visualize results

**Data Format:** Thermo `.raw` files  
**Recommended For:** LC-MS users, metabolomics, lipidomics, DDA workflows

---

### 4. Mass_Recalibration_Tutorial.ipynb
**Improving mass accuracy through calibration**

Master mass recalibration techniques for high-resolution data:
- Manual recalibration with visual inspection
- Segmented recalibration for non-linear errors
- Automatic calibration using reference mass lists
- Assignment-based automatic recalibration

**Recommended For:** Users requiring high mass accuracy, FT-ICR MS users

---

### 5. Noise_Thresholding_Methods.ipynb
**Comparing noise filtering approaches**

Compare different noise thresholding methods for peak detection:
- Log-normal distribution method (recommended)
- Local minima-based thresholding
- Signal-to-noise ratio thresholds
- Relative and absolute abundance thresholds

**Recommended For:** Method optimization, quality control

---

### 6. Setting_MSParameters.ipynb
**Configuring global parameters in CoreMS**

Understand how to control data processing behavior:
- Set global MSParameters before data import
- Configure noise thresholding methods
- Adjust peak picking parameters
- Reset parameters to defaults

**Recommended For:** All users, beginners
