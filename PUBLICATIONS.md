# CoreMS Publications and Citations

This document tracks peer-reviewed publications, conference papers, and other scientific works that reference, use, or review CoreMS. The entries are organized by category to facilitate discovery and understanding of CoreMS applications and methodological contributions.

Current totals: 21 CoreMS citations tracked (Methods: 4, Application: 13, Additional: 4).

## Table of Contents

- [Methods & Evaluation Papers (4)](#methods--evaluation-papers-4)
- [Application & Use Case Papers (13)](#application--use-case-papers-13)
- [Additional Citations & References (4)](#additional-citations--references-4)

---

## Methods & Evaluation Papers (4)

Papers that present methodological contributions, algorithmic improvements, or comprehensive evaluation/benchmarking of CoreMS capabilities.

| No. | Authors | Title | Year | Journal/Conference | DOI/Link | Key Contributions |
|-----|---------|-------|------|-------------------|----------|------------------|
| 1 | Dewey, C., Corilo, Y., Kew, W., Boiteau, R.M. | Data Set Analysis to Reduce Uncertainty in Formula Assignments of Ultrahigh Resolution Mass Spectra | 2025 | Analytical Chemistry, Vol. 97, No. 25, pp. 13031-13039 | https://doi.org/10.1021/acs.analchem.4c06826 | Data set analysis strategy for improving formula assignment confidence, filtering false assignments, and mitigating bias in HRMS formula routines; demonstrated on 21 T FT-ICR MS oceanographic data |
| 2 | Degnan, D.J., Flores, J.E., Brayfindley, E.R., Paurus, V.L., Webb-Robertson, B.J.M., Clendinen, C.S., Bramer, L.M. | Characterizing Families of Spectral Similarity Scores and Their Use Cases for Gas Chromatography-Mass Spectrometry Small Molecule Identification | 2023 | Metabolites, Vol. 13, No. 10, 1101 | https://doi.org/10.3390/metabo13101101 | Benchmarking study of 66 GC-MS spectral similarity metrics and score families to improve identification reproducibility; includes CoreMS in the referenced software ecosystem. |
| 3 | Flores, J.E., Bramer, L.M., Degnan, D.J., Paurus, V.L., Corilo, Y.E., Clendinen, C.S. | Gaussian Mixture Modeling Extensions for Improved False Discovery Rate Estimation in GC-MS Metabolomics | 2023 | Journal of the American Society for Mass Spectrometry, Vol. 34, No. 6, pp. 1096-1104 | https://doi.org/10.1021/jasms.3c00039 | Methodological extension for improved FDR estimation in GC-MS identification workflows; references CoreMS in the computational pipeline context. |
| 4 | Letourneau, D.R., August, D.D., Volmer, D.A. | New algorithms demonstrate untargeted detection of chemically meaningful changing units and formula assignment for HRMS data of polymeric mixtures in the open-source Constellation web application | 2023 | Journal of Cheminformatics, Vol. 15, Article 7 | https://doi.org/10.1186/s13321-023-00680-5 | Introduces unsupervised HRMS algorithms and integrates a molecular formula assignment routine from CoreMS into Constellation for polymeric trend detection and formula annotation. |


---

## Application & Use Case Papers (13)

Papers that apply CoreMS to specific scientific domains, analytical workflows, or real-world problems.

| No. | Authors | Title | Year | Journal/Conference | DOI/Link | Application Domain | Key Results |
|-----|---------|-------|------|-------------------|----------|-------------------|-------------|
| 1 | Ferrer-Gonzalez, F.X., Heal, K.R., Sacks, J.S., Romero-Maysonet, Y., Finch, A.H., Carlson, L.T., Coe, L.S.Y., Bartolek, Z., Luthy, C.B., Gaffney, M.P., Angier, S.K., Flynn, S.E., Bachmann Gomez, C.I., Dunn, J.M.L., Bay, K.J., Yamamoto, L.A., Tien, M.Z., Armbrust, E.V., Durham, B.P., Sosa, O.A., Ingalls, A.E. | Conserved pathway for homarine catabolism in environmental bacteria | 2026 | Nature Microbiology | https://doi.org/10.1038/s41564-026-02313-7 | Marine microbial metabolism and homarine catabolism | Identified the conserved homABCDER operon for homarine degradation across diverse bacterial clades, with N-methylglutamic acid and glutamic acid as key products. CoreMS was used for untargeted feature detection via a persistent homology approach and isotopologue pairing (2H3, 2H2, 13C) to support methyl-transfer feature discovery in amendment experiments. |
| 2 | Timilsina, A., Lokesh, S., Shahriar, A., Basyal, S., Stincone, P., Schramm, T., Berndt, A., Dewey, C., Boiteau, R.M., Petras, D., Yang, Y. | Unveiling Redox-Active Quinones in Pyrogenic Carbon by Nontargeted Metabolomics and Dual Chemical Tagging | 2026 | Environmental Science & Technology | https://doi.org/10.1021/acs.est.6c01199 | Redox-active organic matter in pyrogenic carbon | Nontargeted metabolomics and dual-tagging workflow to characterize quinones in pyrogenic carbon; cites CoreMS in the HRMS processing workflow. |
| 3 | Kim, J., Wilson, S.J., Myers-Pigg, A., Bociu, I., Kemner, K.M., Regier, P., Rich, R., Kew, W., Megonigal, J.P., Bailey, V., Ward, N. | Variations in the Optical and Molecular Composition of Dissolved Organic Matter Exported from Coastal Wetlands | 2026 | Environmental Science & Technology (ASAP) | https://doi.org/10.1021/acs.est.5c11869 | Dissolved Organic Matter (DOM) analysis in coastal wetlands | LC-Orbitrap MS data processing including spectrum averaging, mass recalibration, and formula annotation |
| 4 | Shi, J., Tanentzap, A.J., Sun, Y., Wang, J., Xing, B., Rillig, M.C., Li, C., Jin, L., Wang, F., Adyel, T.M., Shang, J., Wang, X., Wang, J. | Microplastics Generate Less Mineral Protection of Soil Carbon and More CO2 Emissions | 2025 | Advanced Science, Vol. 12, No. 7, 2409585 | https://doi.org/10.1002/advs.202409585 | Soil carbon cycling and microplastic-derived DOM | Microcosm experiments showed MP-DOM increased CO2 emissions and reduced mineral-associated organic carbon relative to NOM; references CoreMS for molecular-level DOM characterization. |
| 5 | Nyarko, L., Dewey, C., Nason, J.A., Boiteau, R.M. | Tracking Changes in Organic-Copper Speciation during Wastewater Treatment Using LC-ICPMS-ESIMS | 2025 | ACS Environmental Au, Vol. 5, No. 2, pp. 230-240 | https://doi.org/10.1021/acsenvironau.4c00114 | Wastewater copper speciation and treatment chemistry | Demonstrates LC-ICPMS-ESIMS workflow for tracking organic-copper complexes through treatment; cites CoreMS in mass feature and formula processing workflow. |
| 6 | Shi, C., Mudunuru, M., Bowman, M., Zhao, Q., Toyoda, J., Kew, W., Corilo, Y., Qafoku, O., Bargar, J.R., Karra, S., Graham, E.B. | Scaling High-Resolution Soil Organic Matter Composition to Improve Predictions of Potential Soil Respiration Across the Continental United States | 2025 | Geophysical Research Letters, Vol. 52, No. 4, e2024GL113091 | https://doi.org/10.1029/2024GL113091 | Soil organic matter composition and carbon cycling | Machine learning approach to extract SOM formula features; demonstrated improved soil respiration predictions when combined with soil physicochemistry |
| 7 | Boiteau, R.M., Corilo, Y.E., Kew, W.R., Dewey, C., Alvarez Rodriguez, M.C., Carlson, C.A., Conway, T.M. | Relating Molecular Properties to the Persistence of Marine Dissolved Organic Matter with Liquid Chromatography–Ultrahigh-Resolution Mass Spectrometry | 2024 | Environmental Science & Technology, Vol. 58, No. 7, pp. 3267-3277 | https://doi.org/10.1021/acs.est.3c08245 | Marine dissolved organic matter (DOM) lability and persistence | Novel CoreMS data pipeline for molecular formula assignments and isomeric complexity metrics; identified molecular properties associated with DOM lability fractions in North Atlantic water column |
| 8 | Coffey, N.R., Dewey, C., Manning, K., Corilo, Y., Kew, W., Babcock-Adams, L., McKenna, A.M., Stuart, R.K., Boiteau, R.M. | Annotation of DOM metabolomes with an ultrahigh resolution mass spectrometry molecular formula library | 2024 | Organic Geochemistry, Vol. 197, p. 104880 | https://doi.org/10.1016/j.orggeochem.2024.104880 | Metabolomics and dissolved organic matter annotation | FT-ICR MS molecular formula library approach achieving 53% annotation rate for LC-MS metabolomic features; demonstrated on marine diatom exometabolome |
| 9 | Kelliher, J.M., Xu, Y., Flynn, M.C., Babinski, M., Canon, S., Cavanna, E., Clum, A., Corilo, Y.E., Fujimoto, G., Giberson, C., et al. | Standardized and accessible multi-omics bioinformatics workflows through the NMDC EDGE resource | 2024 | Computational and Structural Biotechnology Journal, Vol. 23, pp. 3575-3583 | https://doi.org/10.1016/j.csbj.2024.09.018 | Multi-omics workflow standardization and accessibility | Describes accessible standardized workflows across omics modalities, including integration of CoreMS within the NMDC workflow stack. |
| 10 | Sakas, J., Kitson, E., Bell, N.G.A., Uhrin, D. | MS and NMR Analysis of Isotopically Labeled Chloramination Disinfection Byproducts: Hyperlinks and Chemical Reactions | 2024 | Analytical Chemistry, Vol. 96, No. 21, pp. 8263-8272 | https://doi.org/10.1021/acs.analchem.3c03888 | Disinfection byproduct transformation analysis | Combined isotopic labeling, FT-ICR MS, and NMR for DBP characterization; references CoreMS in formula assignment and data processing workflow. |
| 11 | Kew, W., Myers-Pigg, A., Chang, C.H., Colby, S.M., Eder, J., Tfaily, M.M., Hawkes, J., Chu, R.K., Stegen, J.C. | Reviews and syntheses: Opportunities for robust use of peak intensities from high-resolution mass spectrometry in organic matter studies | 2024 | Biogeosciences, Vol. 21, No. 20, pp. 4665-4679 | https://doi.org/10.5194/bg-21-4665-2024 | Review and synthesis of best practices for utilizing peak intensities from HRMS in organic matter studies |
| 12 | Letourneau, D.R., Marzullo, B.P., Alexandridou, A., Barrow, M.P., O'Connor, P.B., Volmer, D.A. | Characterizing lignins from various sources and treatment processes after optimized sample preparation techniques and analysis via ESI-HRMS and custom mass defect software tools | 2023 | Analytical and Bioanalytical Chemistry, Vol. 415, No. 27, pp. 6663-6675 | https://doi.org/10.1007/s00216-023-04942-x | Lignin characterization and sample preparation optimization | Establishes optimized lignin sample-prep and ESI-HRMS conditions, with custom mass-defect analysis software; cites CoreMS among HRMS data-processing tools. |
| 13 | Bahureksa, W., Borch, T., Young, R.B., Weisbrod, C.R., Blakney, G.T., McKenna, A.M. | Improved Dynamic Range, Resolving Power, and Sensitivity Achievable with FT-ICR Mass Spectrometry at 21 T Reveals the Hidden Complexity of Natural Organic Matter | 2022 | Analytical Chemistry, Vol. 94, No. 32, pp. 11382-11389 | https://doi.org/10.1021/acs.analchem.2c02377 | Ultrahigh-field FT-ICR MS for natural organic matter | Demonstrates analytical gains at 21 T for complex NOM characterization and cites CoreMS in the downstream computational workflow. |

---

## Additional Citations & References (4)

Other references including dissertations, technical reports, conference abstracts, or works in progress that mention CoreMS.

| No. | Authors | Title | Year | Type | DOI/Link | Notes |
|-----|---------|-------|------|------|----------|-------|
| 1 | Eloe-Fadrosh, E.A., Ahmed, F., Babinski, M., et al. | The National Microbiome Data Collaborative Data Portal: an integrated multi-omics microbiome data resource | 2022 | Peer-reviewed article | https://doi.org/10.1093/nar/gkab990 | Data resource and portal paper that cites CoreMS in the broader NMDC software and workflow ecosystem. |
| 2 | Degnan, D.J., Claborne, D.M., White, A.M., Akers, S.M., Winans, N.M., Corilo, Y.E., Strauch, C.W., Bailey, V.L., McCue, L.A., Stratton, K.G., Bramer, L.M. | FREDA: A Web Application for the Processing, Analysis, and Visualization of Fourier-Transform Mass Spectrometry Data | 2025 | Peer-reviewed article | https://doi.org/10.1002/rcm.9980 | Compatible with CoreMS outputs and references CoreMS, but does not directly use CoreMS as part of the primary analysis workflow. |
| 3 | Ayala-Ortiz, C., Graf-Grachet, N., Freire-Zapata, V., et al. | MetaboDirect: an analytical pipeline for the processing of FT-ICR MS-based metabolomic data | 2023 | Peer-reviewed article | https://doi.org/10.1186/s40168-023-01476-3 | Cites CoreMS as a compatible data-processing tool, but the paper itself is not a CoreMS-use case or a CoreMS-focused methods paper. |
| 4 | Kitson, E., Kew, W., Ding, W., Bell, N.G.A. | PyKrev: A Python Library for the Analysis of Complex Mixture FT-MS Data | 2021 | Peer-reviewed article | https://doi.org/10.1021/jasms.1c00064 | Describes a related tool that can use CoreMS outputs and cites CoreMS, but does not directly use or evaluate CoreMS in the paper's primary workflow. |

---

## How to Add References

When adding a new publication:

1. Determine which category it belongs to (Methods, Application, or Additional)
2. Fill in the relevant fields in the table
3. Include a direct DOI link or publication URL when available
4. Add brief notes on key contributions or findings related to CoreMS
5. Keep entries ordered by year (most recent first) within each section

## Notes

- This document serves as a living bibliography for the CoreMS project
- Contributions are welcome! Please submit additions via pull requests
- Links should preferably be to peer-reviewed sources or official repositories
