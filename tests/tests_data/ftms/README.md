# tests/tests_data/ftms/

FT-MS/FTICR test fixtures, one subdirectory per dataset. Filenames are kept
as-is from their original acquisition/export — only the containing path
changed when this directory was reorganized (2026-09).

Provenance below is marked **confirmed** when it's established by direct
evidence (a matching filename fragment inside a `.d` folder's nested
metadata, an explicit reference in code/docs) and **inferred** when it's a
plausible guess from naming/context that hasn't been independently verified.
Where nothing links a file to a dataset, it's parked in `unreferenced/`
rather than deleted.

## Known coverage gaps

- No Bruker imaging-mode `.d` sample (i.e. one carrying `ImagingInfo.xml`
  instead of `scan.xml`). `corems/mass_spectra/input/brukerSolarix_utils.py`
  has a code path for it, but nothing in this repo exercises it.
- No 1ω/2ω-labeled acquisition. No code under `corems/transient/` or
  `corems/mass_spectra/input/brukerSolarix*.py` has dedicated 1ω/2ω
  handling today.

## Datasets

### `srfa_bruker_solarix_direct_infusion/`

Bruker solariX FTICR, direct-infusion (single-scan `fid`), Suwannee River
Fulvic Acid (SRFA), negative ESI.

- `ESI_NEG_SRFA.d/` — the acquisition folder (fid, apexAcquisition.method,
  ExciteSweep, etc).
- `SRFA.ref` — calibration reference mass list for this acquisition.
- `ESI_NEG_SRFA_UnCal_Unassign.csv` — uncalibrated/unassigned mass list
  exported from the same spectrum (confirmed by test usage: fixture below).
- `ESI_NEG_SRFA_COREMS_withdupes.csv` — CoreMS-format mass list export with
  duplicate entries, used to exercise dedup handling.

Consumers: the root `conftest.py` fixtures `ftms_file_location` →
`bruker_transient` → `mass_spectrum_ftms`, and `ref_file_location` — used
directly or indirectly by `test_mass_spectrum.py`, `test_mspeak.py`,
`test_calibration.py`, `test_classification.py`,
`test_molecular_formula_search.py`, `test_setting_settings.py`,
`test_output.py`, `test_search_mass_list.py`, and two tests in
`test_input.py`. Also `tests/archive_tests/test_assembly_identification.py`,
the root `README.md` quick-start, and several tutorial notebooks/scripts
under `examples/`.

### `srfa_bruker_solarix_autosampler/`

Bruker solariX FTICR, autosampler/serial-mode acquisition (`ser` +
`scan.xml`) of SRFA, negative ESI — same sample as above, different
acquisition mode (LC/autosampler vs. direct infusion).

- `NEG_ESI_SRFA_Auto.d/` — the acquisition folder.
- `NEG_ESI_SRFA_CoreMS.hdf5` — CoreMS HDF5 export of the processed spectrum.
- `NEG_ESI_SRFA_CoreMS.xlsx` — CoreMS Excel export of the processed
  spectrum.
- `NEG_ESI_SRFA_CoreMS.corems/` — CoreMS text/LCMS export directory.

Consumers: `test_input.py` (`test_import_lcms_from_transient`,
`test_import_corems_hdf5`, `test_import_mass_list`,
`test_import_corems_mass_list`), plus several tutorial notebooks.

### `srfa_maglab_pks/`

MagLab `.pks` peak-list export of SRFA (negative ESI).

- `SRFA.pks`

Consumer: `test_input.py::test_import_maglab_pks`.

### `srfa_thermo_orbitrap/`

Thermo Orbitrap `.raw` file of SRFA (negative ESI) — the only Orbitrap
sample in this directory (see also `tests/tests_data/lcms/` for other
Thermo instrument coverage).

- `SRFA_NEG_ESI_ORB.raw`

Consumer: `test_input.py::test_import_thermo_average`.

### `srfa_xml_masslist/`

Generic XML mass-list export of SRFA (negative ESI, centroid, ~36k peaks).

- `srfa_neg_xml_example.xml`

Consumer: `test_input.py::test_import_xml_mass_list`.

### `thermo_profile_masslist/`

Thermo-style profile-mode text mass list (not tied to SRFA — a generic
Thermo profile export used to test the profile-text parser path).

- `Thermo_Profile_MassList.txt`

Consumer: `test_input.py::test_import_thermo_profile_mass_list`.

### `esfa_booster_hdf5/`

Booster HDF5 export of Elliott Soil Fulvic Acid (ESFA), 100k resolution
setting, negative ESI.

- `ESFA_100k_9767-13548_chB.A_re_pc_CoAddAll_mFT.h5`

Consumers: `test_input.py::test_import_booster_mass_spectrum_hdf` and
`test_import_booster_mass_spectra_hdf`.

### `esfa_ascii_masslist/`

ASCII mass-list export of ESFA (negative ESI).

- `ESI_NEG_ESFA.ascii`

Consumers: `test_input.py::test_import_mass_list`; also referenced as an S3
key (`1/ESI_NEG_ESFA.ascii`) by `tests/archive_tests/s3_test.py`, which has
no local test functions and is unaffected by this reorg.

### `unreferenced/`

Files with no test, example, or notebook currently pointing at them.
Quarantined here rather than deleted so a future pass can decide to wire
them into a test, or remove them. Lineage below is **inferred**, not
verified against file contents.

- `NEG_ESI_LIGNIN.raw` — Thermo raw file of a lignin sample. No known
  consumer.
- `ESFA_15T_3sFID_calibrated.txt` — calibrated text export, presumably ESFA
  on a 15T instrument, 3s FID acquisition. Naming suggests a relationship
  to `esfa_booster_hdf5/`, but resolution/acquisition details don't
  obviously match, so the link is unconfirmed. No known consumer.
- `ESI_NEG_SRFA_COREMS.csv` — CoreMS CSV export of SRFA (note: distinct
  from `srfa_bruker_solarix_direct_infusion/ESI_NEG_SRFA_COREMS_withdupes.csv`,
  which nothing else in this repo currently loads without the file
  changing name). No known consumer.
- `Exploris_SRFA_Example.raw` — Thermo Exploris raw file of SRFA. Used only
  by `examples/archive/scripts/HR-MS Thermo Raw 21T.py` (an archived
  script, not part of the test suite).
- `Auto_SRFA_QC.csv`, `Auto_SRFA_QC II.csv` — QC exports; the nested
  acquisition folder name inside `srfa_bruker_solarix_autosampler/NEG_ESI_SRFA_Auto.d/`
  contains a matching `Auto_SRFA_QC` fragment, so these likely derive from
  that acquisition (inferred, not confirmed). Used only by
  `examples/archive/scripts/Molecular Formula Data Aggreation.py` (archived,
  not part of the test suite).
