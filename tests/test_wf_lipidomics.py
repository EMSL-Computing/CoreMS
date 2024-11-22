# %% Import libs
import sys

from pathlib import Path

import shutil
import numpy as np
import json

from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectra
from corems.mass_spectra.input.rawFileReader import ImportMassSpectraThermoMSFileReader
from corems.mass_spectra.input.mzml import MZMLSpectraParser
from corems.mass_spectra.output.export import LipidomicsExport
from corems.molecular_id.search.database_interfaces import MetabRefLCInterface
from corems.molecular_id.factory.lipid_molecular_metadata import LipidMetadata
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import LCMSParameters, reset_lcms_parameters, reset_ms_parameters


def test_import_lcmsobj_mzml():
    # Delete the "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.corems" directory
    # Instantiate parser based on binary file type
    file_mzml = (
        Path.cwd()
        / "tests/tests_data/lcms/"
        / "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.mzML"
    )

    parser = MZMLSpectraParser(file_mzml)

    # Instatiate lc-ms data object using parser and pull in ms1 spectra into dataframe (without storing as MassSpectrum objects to save memory)
    myLCMSobj = parser.get_lcms_obj(spectra="ms1")
    myLCMSobj.parameters = LCMSParameters(use_defaults=True)

    myLCMSobj.find_mass_features()
    myLCMSobj.add_associated_ms1(
        auto_process=True, use_parser=False, spectrum_mode="profile"
    )
    myLCMSobj.integrate_mass_features()
    mass_features_df = myLCMSobj.mass_features_to_df()
    assert mass_features_df.shape == (197, 12)
    
    # Reset the MSParameters to the original values
    reset_lcms_parameters()
    reset_ms_parameters()


def test_lipidomics_workflow():
    # Delete the "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.corems" directory
    shutil.rmtree(
        "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.corems",
        ignore_errors=True,
    )

    # Instantiate parser based on binary file type
    file_raw = (
        Path.cwd()
        / "tests/tests_data/lcms/"
        / "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.raw"
    )
    parser = ImportMassSpectraThermoMSFileReader(file_raw)

    # Instatiate lc-ms data object using parser and pull in ms1 spectra into dataframe (without storing as MassSpectrum objects to save memory)
    myLCMSobj = parser.get_lcms_obj(spectra="ms1")

    # Set parmaeters to the defaults for reproducible testing
    myLCMSobj.parameters = LCMSParameters(use_defaults=True)

    # Set parameters on the LCMS object that are reasonable for testing
    ## persistent homology parameters
    myLCMSobj.parameters.lc_ms.peak_picking_method = "persistent homology"
    myLCMSobj.parameters.lc_ms.ph_inten_min_rel = 0.0005
    myLCMSobj.parameters.lc_ms.ph_persis_min_rel = 0.05
    myLCMSobj.parameters.lc_ms.ph_smooth_it = 0
    myLCMSobj.parameters.lc_ms.ms2_min_fe_score = 0.3
    myLCMSobj.parameters.lc_ms.ms1_scans_to_average = 5

    ## MSParameters for ms1 mass spectra
    ms1_params = myLCMSobj.parameters.mass_spectrum['ms1']
    ms1_params.mass_spectrum.noise_threshold_method = "relative_abundance"
    ms1_params.mass_spectrum.noise_threshold_min_relative_abundance = 0.1
    ms1_params.mass_spectrum.noise_min_mz, ms1_params.mass_spectrum.min_picking_mz = 0, 0
    ms1_params.mass_spectrum.noise_max_mz, ms1_params.mass_spectrum.max_picking_mz = np.inf, np.inf
    ms1_params.ms_peak.legacy_resolving_power = False
    ms1_params.molecular_search.url_database = "postgresql://coremsdb:coremsmolform@postgres:5432/molformula"
    ms1_params.molecular_search.usedAtoms = {
        'C': (10, 30),
        'H': (18, 200),
        'O': (1, 23),
        'N': (0, 3),
        'P': (0, 1),
        'S': (0, 1),
    }

    ## settings for ms2 data (HCD scans)
    ms2_params_hcd = ms1_params.copy()
    ms2_params_hcd.molecular_search.ion_types_excluded = ["[M+HCOO]-"]
    myLCMSobj.parameters.mass_spectrum['ms2'] = ms2_params_hcd

    ## settings for ms2 data (CID scans)
    ms2_params_cid = ms2_params_hcd.copy()
    ms2_params_cid.molecular_search.max_ppm_error = 200 # wider ppm error for CID scans
    ms2_params_cid.mass_spectrum.noise_threshold_min_relative_abundance = 0.01 # lower noise threshold for CID scans
    myLCMSobj.parameters.mass_spectrum['ms2_cid'] = ms2_params_cid

    ## reporting settings
    myLCMSobj.parameters.lc_ms.search_as_lipids = True
    myLCMSobj.parameters.lc_ms.include_fragment_types = True
    myLCMSobj.parameters.lc_ms.export_eics = True
    myLCMSobj.parameters.lc_ms.export_profile_spectra = True

    # Use persistent homology to find mass features in the lc-ms data
    # Find mass features, cluster, and integrate them.  Then annotate pairs of mass features that are c13 iso pairs.

    myLCMSobj.find_mass_features()
    myLCMSobj.add_associated_ms1(
        auto_process=True, use_parser=False, spectrum_mode="profile"
    )

    myLCMSobj.integrate_mass_features(drop_if_fail=True)
    myLCMSobj.deconvolute_ms1_mass_features()

    mass_spec_decon = myLCMSobj.mass_features[1].mass_spectrum_deconvoluted
    assert len(mass_spec_decon.mspeaks) < len(
        myLCMSobj.mass_features[1].mass_spectrum.mspeaks
    )
    myLCMSobj.add_peak_metrics()
    myLCMSobj.find_c13_mass_features()
    assert len(myLCMSobj.mass_features) == 130

    # Perform a molecular search on a few of the mass features
    mf_df = myLCMSobj.mass_features_to_df() 
    unique_scans = mf_df.apex_scan.unique()
    i = 0
    for scan in unique_scans:
        if i > 1:  # only search first 3 scans for testing
            break
        print("searching mz for scan: ", str(i), " of ", str(len(unique_scans)))
        # gather mass features for this scan
        mf_df_scan = mf_df[mf_df.apex_scan == scan]
        peaks_to_search = [
            myLCMSobj.mass_features[x].ms1_peak for x in mf_df_scan.index.tolist()
        ]
        SearchMolecularFormulas(
            myLCMSobj._ms[scan],
            first_hit=False,
            find_isotopologues=True,
        ).run_worker_ms_peaks(peaks_to_search)
        i += 1

    # Check results of molecular search
    assert myLCMSobj.mass_features[0].ms1_peak[0].string == "C20 H30 O2"
    assert myLCMSobj.mass_features_ms1_annot_to_df().shape[0] == 130
    myLCMSobj.mass_features[0].mass_spectrum.to_dataframe()

    # Add hcd ms2 data to lcms object, using the ms2 mass spectrum parameters
    og_ms_len = len(myLCMSobj._ms)
    myLCMSobj.add_associated_ms2_dda(spectrum_mode="centroid", scan_filter="hcd")
    assert len(myLCMSobj._ms) > og_ms_len

    # Add cid ms2 data to lcms object, using the ms2_cid mass spectrum parameters
    og_ms_len = len(myLCMSobj._ms)
    myLCMSobj.add_associated_ms2_dda(spectrum_mode="centroid", ms_params_key="ms2_cid", scan_filter="cid")
    assert len(myLCMSobj._ms) > og_ms_len

    # Export the mass features to a pandas dataframe
    df = myLCMSobj.mass_features_to_df()
    assert df.shape == (130, 16)

    # Plot a mass feature
    myLCMSobj.mass_features[0].plot(return_fig=False)

    """
    # This code should be left as an example for how to generate example json data
    import dataclasses

    mzs = [i.mz for k, i in myLCMSobj.mass_features.items()]
    metabref = MetabRefLCInterface()
    metabref.set_token("tmp_data/thermo_raw_NMDC/metabref.token")
    spectra_library, lipid_metadata = metabref.get_lipid_library(
        mz_list=mzs[1:10],
        polarity="negative",
        mz_tol_ppm=5,
        mz_tol_da_api=0.01,
        format="json",
        normalize=True
    )
    # Save the json spectra library and lipid metadata to a text file and then load it back in
    import json
    with open('tests/tests_data/lcms/metabref_spec_lib.json', "w") as final:
        json.dump(spectra_library, final)
    lipid_metadata_raw = {
        k: dataclasses.asdict(v) for k, v in lipid_metadata.items()
        }
    with open('tests/tests_data/lcms/metabref_lipid_metadata.json', "w") as final:
        json.dump(lipid_metadata_raw, final)
    """
    metabref = MetabRefLCInterface()

    # Load an example json spectral library and convert to flashentropy format
    with open("tests/tests_data/lcms/metabref_spec_lib.json") as f:
        spectra_library_json = json.load(f)
    spectra_library_fe = metabref._to_flashentropy(
        spectra_library_json,
        normalize=True,
        fe_kwargs={
            "normalize_intensity": True,
            "min_ms2_difference_in_da": 0.02,  # for cleaning spectra
            "max_ms2_tolerance_in_da": 0.01,  # for setting search space
            "max_indexed_mz": 3000,
            "precursor_ions_removal_da": None,
            "noise_threshold": 0,
        },
    )
    
    # Load the associated lipid metadata and convert to correct class
    with open("tests/tests_data/lcms/metabref_lipid_metadata.json") as f:
        lipid_metadata_json = json.load(f)
    lipid_metadata = {
        k: metabref._dict_to_dataclass(v, LipidMetadata)
        for k, v in lipid_metadata_json.items()
    }

    # Perform a spectral search on the mass features
    hcd_ms2_scan_df = myLCMSobj.scan_df[
        myLCMSobj.scan_df.scan_text.str.contains("hcd")
        & (myLCMSobj.scan_df.ms_level == 2)
    ]
    ms2_scans_oi_hr = [
        x for x in hcd_ms2_scan_df.scan.tolist() if x in myLCMSobj._ms.keys()
    ]
    myLCMSobj.fe_search(
        scan_list=ms2_scans_oi_hr, fe_lib=spectra_library_fe, peak_sep_da=0.01
    )
    # Export the lcms object to an hdf5 file using the LipidomicsExport class
    exporter = LipidomicsExport(
        "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801", myLCMSobj
    )
    exporter.to_hdf(overwrite=True)
    exporter.report_to_csv(molecular_metadata=lipid_metadata)
    report = exporter.to_report(molecular_metadata=lipid_metadata)
    assert report['Ion Formula'][1] == 'C24 H47 O2'

    # Import the hdf5 file, assert that its df is same as above and that we can plot a mass feature
    parser = ReadCoreMSHDFMassSpectra(
        "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.corems/Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.hdf5"
    )
    myLCMSobj2 = parser.get_lcms_obj()

    # Check that the parameters match
    assert myLCMSobj2.parameters == myLCMSobj.parameters
    assert myLCMSobj2.spectra_parser_class.__name__ == "ImportMassSpectraThermoMSFileReader"
    df2 = myLCMSobj2.mass_features_to_df()
    assert df2.shape == (130, 16)
    myLCMSobj2.mass_features[0].mass_spectrum.to_dataframe()
    assert myLCMSobj2.mass_features[0].ms1_peak[0].string == "C20 H30 O2"
    assert myLCMSobj2.mass_features_ms1_annot_to_df().shape[0] == 130
    myLCMSobj2.mass_features[0].plot(return_fig=False)

    # Delete the "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.corems" directory
    shutil.rmtree(
        "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.corems",
        ignore_errors=True,
    )

    # Reset the MSParameters to the original values
    reset_lcms_parameters()
    reset_ms_parameters()
