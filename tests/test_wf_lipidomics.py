# %% Import libs
import sys

from pathlib import Path

import shutil
import numpy as np
import json

from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectra
from corems.mass_spectra.input.mzml import MZMLSpectraParser
from corems.mass_spectra.output.export import LipidomicsExport
from corems.molecular_id.search.database_interfaces import MetabRefLCInterface
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulasLC
from corems.encapsulation.factory.parameters import LCMSParameters, reset_lcms_parameters, reset_ms_parameters


def test_import_lcmsobj_mzml():
    # Instantiate parser based on binary file type
    file_mzml = (
        Path.cwd()
        / "tests/tests_data/lcms/"
        / "test_centroid_neg_RP_metab.mzML"
    )

    parser = MZMLSpectraParser(file_mzml)

    # Instatiate lc-ms data object using parser and pull in ms1 spectra into dataframe (without storing as MassSpectrum objects to save memory)
    myLCMSobj = parser.get_lcms_obj(spectra="ms1")
    myLCMSobj.parameters = LCMSParameters(use_defaults=True)

    # Modify parameters to deal with centroid data
    myLCMSobj.parameters.lc_ms.peak_picking_method = "centroided_persistent_homology"
    myLCMSobj.parameters.mass_spectrum[
            "ms1"
        ].mass_spectrum.noise_threshold_method = "relative_abundance"

    myLCMSobj.find_mass_features()
    myLCMSobj.add_associated_ms1(
        auto_process=True, use_parser=True, spectrum_mode="centroid"
    )
    myLCMSobj.integrate_mass_features()
    mass_features_df = myLCMSobj.mass_features_to_df()
    assert mass_features_df.shape == (1395, 12)
    
    # Reset the MSParameters to the original values
    reset_lcms_parameters()
    reset_ms_parameters()


def test_lipidomics_workflow(postgres_database, lcms_obj):
    # Delete the "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.corems" directory
    shutil.rmtree(
        "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.corems",
        ignore_errors=True,
    )

    # Set parmaeters to the defaults for reproducible testing
    lcms_obj.parameters = LCMSParameters(use_defaults=True)

    # Set parameters on the LCMS object that are reasonable for testing
    ## persistent homology parameters
    lcms_obj.parameters.lc_ms.peak_picking_method = "persistent homology"
    lcms_obj.parameters.lc_ms.ph_inten_min_rel = 0.0005
    lcms_obj.parameters.lc_ms.ph_persis_min_rel = 0.05
    lcms_obj.parameters.lc_ms.ph_smooth_it = 0
    lcms_obj.parameters.lc_ms.ms2_min_fe_score = 0.3
    lcms_obj.parameters.lc_ms.ms1_scans_to_average = 5

    ## MSParameters for ms1 mass spectra
    ms1_params = lcms_obj.parameters.mass_spectrum['ms1']
    ms1_params.mass_spectrum.noise_threshold_method = "relative_abundance"
    ms1_params.mass_spectrum.noise_threshold_min_relative_abundance = 0.1
    ms1_params.mass_spectrum.noise_min_mz, ms1_params.mass_spectrum.min_picking_mz = 0, 0
    ms1_params.mass_spectrum.noise_max_mz, ms1_params.mass_spectrum.max_picking_mz = np.inf, np.inf
    ms1_params.ms_peak.legacy_resolving_power = False
    ms1_params.molecular_search.url_database = postgres_database
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
    lcms_obj.parameters.mass_spectrum['ms2'] = ms2_params_hcd

    ## settings for ms2 data (CID scans)
    ms2_params_cid = ms2_params_hcd.copy()
    ms2_params_cid.molecular_search.max_ppm_error = 200 # wider ppm error for CID scans
    ms2_params_cid.mass_spectrum.noise_threshold_min_relative_abundance = 0.01 # lower noise threshold for CID scans
    lcms_obj.parameters.mass_spectrum['ms2_cid'] = ms2_params_cid

    ## reporting settings
    lcms_obj.parameters.lc_ms.search_as_lipids = True
    lcms_obj.parameters.lc_ms.include_fragment_types = True
    lcms_obj.parameters.lc_ms.export_eics = True
    lcms_obj.parameters.lc_ms.export_profile_spectra = True

    # Use persistent homology to find mass features in the lc-ms data
    # Find mass features, cluster, and integrate them.  Then annotate pairs of mass features that are c13 iso pairs.

    lcms_obj.find_mass_features()
    lcms_obj.add_associated_ms1(
        auto_process=True, use_parser=False, spectrum_mode="profile"
    )

    lcms_obj.integrate_mass_features(drop_if_fail=True)
    lcms_obj.deconvolute_ms1_mass_features()

    mass_spec_decon = lcms_obj.mass_features[1].mass_spectrum_deconvoluted
    assert len(mass_spec_decon.mspeaks) < len(
        lcms_obj.mass_features[1].mass_spectrum.mspeaks
    )
    lcms_obj.add_peak_metrics()
    lcms_obj.find_c13_mass_features()
    assert len(lcms_obj.mass_features) == 130

    # Perform a molecular search on all of the mass features' ms1 peaks
    mol_form_search = SearchMolecularFormulasLC(lcms_obj)
    mol_form_search.run_mass_feature_search()

    # Check results of molecular search
    assert lcms_obj.mass_features[0].ms1_peak[0].string == "C20 H30 O2"
    assert lcms_obj.mass_features_ms1_annot_to_df().shape[0] > 130
    lcms_obj.mass_features[0].mass_spectrum.to_dataframe()

    # Add hcd ms2 data to lcms object, using the ms2 mass spectrum parameters
    og_ms_len = len(lcms_obj._ms)
    lcms_obj.add_associated_ms2_dda(spectrum_mode="centroid", scan_filter="hcd")
    assert len(lcms_obj._ms) > og_ms_len

    # Add cid ms2 data to lcms object, using the ms2_cid mass spectrum parameters
    og_ms_len = len(lcms_obj._ms)
    lcms_obj.add_associated_ms2_dda(spectrum_mode="centroid", ms_params_key="ms2_cid", scan_filter="cid")
    assert len(lcms_obj._ms) > og_ms_len

    lcms_obj.plot_composite_mz_features()

    # Export the mass features to a pandas dataframe
    df = lcms_obj.mass_features_to_df()
    assert df.shape == (130, 16)

    # Plot a mass feature
    lcms_obj.mass_features[0].plot(return_fig=False)

    # Query the lipidomics database to prepare a small search library for the mass features
    metabref = MetabRefLCInterface()
    mzs = [i.mz for k, i in lcms_obj.mass_features.items()]
    spectra_library_fe, lipid_metadata = metabref.get_lipid_library(
            mz_list=mzs[1:10],
            polarity="negative",
            mz_tol_ppm=5,
            format="flashentropy",
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

    # Perform a spectral search on the mass features
    hcd_ms2_scan_df = lcms_obj.scan_df[
        lcms_obj.scan_df.scan_text.str.contains("hcd")
        & (lcms_obj.scan_df.ms_level == 2)
    ]
    ms2_scans_oi_hr = [
        x for x in hcd_ms2_scan_df.scan.tolist() if x in lcms_obj._ms.keys()
    ]
    lcms_obj.fe_search(
        scan_list=ms2_scans_oi_hr, fe_lib=spectra_library_fe, peak_sep_da=0.01
    )
    # Export the lcms object to an hdf5 file using the LipidomicsExport class
    exporter = LipidomicsExport(
        "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801", lcms_obj
    )
    exporter.to_hdf(overwrite=True)
    exporter.report_to_csv(molecular_metadata=lipid_metadata)
    report = exporter.to_report(molecular_metadata=lipid_metadata)
    assert report['Ion Formula'][1] == 'C24 H47 O2'
    #assert report['Lipid Species'][1] == 'FA 24:0'

    # Import the hdf5 file, assert that its df is same as above and that we can plot a mass feature
    parser = ReadCoreMSHDFMassSpectra(
        "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.corems/Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.hdf5"
    )
    myLCMSobj2 = parser.get_lcms_obj()

    # Check that the parameters match
    assert myLCMSobj2.parameters == lcms_obj.parameters
    assert myLCMSobj2.spectra_parser_class.__name__ == "ImportMassSpectraThermoMSFileReader"
    df2 = myLCMSobj2.mass_features_to_df()
    assert df2.shape == (130, 16)
    myLCMSobj2.mass_features[0].mass_spectrum.to_dataframe()
    assert myLCMSobj2.mass_features[0].ms1_peak[0].string == "C20 H30 O2"
    assert myLCMSobj2.mass_features_ms1_annot_to_df().shape[0] > 130
    myLCMSobj2.mass_features[0].plot(return_fig=False)

    # Delete the "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.corems" directory
    shutil.rmtree(
        "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.corems",
        ignore_errors=True,
    )

    # Reset the MSParameters to the original values
    reset_lcms_parameters()
    reset_ms_parameters()