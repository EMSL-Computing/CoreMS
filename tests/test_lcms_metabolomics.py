# %% Import libs
import shutil
import numpy as np

from corems.mass_spectra.output.export import LCMSMetabolomicsExport
from corems.molecular_id.search.database_interfaces import MSPInterface
from corems.encapsulation.factory.parameters import LCMSParameters, reset_lcms_parameters, reset_ms_parameters
from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectra


def test_lcms_metabolomics(postgres_database, lcms_obj, msp_file_location):
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
        'C': (5, 30),
        'H': (18, 200),
        'O': (1, 23),
        'N': (0, 3),
        'P': (0, 1),
        'S': (0, 1),
    }

    ## settings for ms2 data (HCD scans)
    ms2_params_hcd = ms1_params.copy()
    lcms_obj.parameters.mass_spectrum['ms2'] = ms2_params_hcd

    ## reporting settings
    lcms_obj.parameters.lc_ms.export_eics = True
    lcms_obj.parameters.lc_ms.export_profile_spectra = True

    ## peak metrics filtering settings - test the new functionality
    lcms_obj.parameters.lc_ms.remove_mass_features_by_peak_metrics = True
    lcms_obj.parameters.lc_ms.mass_feature_attribute_filter_dict = {
        'dispersity_index': {'value': 0.5, 'operator': '<'}
    }

    # Use persistent homology to find mass features in the lc-ms data
    # Find mass features, cluster, and integrate them.  Then annotate pairs of mass features that are c13 iso pairs.
    lcms_obj.find_mass_features()
    lcms_obj.integrate_mass_features(drop_if_fail=True)
    
    # Record number of mass features before filtering for comparison
    num_features_before_filtering = len(lcms_obj.mass_features)
        # Add peak metrics and filter mass features based on the new parameters
    lcms_obj.add_peak_metrics()

    # Record number of mass features after filtering for comparison    
    num_features_after_filtering = len(lcms_obj.mass_features)
    assert num_features_after_filtering < num_features_before_filtering

    # Add associated ms1 after integration and peak shape metrics calculation and filtering
    # This prevents adding unnecessary ms1 data for mass features that will be filtered out
    lcms_obj.add_associated_ms1(
        auto_process=True, use_parser=False, spectrum_mode="profile"
    )

    # Add hcd ms2 data to lcms object, using the ms2 mass spectrum parameters
    og_ms_len = len(lcms_obj._ms)
    lcms_obj.add_associated_ms2_dda(spectrum_mode="centroid", scan_filter="hcd")
    assert len(lcms_obj._ms) > og_ms_len

    # Query the lipidomics database to prepare a small search library for the mass features
    my_msp = MSPInterface(file_path=msp_file_location)
    msp_negative, metabolite_metadata_negative = (
        my_msp.get_metabolomics_spectra_library(
            polarity="negative",
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
        scan_list=ms2_scans_oi_hr, fe_lib=msp_negative, peak_sep_da=0.01
    )

    # Export the lcms object to an hdf5 file using the LipidomicsExport class
    exporter = LCMSMetabolomicsExport(
        "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801_metab", lcms_obj
    )
    exporter.to_hdf(overwrite=True)
    exporter.report_to_csv(molecular_metadata=metabolite_metadata_negative)
    report = exporter.to_report(molecular_metadata=metabolite_metadata_negative)
    assert report['Ion Formula'][1] == 'C24 H47 O2'
    assert report['chebi'][1] == 28866

    # Test parameter re-import by loading the HDF5 file and comparing parameters
    parser = ReadCoreMSHDFMassSpectra(
        "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801_metab.corems/Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801_metab.hdf5"
    )
    lcms_obj_reimported = parser.get_lcms_obj()

    # Check that the parameters match, including the new peak metrics filtering parameters
    assert lcms_obj_reimported.parameters == lcms_obj.parameters, \
        "Re-imported parameters should match original parameters"
    
    # Specifically check the new peak metrics filtering parameters
    assert lcms_obj_reimported.parameters.lc_ms.remove_mass_features_by_peak_metrics
    assert lcms_obj_reimported.parameters.lc_ms.mass_feature_attribute_filter_dict == {
        'dispersity_index': {'value': 0.5, 'operator': '<'}
    }
    
    # Test that the number of mass features is preserved in the re-import
    assert len(lcms_obj_reimported.mass_features) == len(lcms_obj.mass_features), \
        "Re-imported LCMS object should have the same number of mass features as original"
    
    # Delete the "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.corems" directory
    shutil.rmtree(
        "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801_metab.corems",
        ignore_errors=True,
    )

    # Reset the MSParameters to the original values
    reset_lcms_parameters()
    reset_ms_parameters()