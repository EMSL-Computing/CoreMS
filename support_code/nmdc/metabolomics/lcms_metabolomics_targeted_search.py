from pathlib import Path

from corems.mass_spectra.input.rawFileReader import ImportMassSpectraThermoMSFileReader

if __name__ == "__main__":
    # =============================================================================
    # Configuration
    # =============================================================================
    # Paths
    base_path = Path("/Volumes/LaCie/nmdc_data/collection_testing/dev_test/")
    raw_data_path = base_path / "raw"

    # Instantiate LCMS object from the first raw data file
    first_raw_file = next(raw_data_path.glob("*.raw"))
    parser = ImportMassSpectraThermoMSFileReader(first_raw_file)
    lcms_obj = parser.get_lcms_obj(spectra="ms1")
    assert lcms_obj is not None, "Failed to instantiate LCMS object."

    # Set parameters for the LCMS object, only ones we NEED to change from defaults to get things running
    # Ideally, we'd load carefully chosen parameters from a file here for general processing
    lcms_obj.parameters.mass_spectrum['ms2'].mass_spectrum.noise_threshold_method = "relative_abundance"

    # Prepare a search dictionary for a targeted search
    # Here are examples of target m/z and RT values
    # Could derive this from an msp file.
    target_mz_list = [254.2837, 521.325988, 503.3134460449219, 317.2839050292969]
    target_rt_list = [5.83, 6.017, 7.977, 6.8658]
    target_search_dict = {
        "target_mz_list": target_mz_list,
        "target_rt_list": target_rt_list,
        "mz_tolerance_ppm": 5, #QUESTION: per target or bulk?
        "rt_tolerance": 0.5, #QUESTION: per target or bulk?
        "type": "internal standard"
    }
    
    # Look for mass features in targeted search mode
    lcms_obj.find_mass_features(
        targeted_search = True, 
        target_search_dict=target_search_dict
        )
    # Alternatively, or additionally, we could do normal (untargeted) mass feature finding here
    # lcms_obj.find_mass_features(targeted_search=False)

    # Let's integrate, add ms1 and ms2
    lcms_obj.integrate_mass_features()
    lcms_obj.add_associated_ms1(use_parser=False, spectrum_mode="profile")
    lcms_obj.add_associated_ms2_dda(use_parser=True, spectrum_mode="centroid")

    # Here we could do formula searching and MS2 search against a MSP database with MS2s of the search targets (similar to test_lcms_metabolomics.py)
    # BUT, we'd need a compatible msp file, and that's not ultra simple
    
    # Check out the mass feature dataframe
    mf_df = lcms_obj.mass_features_to_df(drop_na_cols=True)

    # Can visualize as well
    lcms_obj.mass_features[0].plot(return_fig = False)

    # Here could do some post-processing to report the "best hit" per target and include confidence metrics
    print("here")