# %% Import libs
import numpy as np
import pytest
import pandas as pd

from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectraCollection
from corems.mass_spectra.output.export import LCMSMetabolomicsExport, LCMSCollectionExport
from corems.encapsulation.factory.parameters import LCMSParameters, LCMSCollectionParameters
from corems.molecular_id.search.database_interfaces import MSPInterface


@pytest.fixture
def lcms_collection_folder(tmp_path, lcms_obj):
    """
    Creates a temporary folder with processed LCMS objects for collection testing.
    
    This fixture creates 3 samples with different levels of mass features:
    - Sample 1: Full set of mass features (all found features)
    - Sample 2: Partial set (first 50 mass features only)
    - Sample 3: No mass features (tests gap filling on completely empty sample)
    
    This setup allows comprehensive testing of gap filling functionality.
    """
    # Create a temporary folder for processed data
    processed_folder = tmp_path / "processed_lcms_collection"
    processed_folder.mkdir()
    
    # Set parameters on the LCMS object that are reasonable for testing
    lcms_obj.parameters = LCMSParameters(use_defaults=True)
    
    # Set persistent homology parameters for fast testing
    lcms_obj.parameters.lc_ms.peak_picking_method = "persistent homology"
    lcms_obj.parameters.lc_ms.ph_inten_min_rel = 0.001
    lcms_obj.parameters.lc_ms.ph_persis_min_rel = 0.05
    lcms_obj.parameters.lc_ms.ph_smooth_it = 0
    lcms_obj.parameters.lc_ms.ms1_scans_to_average = 3
    
    # MS1 parameters for quick testing
    ms1_params = lcms_obj.parameters.mass_spectrum['ms1']
    ms1_params.mass_spectrum.noise_threshold_method = "relative_abundance"
    ms1_params.mass_spectrum.noise_threshold_min_relative_abundance = 0.1
    ms1_params.mass_spectrum.noise_min_mz = 0
    ms1_params.mass_spectrum.min_picking_mz = 0
    ms1_params.mass_spectrum.noise_max_mz = np.inf
    ms1_params.mass_spectrum.max_picking_mz = np.inf
    
    # Process SAMPLE 1: Find and integrate mass features (FULL)
    lcms_obj.find_mass_features(assign_ms2_scans=True)
    lcms_obj.integrate_mass_features(drop_if_fail=True)
    
    # Save all mass features for later use
    all_mass_features = dict(lcms_obj.mass_features)
    all_eics = dict(lcms_obj.eics)
    
    # Export sample 1 with ALL mass features
    sample_name_1 = "test_sample_01"
    exporter1 = LCMSMetabolomicsExport(str(processed_folder / sample_name_1), lcms_obj)
    exporter1.to_hdf(overwrite=True)
    
    # Create SAMPLE 2: Partial set (first 50 mass features only)
    sample_name_2 = "test_sample_02"
    
    # Take only the first 50 mass features
    first_50_mf_ids = list(lcms_obj.mass_features.keys())[:50]
    first_mass_features = {mf_id: lcms_obj.mass_features[mf_id] for mf_id in first_50_mf_ids}
    lcms_obj.mass_features = first_mass_features
    
    # Export sample 2 with partial mass features
    exporter2 = LCMSMetabolomicsExport(str(processed_folder / sample_name_2), lcms_obj)
    exporter2.to_hdf(overwrite=True)
    
    # Create SAMPLE 3: No mass features at all (EMPTY)
    sample_name_3 = "test_sample_03"
    
    # Clear all mass features and EICs
    lcms_obj.mass_features = {}
    lcms_obj.eics = {}
    
    # Export sample 3 with no mass features
    exporter3 = LCMSMetabolomicsExport(str(processed_folder / sample_name_3), lcms_obj)
    exporter3.to_hdf(overwrite=True)
    
    # Create a manifest file to explicitly set sample 1 as the center
    # This is important because sample 1 has all features, making it the best reference
    import csv
    manifest_path = processed_folder / "manifest.csv"
    with open(manifest_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['sample_name', 'batch', 'order', 'center'])
        writer.writerow(['test_sample_01', 1, 1, True])   # Sample 1 is center (has all features)
        writer.writerow(['test_sample_02', 1, 2, False])  # Sample 2 (partial features)
        writer.writerow(['test_sample_03', 1, 3, False])  # Sample 3 (no features)
    
    return processed_folder


@pytest.fixture
def lcms_collection(lcms_collection_folder):
    """
    Creates an LCMSCollection object from processed LCMS data.
    
    Returns a collection with 3 samples for testing collection-level operations:
    - Sample 1: Full set of mass features
    - Sample 2: Partial set (first 50 mass features)
    - Sample 3: No mass features (for gap filling testing)
    """
    # Load the collection from the processed folder
    manifest_file = lcms_collection_folder / "manifest.csv"
    parser = ReadCoreMSHDFMassSpectraCollection(
        folder_location=lcms_collection_folder,
        manifest_file=manifest_file,
        cores=1
    )
    
    # Get the LCMS collection without light loading since sample 3 has no mass features
    # load_light=True would fail on sample 3 because it calls mass_features_to_df()
    collection = parser.get_lcms_collection(load_raw=False, load_light=False)
    
    # Adjust collection parameters to allow clusters with only 1 sample
    # This ensures features aren't filtered out before gap filling can occur
    collection.parameters.lcms_collection.cluster_size_min_samples = 1
    collection.parameters.lcms_collection.cluster_size_min_sample_percentage = 0.0
    
    # Set more lenient anchor feature parameters for alignment
    # Use relative_intensity with threshold 0 to accept all features as anchors
    collection.parameters.lcms_collection.mass_feature_anchor_technique = ['relative_intensity']
    collection.parameters.lcms_collection.mass_feature_anchor_relative_intensity_threshold = 0.0
    
    # Set reasonable alignment tolerances
    collection.parameters.lcms_collection.alignment_mz_tol_ppm = 5  # Tight m/z tolerance  
    collection.parameters.lcms_collection.alignment_rt_tol = 0.2  # 12 second RT tolerance
    
    # Set MS2 processing parameters for centroid data
    for lcms_obj in collection:
        ms2_params = lcms_obj.parameters.mass_spectrum['ms2']
        ms2_params.mass_spectrum.noise_threshold_method = "relative_abundance"
        ms2_params.mass_spectrum.noise_threshold_min_relative_abundance = 1.0
    
    return collection


def test_lcms_collection_creation(lcms_collection):
    """Test that an LCMSCollection can be created and has expected properties."""
    # Check that the collection was created
    assert lcms_collection is not None
    
    # Check number of samples (should be 3)
    assert len(lcms_collection) == 3
    assert len(lcms_collection.samples) == 3
    
    # Check that we can access individual LCMS objects
    assert lcms_collection[0] is not None
    assert lcms_collection[1] is not None
    assert lcms_collection[2] is not None
    
    # Check that manifest was created
    assert lcms_collection.manifest is not None
    assert len(lcms_collection.manifest) == 3
    
    # Check manifest dataframe
    manifest_df = lcms_collection.manifest_dataframe
    assert len(manifest_df) == 3
    assert 'batch' in manifest_df.columns
    assert 'order' in manifest_df.columns


def test_lcms_collection_mass_features_dataframe(lcms_collection):
    """Test that mass features from all samples are combined correctly."""
    # Get mass features dataframe
    mf_df = lcms_collection.mass_features_dataframe
    
    # Check that dataframe exists and has data
    assert mf_df is not None
    assert len(mf_df) > 0
    
    # Check that required columns exist
    required_columns = ['mf_id', 'sample_name', 'mz', 'scan_time', 'intensity']
    for col in required_columns:
        assert col in mf_df.columns, f"Missing required column: {col}"
    
    # Check that we have features from sample 1 (sample 2 has partial, sample 3 has none initially)
    unique_samples = mf_df['sample_name'].unique()
    assert len(unique_samples) >= 1
    
    # Check that coll_mf_id exists (collection-level unique ID) - it's the index
    assert mf_df.index.name == 'coll_mf_id' or 'coll_mf_id' in mf_df.columns


def test_lcms_collection_rt_alignment(lcms_collection):
    """Test retention time alignment across samples in the collection."""
    # Check initial state
    assert not lcms_collection.rt_aligned
    
    # Test plotting TICs before alignment
    try:
        import matplotlib
        matplotlib.use('Agg')  # Use non-interactive backend for testing
        lcms_collection.plot_tics(ms_level=1, type="raw", plot_legend=False)
    except Exception as e:
        pytest.fail(f"plot_tics raised an exception: {e}")
    
    # Perform alignment - this should succeed
    lcms_collection.align_lcms_objects()
    
    # Check that alignment was attempted (it may not use spline if already well-aligned)
    assert lcms_collection.rt_alignment_attempted
    
    # Check that alignment was rejected
    assert not lcms_collection.rt_aligned

    # Check that scan_time_aligned was added to all samples
    for i, lcms_obj in enumerate(lcms_collection):
        sample_name = lcms_collection.samples[i]
        assert 'scan_time_aligned' in lcms_obj.scan_df.columns, f"Missing scan_time_aligned in {sample_name}"
    
    # Test plotting after alignment - both raw and corrected TICs
    try:
        lcms_collection.plot_tics(ms_level=1, type="both", plot_legend=True)
    except Exception as e:
        pytest.fail(f"plot_tics with type='both' raised an exception: {e}")
    
    # Test plotting alignments (shows time differences)
    try:
        lcms_collection.plot_alignments(plot_legend=True)
    except Exception as e:
        pytest.fail(f"plot_alignments raised an exception: {e}")


def test_lcms_collection_consensus_features(lcms_collection):
    """Test generation of consensus mass features (clustering)."""
    # Ensure alignment is done first
    if not lcms_collection.rt_aligned:
        lcms_collection.align_lcms_objects()
    
    # Generate consensus features
    lcms_collection.add_consensus_mass_features()
    
    # Check cluster summary dataframe
    cluster_summary = lcms_collection.cluster_summary_dataframe
    assert cluster_summary is not None
    assert len(cluster_summary) > 0
    
    # Check that cluster column was added to mass features
    mf_df = lcms_collection.mass_features_dataframe
    assert 'cluster' in mf_df.columns
    
    # Check cluster feature dictionary
    cluster_dict = lcms_collection.cluster_feature_dictionary
    assert cluster_dict is not None
    assert len(cluster_dict) > 0
    
    # Verify that all clusters in summary are in the dictionary
    for cluster_id in cluster_summary.index:
        assert cluster_id in cluster_dict


def test_lcms_collection_gap_filling(lcms_collection):
    """Test gap filling to create induced mass features."""
    # Setup: align and cluster first
    if not lcms_collection.rt_aligned:
        lcms_collection.align_lcms_objects()
    lcms_collection.add_consensus_mass_features()
    
    # Verify initial state of each sample before gap filling
    sample_1_mf_count = len(lcms_collection[0].mass_features)
    sample_2_mf_count = len(lcms_collection[1].mass_features)
    sample_3_mf_count = len(lcms_collection[2].mass_features)
    
    # Perform gap filling
    pipeline_results = lcms_collection.process_consensus_features(
        load_representatives=False,
        perform_gap_filling=True,
        add_ms1=False,  
        add_ms2=False,
        molecular_formula_search=False,
        ms2_spectral_search=False,
        spectral_lib=False,
        molecular_metadata=None,
        gather_eics=True,
        keep_raw_data=False
    )
    
    # Check that induced mass features dataframe exists
    induced_df = lcms_collection.induced_mass_features_dataframe
    assert induced_df is not None
    
    # With sample 3 having missing features, gap filling should create induced features
    assert len(induced_df) > 0, "Gap filling should create induced mass features in sample 3"
    
    # Check that induced mass features have proper columns
    assert 'cluster' in induced_df.columns
    assert 'sample_name' in induced_df.columns
    assert 'mf_id' in induced_df.columns
    
    # Check induced features per sample in the dataframe (not individual objects)
    sample_3_induced = len(induced_df[induced_df['sample_id'] == 2])
    
    # Sample 3 should have induced features (started with 0, all 50 clusters are missing)
    assert sample_3_induced == 50, "Sample 3 should have 50 induced mass features (one for each cluster)"
    
    # By design, individual sample objects should have empty induced_mass_features dict
    # because they are collected into the induced_mass_features_dataframe
    assert len(lcms_collection[0].induced_mass_features) == 0
    assert len(lcms_collection[1].induced_mass_features) == 0
    assert len(lcms_collection[2].induced_mass_features) == 0


def test_lcms_collection_pivot_table(lcms_collection):
    """Test creation of pivot tables for collection data."""
    # Setup: ensure we have clustered features
    if not lcms_collection.rt_aligned:
        lcms_collection.align_lcms_objects()
    lcms_collection.add_consensus_mass_features()
    
    # Create pivot table with default attribute (coll_mf_id) before gap-filling
    pivot_df = lcms_collection.collection_pivot_table(verbose=False)
    
    # Check that pivot table was created
    assert pivot_df is not None
    assert isinstance(pivot_df, pd.DataFrame)
    
    # Check that all samples are included (even sample 3 with no features)
    assert len(pivot_df.columns) == len(lcms_collection.samples)
    assert 'test_sample_03' in pivot_df.columns
    
    # Sample 3 should have all NAs before gap-filling
    assert pivot_df['test_sample_03'].isna().all(), "Sample 3 should have all NAs before gap-filling"
    
    # Perform gap-filling
    lcms_collection.process_consensus_features(
        load_representatives=False,
        perform_gap_filling=True,
        add_ms1=False,
        add_ms2=False,
        molecular_formula_search=False,
        ms2_spectral_search=False,
        spectral_lib=False,
        molecular_metadata=None,
        gather_eics=True,
        keep_raw_data=False
    )
    
    # Create pivot table again after gap-filling
    pivot_df_after = lcms_collection.collection_pivot_table(verbose=False)
    
    # Sample 3 should no longer have all NAs after gap-filling
    assert not pivot_df_after['test_sample_03'].isna().all(), "Sample 3 should have filled features after gap-filling"
    
    # Check that sample 3 has exactly 50 non-NA values (one for each cluster)
    assert pivot_df_after['test_sample_03'].notna().sum() == 50, "Sample 3 should have 50 gap-filled features"
    
    # Create pivot table with intensity attribute
    pivot_intensity = lcms_collection.collection_pivot_table(
        attribute='intensity', 
        verbose=False
    )
    assert pivot_intensity is not None
    assert len(pivot_intensity.columns) == len(lcms_collection.samples)


def test_lcms_collection_cluster_representatives(lcms_collection):
    """Test extraction of representative features for each cluster."""
    # Setup: ensure we have clustered features
    if not lcms_collection.rt_aligned:
        lcms_collection.align_lcms_objects()
    lcms_collection.add_consensus_mass_features()
    
    # Get cluster representatives table
    reps_table = lcms_collection.cluster_representatives_table()
    
    # Check that table was created
    assert reps_table is not None
    assert isinstance(reps_table, pd.DataFrame)
    assert len(reps_table) > 0
    
    # Check required columns
    required_cols = ['cluster', 'coll_mf_id', 'mz', 'scan_time', 'intensity']
    for col in required_cols:
        assert col in reps_table.columns, f"Missing column: {col}"
    
    # Check that each cluster has exactly one representative
    cluster_counts = reps_table['cluster'].value_counts()
    assert all(count == 1 for count in cluster_counts)


def test_lcms_collection_export_import_hdf5(lcms_collection, tmp_path):
    """Test exporting and re-importing a collection from HDF5."""
    # Setup: align and cluster
    if not lcms_collection.rt_aligned:
        lcms_collection.align_lcms_objects()
    lcms_collection.add_consensus_mass_features()
    
    # Export collection to HDF5
    export_path = tmp_path / "test_collection"
    exporter = LCMSCollectionExport(
        out_file_path=str(export_path),
        mass_spectra_collection=lcms_collection
    )
    exporter.export_to_hdf5(overwrite=True, save_parameters=True)
    
    # Check that HDF5 file was created
    hdf5_path = export_path.with_suffix('.hdf5')
    assert hdf5_path.exists()
    
    # Re-import the collection
    from corems.mass_spectra.input.corems_hdf5 import ReadSavedLCMSCollection
    reader = ReadSavedLCMSCollection(
        collection_hdf5_path=str(hdf5_path),
        cores=1
    )
    collection2 = reader.get_lcms_collection(
        load_raw=False,
        load_light=True
    )
    
    # Verify the reloaded collection
    assert len(collection2) == len(lcms_collection)
    
    # Check that mass features match
    mf_count_1 = len(lcms_collection.mass_features_dataframe)
    mf_count_2 = len(collection2.mass_features_dataframe)
    assert mf_count_1 == mf_count_2
    
    # Check that cluster count matches if clustering was done
    if len(lcms_collection.cluster_summary_dataframe) > 0:
        cluster_count_1 = len(lcms_collection.cluster_summary_dataframe)
        cluster_count_2 = len(collection2.cluster_summary_dataframe)
        assert cluster_count_1 == cluster_count_2


def test_lcms_collection_drop_isotopologues(lcms_collection):
    """Test dropping isotopologues from the collection."""
    # Get initial mass features count
    initial_mf_count = len(lcms_collection.mass_features_dataframe)
    
    # Drop isotopologues
    lcms_collection._drop_isotopologues()
    
    # Check that flag was set
    assert lcms_collection.isotopes_dropped
    
    # Mass features count should be less than or equal to initial
    # (equal if no isotopologues were found)
    final_mf_count = len(lcms_collection.mass_features_dataframe)
    assert final_mf_count <= initial_mf_count


def test_lcms_collection_feature_annotations_table(lcms_collection, msp_file_location):
    """Test creation of feature annotations table with molecular metadata."""
    # Setup: align and cluster
    if not lcms_collection.rt_aligned:
        lcms_collection.align_lcms_objects()
    lcms_collection.add_consensus_mass_features()
    
    # Load molecular metadata from MSP file with same parameters as test_lcms_metabolomics
    my_msp = MSPInterface(file_path=msp_file_location)
    msp_lib, molecular_metadata = my_msp.get_metabolomics_spectra_library(
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
    
    # Set MS2 score threshold to match test_lcms_metabolomics
    for lcms_obj in lcms_collection:
        lcms_obj.parameters.lc_ms.ms2_min_fe_score = 0.3
    
    # Process consensus features with MS1, MS2, and spectral search
    # This should add molecular annotations before creating the annotations table
    pipeline_results = lcms_collection.process_consensus_features(
        load_representatives=True,
        perform_gap_filling=False,
        add_ms1=True,
        add_ms2=True,
        molecular_formula_search=False,
        ms2_spectral_search=True,
        spectral_lib=msp_lib,
        molecular_metadata=molecular_metadata,
        gather_eics=False,
        keep_raw_data=False
    )
    
    # Create annotations table with metadata
    annotations_table = lcms_collection.feature_annotations_table(
        molecular_metadata=molecular_metadata,
        drop_unannotated=False
    )
    
    assert annotations_table is not None
    assert isinstance(annotations_table, pd.DataFrame)
    assert len(annotations_table) > 0
    
    # Check that cluster information is present
    assert 'cluster' in annotations_table.columns
    
    # Check that we got some spectral matches after processing
    # Look for Entropy Similarity column which indicates MS2 spectral search results
    if 'Entropy Similarity' in annotations_table.columns:
        matched_features = annotations_table[annotations_table['Entropy Similarity'].notna()]
        assert len(matched_features) > 0, "Should have at least some MS2 spectral matches after search"
    else:
        # If column doesn't exist, the test should fail
        raise AssertionError("Expected 'Entropy Similarity' column in annotations table after MS2 spectral search")


def test_lcms_collection_molecular_formula_search(lcms_collection, postgres_database):
    """Test molecular formula search on consensus features."""
    # Setup: align and cluster
    if not lcms_collection.rt_aligned:
        lcms_collection.align_lcms_objects()
    lcms_collection.add_consensus_mass_features()
    
    # Set molecular search parameters for all samples (negative mode)
    for lcms_obj in lcms_collection:
        ms1_params = lcms_obj.parameters.mass_spectrum['ms1']
        ms1_params.molecular_search.url_database = postgres_database
        ms1_params.molecular_search.error_method = "None"
        ms1_params.molecular_search.min_ppm_error = -5
        ms1_params.molecular_search.max_ppm_error = 5
        ms1_params.molecular_search.mz_error_range = 1
        ms1_params.molecular_search.isProtonated = True  # Deprotonated in negative mode
        ms1_params.molecular_search.isRadical = False
        ms1_params.molecular_search.isAdduct = False
        ms1_params.molecular_search.usedAtoms = {
            'C': (1, 90),
            'H': (4, 200),
            'O': (0, 30),
            'N': (0, 3),
            'P': (0, 2),
            'S': (0, 2),
        }
    
    # Process consensus features with molecular formula search
    pipeline_results = lcms_collection.process_consensus_features(
        load_representatives=True,
        perform_gap_filling=False,
        add_ms1=True,
        add_ms2=False,
        molecular_formula_search=True,
        ms2_spectral_search=False,
        spectral_lib=False,
        molecular_metadata=None,
        gather_eics=False,
        keep_raw_data=False
    )
    
    # Get annotations table
    annotations_table = lcms_collection.feature_annotations_table(
        molecular_metadata=None,
        drop_unannotated=False
    )
    
    assert annotations_table is not None
    assert isinstance(annotations_table, pd.DataFrame)
    assert len(annotations_table) > 0
    
    # Check that molecular formula columns are present (column is called 'Ion Formula')
    assert 'Ion Formula' in annotations_table.columns
    assert 'Calculated m/z' in annotations_table.columns
    
    # Check that at least some features got molecular formula assignments
    assigned_formulas = annotations_table[annotations_table['Ion Formula'].notna()]
    assert len(assigned_formulas) > 0, "Should have at least some molecular formula assignments"
    
    # Verify the formulas are reasonable (contain expected elements)
    first_formula = assigned_formulas['Ion Formula'].iloc[0]
    assert 'C' in first_formula or 'H' in first_formula, "Molecular formulas should contain C or H"
    
    # Check that m/z error columns exist (indicates matching happened)
    assert 'm/z Error (ppm)' in annotations_table.columns
    assert 'm/z Error Score' in annotations_table.columns


def test_lcms_collection_update_raw_file_locations(lcms_collection, tmp_path):
    """Test updating raw file locations in the collection."""
    import shutil
    
    # Create a new path for raw files
    new_raw_folder = tmp_path / "new_raw_location"
    new_raw_folder.mkdir()
    
    # Copy the original raw file to the new location for each sample
    # The samples all use the same original raw file but have different sample names
    for lcms_obj in lcms_collection:
        original_raw = lcms_obj.raw_file_location
        # Create a copy with the sample name in the new location
        new_raw_file = new_raw_folder / f"{lcms_obj.sample_name}.raw"
        shutil.copy2(original_raw, new_raw_file)
    
    # Update raw file locations
    lcms_collection.update_raw_file_locations(str(new_raw_folder))
    
    # Check that paths were updated
    for lcms_obj in lcms_collection:
        assert str(new_raw_folder) in str(lcms_obj.raw_file_location)
    
    # Check that flag was set
    assert lcms_collection.raw_files_relocated


def test_lcms_collection_minimal_workflow(lcms_collection):
    """
    Test a minimal end-to-end workflow with the collection.
    
    This test mirrors the workflow in metabolomics_collection.py:
    1. Load collection
    2. Align retention times
    3. Generate consensus features
    4. Perform gap filling
    5. Create reports
    """
    # Step 1: Collection is loaded via fixture
    assert len(lcms_collection) > 0
    
    # Step 2: Align retention times
    lcms_collection.align_lcms_objects()
    
    # Step 3: Generate consensus features
    lcms_collection.add_consensus_mass_features()
    cluster_count = len(lcms_collection.cluster_summary_dataframe)
    assert cluster_count > 0
    
    # Step 4: Perform gap filling using process_consensus_features
    # Load representatives and add MS1 so we can create annotations table
    lcms_collection.process_consensus_features(
        load_representatives=True,
        perform_gap_filling=True,
        add_ms1=True,
        add_ms2=False,
        molecular_formula_search=False,
        ms2_spectral_search=False,
        spectral_lib=False,
        molecular_metadata=None,
        gather_eics=True,
        keep_raw_data=False
    )
    
    # Step 5: Create reports
    pivot_table = lcms_collection.collection_pivot_table(verbose=False)
    assert pivot_table is not None
    
    cluster_reps = lcms_collection.cluster_representatives_table()
    assert len(cluster_reps) == cluster_count
    
    # Create annotations table (requires load_representatives=True and add_ms1=True)
    annotations = lcms_collection.feature_annotations_table(
        molecular_metadata=None,
        drop_unannotated=False
    )
    assert annotations is not None
    assert len(annotations) > 0
    
    # Verify we have cluster information in the annotations
    assert 'cluster' in annotations.columns
    
    # Verify workflow completed successfully
    assert lcms_collection.rt_aligned or lcms_collection.rt_alignment_attempted
    assert cluster_count > 0
    print(f"\nWorkflow completed: {cluster_count} consensus clusters from {len(lcms_collection)} samples")

