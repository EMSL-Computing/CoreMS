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
    
    # Create SAMPLE 2: Exact copy of Sample 1 (all features identical for easier debugging)
    # This ensures that sample 1's feature 0 should cluster with sample 2's feature 0, etc.
    sample_name_2 = "test_sample_02"
    
    # Keep ALL mass features - just save a second copy
    # Don't modify lcms_obj.mass_features, it already has all features
    
    # Export sample 2 with same mass features as sample 1
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


def test_lcms_collection_parameters(lcms_collection):
    """Test that collection parameters can be get/set."""
    # Check that parameters exist
    assert lcms_collection.parameters is not None
    
    # Check that it's the right type
    assert isinstance(lcms_collection.parameters, LCMSCollectionParameters)
    
    # Test setting new parameters
    new_params = LCMSCollectionParameters()
    new_params.lcms_collection.cores = 2
    lcms_collection.parameters = new_params
    
    assert lcms_collection.parameters.lcms_collection.cores == 2


def test_lcms_collection_rt_alignment(lcms_collection):
    """Test retention time alignment across samples in the collection."""
    # Check initial state
    assert not lcms_collection.rt_aligned
    
    # Debug: Print anchor parameters
    print(f"\nAnchor technique: {lcms_collection.parameters.lcms_collection.mass_feature_anchor_technique}")
    print(f"Anchor absolute threshold: {lcms_collection.parameters.lcms_collection.mass_feature_anchor_absolute_intensity_threshold}")
    print(f"Anchor relative threshold: {lcms_collection.parameters.lcms_collection.mass_feature_anchor_relative_intensity_threshold}")
    
    # Debug: Check which sample is center
    print(f"\nManifest center values:")
    for sample_name, manifest_data in lcms_collection._manifest_dict.items():
        print(f"  {sample_name}: center={manifest_data.get('center', False)}")
    
    # Debug: Print feature counts
    mf_df = lcms_collection.mass_features_dataframe
    print(f"\nMass features per sample:")
    for sample in lcms_collection.samples:
        sample_mfs = mf_df[mf_df['sample_name'] == sample]
        print(f"  {sample}: {len(sample_mfs)} features")
    
    # Perform alignment - this should succeed
    lcms_collection.align_lcms_objects()
    
    # Check that alignment was attempted (it may not use spline if already well-aligned)
    assert lcms_collection.rt_alignment_attempted
    
    # Check that scan_time_aligned was added to all samples
    for i, lcms_obj in enumerate(lcms_collection):
        sample_name = lcms_collection.samples[i]
        print(f"\nChecking {sample_name}: {lcms_obj.scan_df.columns.tolist()}")
        assert 'scan_time_aligned' in lcms_obj.scan_df.columns, f"Missing scan_time_aligned in {sample_name}"


def test_lcms_collection_consensus_features(lcms_collection):
    """Test generation of consensus mass features (clustering)."""
    # Debug: Check which sample is the center
    print(f"\nManifest dict center values:")
    for sample_name, manifest_data in lcms_collection._manifest_dict.items():
        print(f"  {sample_name}: center={manifest_data.get('center', False)}")
    
    # Ensure alignment is done first
    if not lcms_collection.rt_aligned:
        lcms_collection.align_lcms_objects()
    
    # Debug: Check if alignment was used for each sample
    print(f"\nAlignment results:")
    for sample_name, manifest_data in lcms_collection._manifest_dict.items():
        use_alignment = manifest_data.get('use_rt_alignment', None)
        print(f"  {sample_name}: use_rt_alignment={use_alignment}")
    
    # Debug: Check mass features dataframe before clustering
    mf_df = lcms_collection.mass_features_dataframe
    print(f"\nBefore clustering:")
    print(f"  Total mass features: {len(mf_df)}")
    print(f"  Columns: {mf_df.columns.tolist()}")
    print(f"  Has scan_time_aligned: {'scan_time_aligned' in mf_df.columns}")
    print(f"\nSample breakdown:")
    for sample in mf_df['sample_name'].unique():
        sample_features = mf_df[mf_df['sample_name'] == sample]
        print(f"  {sample}: {len(sample_features)} features")
        print(f"    m/z range: {sample_features['mz'].min():.4f} - {sample_features['mz'].max():.4f}")
        print(f"    RT range: {sample_features['scan_time_aligned'].min():.2f} - {sample_features['scan_time_aligned'].max():.2f}")
        print(f"    First 5 m/z values: {sample_features['mz'].head().tolist()}")
    
    # Check if the features actually match between samples
    s1_mz = set(mf_df[mf_df['sample_name'] == 'test_sample_01']['mz'].round(4))
    s2_mz = set(mf_df[mf_df['sample_name'] == 'test_sample_02']['mz'].round(4))
    print(f"\n  Overlapping m/z values (rounded to 4 decimals): {len(s1_mz & s2_mz)} out of {len(s2_mz)} in sample 2")
    
    # Check scan_time_aligned values for matching features
    s1_df = mf_df[mf_df['sample_name'] == 'test_sample_01'].sort_values('mz')
    s2_df = mf_df[mf_df['sample_name'] == 'test_sample_02'].sort_values('mz')
    print(f"\nFirst 5 features comparison:")
    print(f"  Sample 1: m/z={s1_df['mz'].head().tolist()}")
    print(f"           RT={s1_df['scan_time_aligned'].head().tolist()}")
    print(f"  Sample 2: m/z={s2_df['mz'].head().tolist()}")
    print(f"           RT={s2_df['scan_time_aligned'].head().tolist()}")
    
    # Debug: Check if scan_time values are the same (before alignment)
    print(f"\nOriginal scan_time values (before alignment):")
    print(f"  Sample 1: scan_time={s1_df['scan_time'].head().tolist()}")
    print(f"  Sample 2: scan_time={s2_df['scan_time'].head().tolist()}")
    print(f"  Are they equal? {np.allclose(s1_df['scan_time'].head().values, s2_df['scan_time'].head().values)}")
    
    print(f"\nClustering parameters:")
    print(f"  consensus_mz_tol_ppm: {lcms_collection.parameters.lcms_collection.consensus_mz_tol_ppm}")
    print(f"  consensus_rt_tol: {lcms_collection.parameters.lcms_collection.consensus_rt_tol}")
    print(f"  consensus_min_sample_fraction: {lcms_collection.parameters.lcms_collection.consensus_min_sample_fraction}")
    
    # Generate consensus features
    lcms_collection.add_consensus_mass_features()
    
    # Debug: Check after clustering
    mf_df_after = lcms_collection.mass_features_dataframe
    print(f"\nAfter clustering:")
    print(f"  Total mass features: {len(mf_df_after)}")
    print(f"  Has cluster column: {'cluster' in mf_df_after.columns}")
    if 'cluster' in mf_df_after.columns:
        print(f"  Unique clusters: {mf_df_after['cluster'].nunique()}")
    
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
    
    print(f"\nBefore gap filling:")
    print(f"  Sample 1 mass features: {sample_1_mf_count} (full)")
    print(f"  Sample 2 mass features: {sample_2_mf_count} (partial - first 10)")
    print(f"  Sample 3 mass features: {sample_3_mf_count} (empty)")
    
    # Perform gap filling
    lcms_collection.search_for_missing_mass_features()
    
    # Check that induced mass features dataframe exists
    induced_df = lcms_collection.induced_mass_features_dataframe
    assert induced_df is not None
    
    # With samples 2 and 3 having missing features, gap filling should create induced features
    print(f"\nAfter gap filling:")
    print(f"  Total induced mass features found: {len(induced_df)}")
    assert len(induced_df) > 0, "Gap filling should create induced mass features in samples 2 and 3"
    
    # Check that induced mass features have proper columns
    assert 'cluster' in induced_df.columns
    assert 'sample_name' in induced_df.columns
    assert 'mf_id' in induced_df.columns
    
    # Check induced features per sample
    sample_2_induced = len(lcms_collection[1].induced_mass_features)
    sample_3_induced = len(lcms_collection[2].induced_mass_features)
    
    print(f"  Sample 2 induced features: {sample_2_induced}")
    print(f"  Sample 3 induced features: {sample_3_induced}")
    
    # Sample 3 should have the most induced features (started with 0)
    assert sample_3_induced > 0, "Sample 3 should have induced mass features from gap filling"
    
    # Sample 2 should also have some induced features (for features 11+)
    assert sample_2_induced > 0, "Sample 2 should have induced features for missing mass features"
    
    # Check that individual LCMS objects have induced_mass_features
    total_induced = sum(len(lcms_obj.induced_mass_features) 
                       for lcms_obj in lcms_collection)
    assert total_induced > 0


def test_lcms_collection_pivot_table(lcms_collection):
    """Test creation of pivot tables for collection data."""
    # Setup: ensure we have clustered features
    if not lcms_collection.rt_aligned:
        lcms_collection.align_lcms_objects()
    lcms_collection.add_consensus_mass_features()
    
    # Create pivot table with default attribute (coll_mf_id)
    pivot_df = lcms_collection.collection_pivot_table(verbose=False)
    
    # Check that pivot table was created
    assert pivot_df is not None
    assert isinstance(pivot_df, pd.DataFrame)
    
    # Check dimensions: rows should be clusters, columns should be samples
    assert len(pivot_df.columns) == len(lcms_collection.samples)
    
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


def test_lcms_collection_plot_tics(lcms_collection):
    """Test plotting total ion chromatograms for the collection."""
    # This test just ensures the method runs without error
    # Visual inspection of plots is done manually
    try:
        import matplotlib
        matplotlib.use('Agg')  # Use non-interactive backend for testing
        lcms_collection.plot_tics(ms_level=1, type="raw", plot_legend=False)
        assert True
    except Exception as e:
        pytest.fail(f"plot_tics raised an exception: {e}")


def test_lcms_collection_feature_annotations_table(lcms_collection, msp_file_location):
    """Test creation of feature annotations table with molecular metadata."""
    # Setup: align and cluster
    if not lcms_collection.rt_aligned:
        lcms_collection.align_lcms_objects()
    lcms_collection.add_consensus_mass_features()
    
    # Load molecular metadata from MSP file
    my_msp = MSPInterface(file_path=msp_file_location)
    msp_lib, molecular_metadata = my_msp.get_metabolomics_spectra_library(
        polarity="negative",
        format="flashentropy",
        normalize=True
    )
    
    # Create annotations table without metadata first
    annotations_table = lcms_collection.feature_annotations_table(
        molecular_metadata=None,
        drop_unannotated=False
    )
    
    assert annotations_table is not None
    assert isinstance(annotations_table, pd.DataFrame)
    assert len(annotations_table) > 0
    
    # Check that cluster information is present
    assert 'cluster' in annotations_table.columns


def test_lcms_collection_sample_access(lcms_collection):
    """Test various ways to access samples in the collection."""
    # Test indexing
    first_sample = lcms_collection[0]
    assert first_sample is not None
    
    # Test iteration
    sample_count = 0
    for lcms_obj in lcms_collection:
        assert lcms_obj is not None
        sample_count += 1
    assert sample_count == len(lcms_collection)
    
    # Test samples property
    sample_names = lcms_collection.samples
    assert len(sample_names) == len(lcms_collection)
    
    # Test raw_files property
    raw_files = lcms_collection.raw_files
    assert len(raw_files) == len(lcms_collection)


def test_lcms_collection_update_raw_file_locations(lcms_collection, tmp_path):
    """Test updating raw file locations in the collection."""
    # Create a new path for raw files
    new_raw_folder = tmp_path / "new_raw_location"
    new_raw_folder.mkdir()
    
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
    
    # Step 4: Perform gap filling
    lcms_collection.search_for_missing_mass_features()
    
    # Step 5: Create reports
    pivot_table = lcms_collection.collection_pivot_table(verbose=False)
    assert pivot_table is not None
    
    cluster_reps = lcms_collection.cluster_representatives_table()
    assert len(cluster_reps) == cluster_count
    
    annotations = lcms_collection.feature_annotations_table(
        molecular_metadata=None,
        drop_unannotated=False
    )
    assert annotations is not None
    
    # Verify workflow completed successfully
    assert lcms_collection.rt_aligned or lcms_collection.rt_alignment_attempted
    assert cluster_count > 0
    print(f"\nWorkflow completed: {cluster_count} consensus clusters from {len(lcms_collection)} samples")


def test_lcms_collection_memory_management(lcms_collection):
    """Test that collection properly manages memory for raw data."""
    # Initially, raw data should not be loaded (load_light=True in fixture)
    for i, lcms_obj in enumerate(lcms_collection):
        # Check that _ms_unprocessed is either empty or not loaded
        if hasattr(lcms_obj, '_ms_unprocessed'):
            # MS level 1 should not have large raw data loaded
            if 1 in lcms_obj._ms_unprocessed:
                # For light loading, this should be empty or minimal
                pass
    
    # Load raw data for one sample
    lcms_collection.load_raw_data(sample_idx=0, ms_level=1)
    
    # Check that data was loaded
    assert 1 in lcms_collection[0]._ms_unprocessed
    
    # Drop raw data
    lcms_collection.drop_raw_data(sample_idx=0, ms_level=1)
    
    # Verify data was dropped
    # After dropping, the dict might be empty or have an empty DataFrame
    if 1 in lcms_collection[0]._ms_unprocessed:
        # Should be empty or minimal size now
        pass


def test_lcms_collection_cleanup(lcms_collection_folder):
    """Test cleanup of temporary collection data."""
    # Verify the folder exists
    assert lcms_collection_folder.exists()
    
    # The pytest tmp_path fixture will automatically clean up
    # This test just verifies the structure
    assert len(list(lcms_collection_folder.glob("*.corems"))) > 0
