"""
Sample-level operations for LCMS collection processing pipelines.

This module provides a framework for defining reusable, composable operations
that can be executed on individual samples in a parallelized manner.

Classes
-------
SampleOperation
    Base class for all sample-level operations
GapFillOperation
    Gap-fill missing cluster features for a sample
ReloadFeaturesOperation
    Reload mass features from HDF5 for a sample

"""

import pandas as pd


class SampleOperation:
    """
    Base class for operations that can be performed on a sample.
    
    All sample operations should inherit from this class and implement
    the execute() method. Optionally override can_execute() for conditional
    execution and collect_results() for custom result collection.
    
    Parameters
    ----------
    name : str
        Name of the operation (for logging and identification)
    **kwargs
        Additional keyword arguments stored as operation parameters
        
    Attributes
    ----------
    name : str
        Operation name
    params : dict
        Dictionary of operation parameters
    """
    
    def __init__(self, name, **kwargs):
        self.name = name
        self.params = kwargs
        
    def needs_raw_ms_data(self):
        """
        Declare whether this operation needs raw MS data loaded.
        
        Override this method to specify raw data requirements. The pipeline
        executor will ensure raw data is loaded before executing operations
        that need it, and can clean it up afterwards.
        
        Returns
        -------
        tuple of (bool, int or None)
            (needs_raw_data, ms_level)
            - needs_raw_data: True if operation needs raw MS data
            - ms_level: MS level needed (1 for MS1, 2 for MS2, etc.) or None
        """
        return False, None
        
    def can_execute(self, sample, collection):
        """
        Check if this operation can be executed on the sample.
        
        Parameters
        ----------
        sample : LCMSBase
            The sample to check
        collection : LCMSBaseCollection
            The collection containing the sample
            
        Returns
        -------
        bool
            True if operation can execute, False otherwise
        """
        return True
        
    def execute(self, sample_id, collection, **runtime_params):
        """
        Execute the operation on a sample.
        
        This method must be implemented by subclasses.
        
        Parameters
        ----------
        sample_id : int
            Sample ID to process
        collection : LCMSBaseCollection
            The collection containing the sample
        **runtime_params
            Runtime parameters passed from the pipeline
            
        Returns
        -------
        result
            Operation result (can be None if operation modifies sample in place)
        """
        raise NotImplementedError(f"execute() not implemented for {self.__class__.__name__}")
        
    def collect_results(self, sample_id, result, collection):
        """
        Collect results back into collection after parallel execution.
        
        Override this method if the operation returns results that need
        to be collected back into the collection object.
        
        Parameters
        ----------
        sample_id : int
            Sample ID that was processed
        result
            Result returned from execute()
        collection : LCMSBaseCollection
            The collection to update
        """
        pass
        
    def __repr__(self):
        return f"{self.__class__.__name__}(name='{self.name}')"


class GapFillOperation(SampleOperation):
    """
    Gap-fill missing cluster features for a sample.
    
    Searches raw MS1 data to find peaks in expected m/z and retention time
    windows for clusters that are present in other samples but missing from
    this sample.
    
    Parameters
    ----------
    name : str
        Operation name
    expand_on_miss : bool, optional
        If True, expands search window when no peak is found. Default is False.
        
    Notes
    -----
    Requires that add_consensus_mass_features() has been run on the collection.
    This operation loads raw MS1 data which will be available for subsequent operations.
    """
    
    def needs_raw_ms_data(self):
        """This operation needs raw MS1 data."""
        return True, 1
    
    def can_execute(self, sample, collection):
        """Check if cluster summary exists."""
        return hasattr(collection, 'cluster_summary_dataframe') and \
               collection.cluster_summary_dataframe is not None
    
    def execute(self, sample_id, collection, missingdf, cluster_dict, expand_on_miss, **runtime_params):
        """
        Execute gap-filling for a single sample.
        
        Parameters
        ----------
        sample_id : int
            Sample index to process
        collection : LCMSBaseCollection
            The collection
        missingdf : pd.DataFrame
            DataFrame with cluster information and missing samples
        cluster_dict : dict
            Cluster feature dictionary
        expand_on_miss : bool
            Whether to expand search window on miss
        **runtime_params
            Additional runtime parameters (ignored)
            
        Returns
        -------
        dict
            Dictionary of induced mass features
        """
        # This is essentially the same logic as _search_for_targeted_mass_features_in_sample
        # but extracted into an operation
        
        # Get clusters missing data for this sample
        sampledf = missingdf[
            missingdf.missing_samples.apply(lambda x: sample_id in x)
        ].reset_index(drop=True).copy()

        # Skip if no missing features for this sample
        if len(sampledf) == 0:
            return {}

        # Load raw data for this sample
        collection.load_raw_data(sample_id, 1)
        
        # Get MS1 data
        ms1df = collection[sample_id]._ms_unprocessed[1].copy()
        scan_df = collection[sample_id].scan_df[['scan', 'scan_time_aligned']]
        ms1df = pd.merge(ms1df, scan_df, on='scan')

        # Pre-extract all values from sampledf
        clusters = sampledf.cluster.values
        mz_mins = sampledf.mz_min.values
        mz_maxs = sampledf.mz_max.values
        st_mins = sampledf.scan_time_aligned_min.values
        st_maxs = sampledf.scan_time_aligned_max.values
        
        if expand_on_miss:
            mz_mins_allowed = sampledf.mz_min_allowed.values
            mz_maxs_allowed = sampledf.mz_max_allowed.values
            st_mins_allowed = sampledf.sta_min_allowed.values
            st_maxs_allowed = sampledf.sta_max_allowed.values

        # Pre-filter ms1df to reduce search space
        mz_global_min = mz_mins.min()
        mz_global_max = mz_maxs.max()
        st_global_min = st_mins.min()
        st_global_max = st_maxs.max()
        
        if expand_on_miss:
            mz_global_min = min(mz_global_min, mz_mins_allowed.min())
            mz_global_max = max(mz_global_max, mz_maxs_allowed.max())
            st_global_min = min(st_global_min, st_mins_allowed.min())
            st_global_max = max(st_global_max, st_maxs_allowed.max())
        
        ms1df_filtered = ms1df[
            (ms1df.mz >= mz_global_min) & 
            (ms1df.mz <= mz_global_max) &
            (ms1df.scan_time_aligned >= st_global_min) &
            (ms1df.scan_time_aligned <= st_global_max)
        ].copy()

        # Generate set_ids for all features
        set_ids = [f'c{clusters[i]}_{i}_i' for i in range(len(sampledf))]
        
        # Use batch method to process all features at once
        if expand_on_miss:
            # First try with normal bounds
            peaks_dict = collection[sample_id].search_for_targeted_mass_features_batch(
                ms1df_filtered,
                mz_mins,
                mz_maxs,
                st_mins,
                st_maxs,
                set_ids,
                obj_idx=sample_id,
                st_aligned=True
            )
            
            # Retry failed features with expanded bounds
            failed_indices = [i for i, sid in enumerate(set_ids) if peaks_dict[sid].apex_scan == -99]
            if failed_indices:
                failed_ids = [set_ids[i] for i in failed_indices]
                retry_peaks = collection[sample_id].search_for_targeted_mass_features_batch(
                    ms1df_filtered,
                    mz_mins_allowed[failed_indices],
                    mz_maxs_allowed[failed_indices],
                    st_mins_allowed[failed_indices],
                    st_maxs_allowed[failed_indices],
                    failed_ids,
                    obj_idx=sample_id,
                    st_aligned=True
                )
                peaks_dict.update(retry_peaks)
        else:
            peaks_dict = collection[sample_id].search_for_targeted_mass_features_batch(
                ms1df_filtered,
                mz_mins,
                mz_maxs,
                st_mins,
                st_maxs,
                set_ids,
                obj_idx=sample_id,
                st_aligned=True
            )
        
        # Build induced_mass_features dict and update cluster_dict
        induced_mass_features = {}
        for i in range(len(sampledf)):
            peak = peaks_dict[set_ids[i]]
            induced_mass_features[peak.id] = peak
            cluster_dict[clusters[i]] += [set_ids[i]]
        
        # Integrate mass features (don't fail on bad integration)
        collection[sample_id].induced_mass_features = induced_mass_features
        collection[sample_id].integrate_mass_features(drop_if_fail=False, induced_features=True)
        
        # Return the induced features
        return collection[sample_id].induced_mass_features
        
    def collect_results(self, sample_id, result, collection):
        """Collect induced mass features back into sample."""
        collection[sample_id].induced_mass_features = result


class ReloadFeaturesOperation(SampleOperation):
    """
    Reload mass features from HDF5 and optionally add MS1/MS2 spectra.
    
    This is useful when the collection was loaded with load_light=True,
    which stores mass features only in the collection dataframe and not
    as LCMSMassFeature objects in individual samples.
    
    Parameters
    ----------
    name : str
        Operation name
    add_ms1 : bool, optional
        If True, adds MS1 spectra to mass features. Automatically uses raw MS1 data
        if available (e.g., from gap-filling), otherwise uses parser. Spectrum mode
        is auto-detected. Default is False.
    add_ms2 : bool, optional
        If True, also loads and associates MS2 spectra. Spectrum mode is auto-detected.
        Default is False.
    auto_process_ms2 : bool, optional
        If True and add_ms2=True, auto-processes MS2 spectra. Default is True.
    ms2_scan_filter : str or None, optional
        Filter string for MS2 scans. Default is None.
        
    Notes
    -----
    MS1 spectra association automatically uses raw MS1 data if loaded by a previous
    operation (e.g., GapFillOperation). This is efficient when multiple operations
    need MS1 data in the same pipeline. All spectrum modes are auto-detected from
    the data.
    """
    
    def can_execute(self, sample, collection):
        """Check if collection parser is available."""
        return hasattr(collection, 'collection_parser') and \
               collection.collection_parser is not None
    
    def execute(self, sample_id, collection, mf_ids_to_load=None, **runtime_params):
        """
        Execute feature reloading for a single sample.
        
        Parameters
        ----------
        sample_id : int
            Sample ID to reload features for
        collection : LCMSBaseCollection
            The collection
        mf_ids_to_load : list of str, optional
            List of collection-level mf_ids to load
        **runtime_params
            Additional runtime parameters (ignored)
            
        Returns
        -------
        dict
            Dictionary of reloaded mass features
        """
        # Get parameters
        add_ms1 = self.params.get('add_ms1', False)
        add_ms2 = self.params.get('add_ms2', False)
        auto_process_ms2 = self.params.get('auto_process_ms2', True)
        ms2_scan_filter = self.params.get('ms2_scan_filter', None)
        
        sample = collection[sample_id]
        sample_name = collection.samples[sample_id]
        
        # Auto-determine if we should use parser for MS1 (check if raw data is available)
        has_raw_ms1 = 1 in sample._ms_unprocessed and not sample._ms_unprocessed[1].empty
        use_parser_for_ms1 = not has_raw_ms1  # Use parser only if raw data not available
        
        # Spectrum modes will be auto-detected (None = auto-detect)
        spectrum_mode_ms1 = None
        ms2_spectrum_mode = None
        
        # Check if we have a collection parser
        if not hasattr(collection, 'collection_parser') or collection.collection_parser is None:
            print(f"Warning: Cannot reload mass features for {sample_name} - no collection_parser available")
            return {}
        
        # Get the HDF5 file for this sample
        hdf5_file = collection.collection_parser.folder_location / f"{sample_name}.corems/{sample_name}.hdf5"
        
        if not hdf5_file.exists():
            print(f"Warning: HDF5 file not found for sample {sample_name}: {hdf5_file}")
            return {}
        
        # Import here to avoid circular imports
        from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectra
        
        # If specific mf_ids requested, extract the local mf_ids we need
        local_mf_ids_to_load = None
        if mf_ids_to_load is not None:
            local_mf_ids_to_load = set()
            for coll_mf_id in mf_ids_to_load:
                # Parse collection-level ID to get local ID
                parts = str(coll_mf_id).split('_', 1)
                if len(parts) == 2:
                    try:
                        local_mf_ids_to_load.add(int(parts[1]))
                    except ValueError:
                        local_mf_ids_to_load.add(parts[1])
        
        # Reload mass features from HDF5
        with ReadCoreMSHDFMassSpectra(hdf5_file) as parser:
            parser.import_mass_features(sample, mf_ids=local_mf_ids_to_load)
        
        # If add_ms1, associate MS1 spectra with the loaded mass features
        if add_ms1 and len(sample.mass_features) > 0:
            # Check if raw MS1 data is already loaded (e.g., from gap-filling)
            has_raw_ms1 = 1 in sample._ms_unprocessed and not sample._ms_unprocessed[1].empty
            
            if has_raw_ms1 and not use_parser_for_ms1:
                # Use already-loaded raw data (more efficient)
                sample.add_associated_ms1(
                    auto_process=True,
                    use_parser=False,
                    spectrum_mode=spectrum_mode_ms1
                )
            else:
                # Use parser to get MS1 spectra
                sample.add_associated_ms1(
                    auto_process=True,
                    use_parser=True,
                    spectrum_mode=spectrum_mode_ms1
                )
        
        # If add_ms2, associate MS2 spectra with the loaded mass features
        if add_ms2 and local_mf_ids_to_load is not None:
            collection._associate_ms2_with_mass_features(
                sample, 
                list(local_mf_ids_to_load),
                auto_process=auto_process_ms2,
                spectrum_mode=ms2_spectrum_mode,
                scan_filter=ms2_scan_filter
            )
        
        return sample.mass_features
        
    def collect_results(self, sample_id, result, collection):
        """Collect reloaded mass features back into sample."""
        collection[sample_id].mass_features = result


class MolecularFormulaSearchOperation(SampleOperation):
    """
    Perform molecular formula search on mass features using associated MS1 spectra.
    
    This operation runs molecular formula search on all mass features in a sample
    that have associated MS1 spectra. Requires MS1 spectra to be loaded and
    processed before execution.
    
    Parameters
    ----------
    name : str
        Operation name (for logging)
    **kwargs
        Additional parameters passed to parent class
        
    Examples
    --------
    >>> op = MolecularFormulaSearchOperation('mf_search')
    >>> # Use in pipeline
    >>> results = collection.process_samples_pipeline([op])
    
    Notes
    -----
    This operation requires that MS1 spectra have been associated with mass
    features (e.g., via ReloadFeaturesOperation with add_ms1=True). The
    molecular formula search uses parameters from the collection's 
    parameters.mass_spectrum["ms1"].molecular_search settings.
    """
    
    def __init__(self, name='molecular_formula_search', **kwargs):
        super().__init__(name, **kwargs)
    
    def needs_raw_ms_data(self):
        """
        This operation doesn't need raw data - it works on processed MS1 spectra
        that are already associated with mass features.
        
        Returns
        -------
        tuple
            (False, None) - no raw data needed
        """
        return False, None
    
    def can_execute(self, sample, collection, **runtime_params):
        """
        Check if molecular formula search can be executed.
        
        Requires that the sample has mass features with associated MS1 spectra.
        
        Parameters
        ----------
        sample : LCMSObject
            The sample object
        collection : LCMSCollection
            The collection containing the sample
        **runtime_params
            Runtime parameters (not used)
            
        Returns
        -------
        bool
            True if sample has mass features with MS1 spectra
        """        
        # Check if sample has mass features
        if not hasattr(sample, 'mass_features') or not sample.mass_features:
            return False
        
        # Check if at least some mass features have MS1 spectra
        has_ms1 = any(
            hasattr(mf, 'mass_spectrum') and mf.mass_spectrum is not None
            for mf in sample.mass_features.values()
        )
        
        return has_ms1
    
    def execute(self, sample_id, collection, **runtime_params):
        """
        Execute molecular formula search on a sample.
        
        Creates a SearchMolecularFormulasLC object and runs mass feature search,
        which annotates mass features with molecular formula assignments.
        
        Parameters
        ----------
        sample_id : str
            Sample identifier
        collection : LCMSCollection
            The collection containing the sample
        **runtime_params
            Runtime parameters (not used)
            
        Returns
        -------
        int
            Number of mass features that were searched
        """
        from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulasLC
        import time
        import sqlalchemy.exc
        import sqlite3
        
        sample = collection[sample_id]
        
        # Verify that mass features exist
        if not hasattr(sample, 'mass_features') or not sample.mass_features:
            return 0  # No mass features to search
        
        # Verify that mass features have MS1 spectra associated
        if not hasattr(sample, '_ms') or not sample._ms:
            raise RuntimeError(
                f"Sample {sample_id} does not have MS1 spectra loaded in _ms dictionary. "
                "Molecular formula search requires MS1 spectra to be associated with mass features. "
                "Ensure add_ms1=True when reloading features."
            )
        
        # Prepare data for bulk molecular formula search
        # Group mass features by their apex scan
        scan_to_mf = {}
        for mf_id, mf in sample.mass_features.items():
            apex_scan = mf.apex_scan
            if apex_scan not in scan_to_mf:
                scan_to_mf[apex_scan] = []
            scan_to_mf[apex_scan].append(mf)
        
        # Build lists of mass spectra and corresponding peaks
        mass_spectrum_list = []
        ms_peaks_list = []
        
        for scan_num, mf_list in scan_to_mf.items():
            # Get the mass spectrum for this scan
            if scan_num not in sample._ms:
                continue  # Skip if spectrum not loaded
                
            mass_spectrum = sample._ms[scan_num]
            
            # Verify spectrum is processed (has peaks)
            if not hasattr(mass_spectrum, '_mspeaks') or not mass_spectrum._mspeaks:
                continue  # Skip unprocessed spectra
            
            # Get the MS1 peaks for each mass feature at this scan
            peaks_for_scan = []
            for mf in mf_list:
                try:
                    # Use the ms1_peak property which finds the closest peak
                    ms1_peak = mf.ms1_peak
                    peaks_for_scan.append(ms1_peak)
                except (AttributeError, IndexError):
                    # Skip if ms1_peak can't be determined
                    continue
            
            if peaks_for_scan:
                mass_spectrum_list.append(mass_spectrum)
                ms_peaks_list.append(peaks_for_scan)
        
        # Run molecular formula search if we have data, with retry logic for database locks
        if mass_spectrum_list and ms_peaks_list:
            max_retries = 10
            retry_delay = 2  # seconds
            
            for attempt in range(max_retries):
                try:
                    mol_search = SearchMolecularFormulasLC(sample)
                    mol_search.bulk_run_molecular_formula_search(mass_spectrum_list, ms_peaks_list)
                    break  # Success, exit retry loop
                except (sqlalchemy.exc.OperationalError, sqlite3.OperationalError) as e:
                    if attempt < max_retries - 1:
                        # Database is locked, retry after delay
                        print(f"Sample {sample_id}: Database locked during molecular formula search, retrying in {retry_delay}s (attempt {attempt + 1}/{max_retries})...")
                        time.sleep(retry_delay)
                    else:
                        # Max retries exceeded, re-raise the exception
                        raise RuntimeError(
                            f"Sample {sample_id}: Molecular formula search failed after {max_retries} attempts due to database lock. "
                            "Try reducing parallel cores or increasing database timeout."
                        ) from e
        
        # Return count of features searched
        return len(sample.mass_features)
    
    def collect_results(self, sample_id, result, collection):
        """
        Collect results (no-op as search modifies mass features in place).
        
        The molecular formula search modifies mass features in place by adding
        molecular formula assignments, so no explicit result collection is needed.
        
        Parameters
        ----------
        sample_id : str
            Sample identifier
        result : int
            Number of features searched
        collection : LCMSCollection
            The collection containing the sample
        """
        # Search modifies mass features in place, nothing to collect
        pass
