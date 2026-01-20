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

from abc import ABC, abstractmethod
import pandas as pd


class SampleOperation(ABC):
    """
    Base class for operations that can be performed on a sample.
    
    All sample operations must inherit from this class and implement all
    abstract methods. This ensures proper integration with the pipeline framework.
    
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
    description : str
        Human-readable description for progress messages (must override in subclasses)
    """
    
    def __init__(self, name, **kwargs):
        self.name = name
        self.params = kwargs
    
    @property
    @abstractmethod
    def description(self):
        """
        Human-readable description for progress messages.
        
        This property must be overridden in subclasses to provide a meaningful
        description that will be shown in progress bars (e.g., "gap-filling",
        "reloading features", etc.).
        
        Returns
        -------
        str
            Brief description of what this operation does
        """
        pass
    
    @abstractmethod
    def needs_raw_ms_data(self):
        """
        Declare whether this operation needs raw MS data loaded.
        
        Subclasses must implement this method to specify raw data requirements.
        The pipeline executor will ensure raw data is loaded before executing
        operations that need it, and can clean it up afterwards.
        
        Returns
        -------
        tuple of (bool, int or None)
            (needs_raw_data, ms_level)
            - needs_raw_data: True if operation needs raw MS data
            - ms_level: MS level needed (1 for MS1, 2 for MS2, etc.) or None
            
        Examples
        --------
        >>> def needs_raw_ms_data(self):
        ...     return True, 1  # Needs MS1 data
        >>> def needs_raw_ms_data(self):
        ...     return False, None  # No raw data needed
        """
        pass
    
    @abstractmethod
    def can_execute(self, sample, collection):
        """
        Check if this operation can be executed on the sample.
        
        Subclasses must implement this method to define prerequisites.
        Return True if the operation can execute, False otherwise.
        
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
            
        Examples
        --------
        >>> def can_execute(self, sample, collection):
        ...     return True  # Can always execute
        >>> def can_execute(self, sample, collection):
        ...     return hasattr(sample, 'mass_features') and sample.mass_features
        """
        pass
    
    @abstractmethod
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
        pass
    
    @abstractmethod
    def collect_results(self, sample_id, result, collection):
        """
        Collect results back into collection after parallel execution.
        
        Subclasses must implement this method to handle result collection.
        If the operation modifies samples in place and doesn't need to collect
        results, simply implement as `pass`.
        
        Parameters
        ----------
        sample_id : int
            Sample ID that was processed
        result
            Result returned from execute()
        collection : LCMSBaseCollection
            The collection to update
            
        Examples
        --------
        >>> def collect_results(self, sample_id, result, collection):
        ...     pass  # Operation modifies sample in place
        >>> def collect_results(self, sample_id, result, collection):
        ...     collection[sample_id].induced_mass_features = result
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
    
    @property
    def description(self):
        """Human-readable description for progress messages."""
        return "gap-filling"
    
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
    
    @property
    def description(self):
        """Human-readable description for progress messages."""
        return "reloading features"
    
    def needs_raw_ms_data(self):
        """This operation doesn't need raw data."""
        return False, None
    
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
            # mf_ids_to_load is already a list of sample-level mf_ids (integers)
            # No parsing needed - they come from the mf_id column in the dataframe
            if len(mf_ids_to_load) == 0:
                # No features to load for this sample - return empty dict
                return {}
            local_mf_ids_to_load = set(mf_ids_to_load)
        
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
        """
        Collect reloaded mass features back into sample.
        
        This operation loads a subset of mass features (e.g., representatives)
        into sample.mass_features for processing, while preserving the full
        mass_features_dataframe at the collection level. Sets a lock flag to
        prevent automatic rebuilding of the collection dataframe from individual
        samples.
        
        Parameters
        ----------
        sample_id : int
            Sample ID that was processed
        result : dict
            Dictionary of reloaded mass features
        collection : LCMSBaseCollection
            The collection
        """
        # Update sample.mass_features with loaded features
        collection[sample_id].mass_features = result
        # Lock the collection dataframe to prevent rebuilding from individual samples
        # (since we've only loaded a subset, rebuilding would lose data)
        collection._mass_features_locked = True


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
    
    @property
    def description(self):
        """Human-readable description for progress messages."""
        return "molecular formula search"
    
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


class MS2SpectralSearchOperation(SampleOperation):
    """
    Perform MS2 spectral search using entropy-based matching.
    
    This operation performs spectral library search on MS2 spectra associated
    with mass features using FlashEntropy for fast similarity scoring. Requires
    MS2 spectra to be loaded and processed before execution.
    
    Parameters
    ----------
    name : str
        Operation name (for logging)
    ms2_scan_filter : str or None, optional
        Filter string for MS2 scans (e.g., 'hcd'). If None, uses all MS2 scans.
        Default is None.
    peak_sep_da : float, optional
        Peak separation in Daltons for spectral matching. Default is 0.01.
    **kwargs
        Additional parameters passed to parent class
        
    Examples
    --------
    >>> op = MS2SpectralSearchOperation('ms2_search', ms2_scan_filter='hcd')
    >>> # Use in pipeline - requires fe_lib in runtime_params
    >>> results = collection.process_samples_pipeline([op])
    
    Notes
    -----
    This operation requires:
    - MS2 spectra to be associated with mass features
    - FlashEntropy library (fe_lib) to be provided in runtime_params
    - MS2 spectra must be processed (centroided)
    
    The spectral search modifies mass features in place by adding spectral
    match scores and metadata.
    """
    
    @property
    def description(self):
        """Human-readable description for progress messages."""
        return "MS2 spectral search"
    
    def __init__(self, name='ms2_spectral_search', ms2_scan_filter=None, **kwargs):
        super().__init__(name, **kwargs)
        self.params['ms2_scan_filter'] = ms2_scan_filter
    
    def needs_raw_ms_data(self):
        """
        This operation doesn't need raw data - it works on processed MS2 spectra
        that are already associated with mass features.
        
        Returns
        -------
        tuple
            (False, None) - no raw data needed
        """
        return False, None
    
    def can_execute(self, sample, collection, **runtime_params):
        """
        Check if MS2 spectral search can be executed.
        
        Requires that the sample has mass features with MS2 spectra associated.
        
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
            True if sample has mass features with MS2 spectra
        """        
        # Check if sample has mass features
        if not hasattr(sample, 'mass_features') or not sample.mass_features:
            return False
        
        # Check if any mass features have MS2 spectra associated
        has_ms2 = any(
            hasattr(mf, 'ms2_mass_spectra') and mf.ms2_mass_spectra
            for mf in sample.mass_features.values()
        )
        
        return has_ms2
    
    def execute(self, sample_id, collection, fe_lib=None, molecular_metadata=None, **runtime_params):
        """
        Execute MS2 spectral search on a sample.
        
        Performs entropy-based spectral library search on all MS2 spectra
        in the sample that match the scan filter criteria.
        
        Parameters
        ----------
        sample_id : str
            Sample identifier
        collection : LCMSCollection
            The collection containing the sample
        fe_lib : FlashEntropy library
            Pre-computed FlashEntropy library for spectral matching
        molecular_metadata : pd.DataFrame, optional
            Metadata for molecules in the spectral library
        **runtime_params
            Runtime parameters (not used)
            
        Returns
        -------
        int
            Number of MS2 spectra searched
        """
        sample = collection[sample_id]
        
        # Get parameters
        ms2_scan_filter = self.params.get('ms2_scan_filter')
        
        # Verify that we have a spectral library
        if fe_lib is None:
            raise ValueError(
                f"Sample {sample_id}: MS2 spectral search requires fe_lib (FlashEntropy library) "
                "to be provided in runtime parameters. Create the library at the collection level "
                "and pass it to the pipeline."
            )
        
        # Extract peak_sep_da from FlashEntropy library configuration
        peak_sep_da = fe_lib.entropy_search.max_ms2_tolerance_in_da
        if peak_sep_da is None:
            raise ValueError(
                f"Sample {sample_id}: Could not extract max_ms2_tolerance_in_da from FlashEntropy library. "
                "Ensure the library was created with this parameter specified."
            )
        
        # Verify that sample has _ms dictionary
        if not hasattr(sample, '_ms') or not sample._ms:
            return 0  # No MS2 spectra to search
        
        # Get MS2 scan numbers based on filter
        if ms2_scan_filter is not None:
            # Filter by scan text
            ms2_scan_df = sample.scan_df[
                sample.scan_df.scan_text.str.contains(ms2_scan_filter) &
                (sample.scan_df.ms_level == 2)
            ]
        else:
            # All MS2 scans
            ms2_scan_df = sample.scan_df[sample.scan_df.ms_level == 2]
        
        # Get scans that are actually loaded in _ms
        ms2_scans_to_search = [
            scan for scan in ms2_scan_df.scan.tolist()
            if scan in sample._ms.keys()
        ]
        
        if not ms2_scans_to_search:
            return 0  # No MS2 spectra to search
        
        # Perform spectral search using the sample's fe_search method
        sample.fe_search(
            scan_list=ms2_scans_to_search,
            fe_lib=fe_lib,
            peak_sep_da=peak_sep_da
        )
        
        # Return the spectral search results for collection
        # (needed for multiprocessing - results populated in worker need to be returned)
        return sample.spectral_search_results
    
    def collect_results(self, sample_id, result, collection):
        """
        Collect spectral search results back into the sample.
        
        In multiprocessing, the worker's modifications don't persist to the
        main process, so we need to explicitly collect and reassign the results.
        This also re-associates the results with mass features.
        
        Parameters
        ----------
        sample_id : str
            Sample identifier
        result : dict
            Dictionary of spectral search results from execute()
        collection : LCMSCollection
            The collection containing the sample
        """
        # Assign the spectral search results back to the sample
        if result:
            collection[sample_id].spectral_search_results.update(result)
            
            # Re-associate results with mass features (same logic as fe_search)
            sample = collection[sample_id]
            if len(sample.mass_features) > 0:
                for mass_feature_id, mass_feature in sample.mass_features.items():
                    scan_ids = mass_feature.ms2_scan_numbers
                    for ms2_scan_id in scan_ids:
                        precursor_mz = mass_feature.mz
                        try:
                            sample.spectral_search_results[ms2_scan_id][precursor_mz]
                        except KeyError:
                            pass
                        else:
                            sample.mass_features[
                                mass_feature_id
                            ].ms2_similarity_results.append(
                                sample.spectral_search_results[ms2_scan_id][precursor_mz]
                            )


class LoadEICsOperation(SampleOperation):
    """
    Load extracted ion chromatograms (EICs) from HDF5 for regular mass features.
    
    Loads EICs for regular mass features that belong to consensus clusters from HDF5.
    Induced (gap-filled) features already have EICs from integrate_mass_features,
    so no additional loading is needed for them.
    
    This operation enables downstream visualization and analysis of chromatographic
    peaks across all samples in a cluster.
    
    Notes
    -----
    Requires that mass features have been loaded and cluster_index assigned.
    Regular mass feature EICs must have been previously saved to HDF5 with export_eics=True.
    Induced mass features already have EICs populated during gap-filling.
    """
    
    @property
    def description(self):
        """Human-readable description for progress messages."""
        return "loading EICs"
    
    def needs_raw_ms_data(self):
        """This operation doesn't need raw data - induced features already have EICs."""
        return False, None
    
    def can_execute(self, sample, collection):
        """
        Check if EIC loading can be executed.
        
        This operation can always execute if the sample exists - the actual work
        is determined by cluster_mz_dict in runtime_params. If cluster_mz_dict is
        empty or None, execute() will simply return 0 (no EICs loaded).
        
        Returns
        -------
        bool
            True (always executable - runtime_params control actual work)
        """
        return True
    
    def execute(self, sample_id, collection, cluster_mz_dict=None, **runtime_params):
        """
        Load EICs from HDF5 for a single sample.
        
        Loads EICs for regular mass features that belong to consensus clusters.
        Induced (gap-filled) mass features already have EICs from integrate_mass_features,
        so no additional loading is needed for them.
        
        The cluster_mz_dict parameter (passed from collection level) maps sample_id
        to a list of m/z values that belong to clusters for that sample.
        
        Parameters
        ----------
        sample_id : int
            Sample index to process
        collection : LCMSBaseCollection
            The collection
        cluster_mz_dict : dict, optional
            Dictionary mapping sample_id to list of m/z values in clusters for that sample.
            If None, will not load any EICs. Default is None.
        **runtime_params
            Additional runtime parameters (ignored)
            
        Returns
        -------
        dict
            Dictionary of loaded EIC_Data objects, keyed by m/z value
        """
        from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectra
        
        sample = collection[sample_id]
        sample_name = collection.samples[sample_id]
        
        # If no cluster info provided or no m/z values for this sample, return early
        if cluster_mz_dict is None or sample_id not in cluster_mz_dict:
            return {}
        
        # Get m/z values for this sample that belong to clusters
        sample_cluster_mz = set(cluster_mz_dict[sample_id])
        
        # Load EICs for each of the sample_cluster_mz
        hdf5_path = sample.file_location
        if hdf5_path and hdf5_path.exists():
            try:
                reader = ReadCoreMSHDFMassSpectra(str(hdf5_path))
                reader.import_eics(sample, mz_list=list(sample_cluster_mz))
                # Return the loaded EICs for multiprocessing collection
                # (modifications in worker process don't persist to main process)
                return sample.eics.copy()
            except (KeyError, AttributeError):
                # No EIC data in HDF5 file for these m/z values
                return {}
        
        return {}
    
    def collect_results(self, sample_id, result, collection):
        """
        Collect loaded EICs back into sample.
        
        In multiprocessing, the worker's modifications don't persist to the
        main process, so we need to explicitly collect and reassign the EICs.
        This also re-associates EICs with mass features.
        
        Parameters
        ----------
        sample_id : int
            Sample ID that was processed
        result : dict
            Dictionary of EIC_Data objects keyed by m/z, returned from execute()
        collection : LCMSBaseCollection
            The collection
        """
        if result:
            # Update sample.eics with loaded EICs
            collection[sample_id].eics.update(result)
            
            # Re-associate EICs with mass features (same logic as import_eics)
            sample = collection[sample_id]
            for idx in sample.mass_features.keys():
                mz = sample.mass_features[idx].mz
                if mz in sample.eics.keys():
                    sample.mass_features[idx]._eic_data = sample.eics[mz]
