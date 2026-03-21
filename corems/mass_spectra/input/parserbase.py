from abc import ABC, abstractmethod
import datetime
import numbers
from typing import Optional, Union, List, Tuple


class SpectraParserInterface(ABC):
    """
    Interface for parsing mass spectra data into MassSpectraBase objects.

    Methods
    -------
    * load().
        Load mass spectra data.
    * run().
        Parse mass spectra data.
    * get_mass_spectra_obj().
        Return MassSpectraBase object with several attributes populated
    * get_mass_spectrum_from_scan(scan_number).
        Return MassSpecBase data object from scan number.
    * get_scans_in_time_range(time_range).
        Return scan numbers within specified retention time range(s).

    Notes
    -----
    This is an abstract class and should not be instantiated directly.
    
    Time Range Filtering
    --------------------
    Many methods support optional time_range parameter to load only scans within
    specific retention time windows. This significantly improves performance for
    targeted workflows. Time ranges can be specified as:
    - Single range: (start_time, end_time) in minutes
    - Multiple ranges: [(start1, end1), (start2, end2), ...] in minutes
    """

    @abstractmethod
    def load(self):
        """
        Load mass spectra data.
        """
        pass

    @abstractmethod
    def run(self):
        """
        Parse mass spectra data.
        """
        pass

    @abstractmethod
    def get_scan_df(self, time_range: Optional[Union[Tuple[float, float], List[Tuple[float, float]]]] = None):
        """
        Return scan data as a pandas DataFrame.
        
        Parameters
        ----------
        time_range : tuple or list of tuples, optional
            Retention time range(s) to filter scans. Can be:
            - Single range: (start_time, end_time) in minutes
            - Multiple ranges: [(start1, end1), (start2, end2), ...] in minutes
            If None, returns all scans.
        
        Returns
        -------
        pd.DataFrame
            DataFrame containing scan information, optionally filtered by time range.
        """
        pass

    @abstractmethod
    def get_ms_raw(self, spectra, scan_df, time_range: Optional[Union[Tuple[float, float], List[Tuple[float, float]]]] = None):
        """
        Return a dictionary of mass spectra data as pandas DataFrames.
        
        Parameters
        ----------
        spectra : str or dict
            Specifies which spectra to load (e.g., 'ms1', 'ms2', or custom dict)
        scan_df : pd.DataFrame
            Scan information DataFrame
        time_range : tuple or list of tuples, optional
            Retention time range(s) to filter scans. Can be:
            - Single range: (start_time, end_time) in minutes
            - Multiple ranges: [(start1, end1), (start2, end2), ...] in minutes
            If None, returns all scans.
        
        Returns
        -------
        dict
            Dictionary of raw mass spectra data, optionally filtered by time range.
        """
        pass

    @abstractmethod
    def get_mass_spectra_obj(self, time_range: Optional[Union[Tuple[float, float], List[Tuple[float, float]]]] = None):
        """
        Return mass spectra data object.
        
        Parameters
        ----------
        time_range : tuple or list of tuples, optional
            Retention time range(s) to load. Can be:
            - Single range: (start_time, end_time) in minutes
            - Multiple ranges: [(start1, end1), (start2, end2), ...] in minutes
            If None, loads all scans. Useful for targeted workflows to improve performance.
        
        Returns
        -------
        MassSpectraBase
            Mass spectra object, optionally filtered to specified time range(s).
        """
        pass

    @abstractmethod
    def get_mass_spectrum_from_scan(
        self, scan_number, spectrum_mode, auto_process=True
    ):
        """
        Return mass spectrum data object from scan number.
        """
        pass

    @abstractmethod
    def get_mass_spectra_from_scan_list(
        self, scan_list, spectrum_mode, auto_process=True
    ):
        """
        Return a list of mass spectrum data objects from a list of scan numbers.
        """
        pass

    @abstractmethod
    def get_scans_in_time_range(
        self, 
        time_range: Union[Tuple[float, float], List[Tuple[float, float]]],
        ms_level: Optional[int] = None
    ) -> List[int]:
        """
        Return scan numbers within specified retention time range(s).
        
        This method provides efficient filtering of scans by retention time,
        which is particularly useful for targeted workflows where only specific
        time windows are of interest.
        
        Parameters
        ----------
        time_range : tuple or list of tuples
            Retention time range(s) in minutes. Can be:
            - Single range: (start_time, end_time)
            - Multiple ranges: [(start1, end1), (start2, end2), ...]
        ms_level : int, optional
            If specified, only return scans of this MS level (e.g., 1 for MS1, 2 for MS2).
            If None, returns scans of all MS levels.
        
        Returns
        -------
        list of int
            List of scan numbers within the specified time range(s) and MS level.
        
        Examples
        --------
        Get MS1 scans between 1.0 and 2.0 minutes:
        
        >>> scans = parser.get_scans_in_time_range((1.0, 2.0), ms_level=1)
        
        Get scans in multiple time windows:
        
        >>> scans = parser.get_scans_in_time_range([(0.5, 1.5), (3.0, 4.0)])
        """
        pass

    @abstractmethod
    def get_instrument_info(self):
        """
        Return instrument information.

        Returns
        -------
        dict
            A dictionary with the keys 'model', and 'serial_number'.
        """
        pass

    @abstractmethod
    def get_creation_time(self) -> datetime.datetime:
        """
        Return the creation time of the mass spectra data.

        Returns
        -------
        datetime.datetime
            The creation time of the mass spectra data.
        """
        pass
    
    @staticmethod
    def _normalize_time_range(
        time_range: Optional[Union[Tuple[float, float], List[Tuple[float, float]]]]
    ) -> Optional[List[Tuple[float, float]]]:
        """
        Normalize time range input to a list of tuples.
        
        Helper method for implementations to standardize time_range parameter.
        Converts single tuple to list of tuples for consistent handling.
        
        Parameters
        ----------
        time_range : tuple, list of tuples, or None
            Input time range(s)
        
        Returns
        -------
        list of tuples or None
            Normalized time ranges as list of (start, end) tuples, or None if input is None.
        
        Examples
        --------
        >>> SpectraParserInterface._normalize_time_range((1.0, 2.0))
        [(1.0, 2.0)]
        
        >>> SpectraParserInterface._normalize_time_range([(1.0, 2.0), (3.0, 4.0)])
        [(1.0, 2.0), (3.0, 4.0)]
        
        >>> SpectraParserInterface._normalize_time_range(None)
        None
        """
        if time_range is None:
            return None
        
        # Check if it's a single tuple (two numbers)
        if isinstance(time_range, tuple) and len(time_range) == 2:
            # Use numbers.Number to catch int, float, and numpy scalar types
            if isinstance(time_range[0], numbers.Number) and isinstance(time_range[1], numbers.Number):
                # Convert to float to ensure consistency (handles numpy scalars)
                return [(float(time_range[0]), float(time_range[1]))]
        
        # Otherwise assume it's already a list of tuples
        return list(time_range)
