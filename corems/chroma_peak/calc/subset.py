# This file contains functions for subsetting dataframes that contain mass feature data.
# This is based on the deimos package, found here: https://github.com/pnnl/deimos/blob/master/deimos/subset.py with some modifications.

import multiprocessing as mp
from functools import partial

import numpy as np
import pandas as pd

class MultiSamplePartitions:
    '''
    Generator object that will lazily build and return each partition constructed
    from multiple samples.

    Attributes
    ----------
    features : :obj:`~pandas.DataFrame`
        Input feature coordinates and intensities.
    split_on : str
        Dimension to partition the data.
    size : int
        Target partition size.
    tol : float
        Largest allowed distance between unique `split_on` observations.
    n_partitions : int
        Number of partitions in the data.

    '''

    def __init__(self,
                 features, 
                 split_on: str = 'mz', 
                 size: int = 500, 
                 tol: float = 25E-6, 
                 relative: bool = False):
        '''
        Initialize :obj:`~deimos.subset.Partitions` instance.

        Parameters
        ----------
        features : :obj:`~pandas.DataFrame`
            Input feature coordinates and intensities.
        split_on : str
            Dimension to partition the data.
        size : int
            Target partition size.
        tol : float
            Largest allowed distance between unique `split_on` observations.

        '''
        if not isinstance(split_on, str):
            raise TypeError(f"Expected 'split_on' to be a string, got {type(split_on).__name__}")
        if not isinstance(size, int):
            raise TypeError(f"Expected 'size' to be an integer, got {type(size).__name__}")
        if not isinstance(tol, float):
            raise TypeError(f"Expected 'tol' to be a float, got {type(tol).__name__}")
        if not isinstance(relative, bool):
            raise TypeError(f"Expected 'relative' to be a boolean, got {type(relative).__name__}")

        self.features = features
        self.split_on = split_on
        self.size = size
        self.tol = tol
        self.relative = relative

        self._compute_splits()

    def _compute_splits(self):
        '''
        Determines data splits for partitioning.

        '''

        self.counter = 0

        idx = self.features.groupby(by=self.split_on).size().sort_index()

        counts = idx.values
        idx = idx.index

        if self.relative:
            dxs = np.diff(idx) / idx[:-1]
        else:
            dxs = np.diff(idx)

        # if relative, convert tol to absolute
        bins = []
        current_count = counts[0]
        current_bin = [idx[0]]
        self._counts = []

        for i, dx in zip(range(1, len(idx)), dxs):
            if (current_count + counts[i] <= self.size) or (dx <= self.tol):
                current_bin.append(idx[i])
                current_count += counts[i]

            else:
                bins.append(np.array(current_bin))
                self._counts.append(current_count)

                current_bin = [idx[i]]
                current_count = counts[i]

        # Add last unadded bin
        bins.append(np.array(current_bin))
        self._counts.append(current_count)

        self.bounds = np.array([[x.min(), x.max()] for x in bins])

        # Number of partitions in the data
        self.n_partitions = len(bins)

    def __iter__(self):
        return self

    def __next__(self):
        if self.counter < len(self.bounds):
            q = '({} >= {}) & ({} <= {})'.format(self.split_on,
                                                 self.bounds[self.counter][0],
                                                 self.split_on,
                                                 self.bounds[self.counter][1])

            subset = self.features.query(q)

            self.counter += 1
            if len(subset.index) > 1:
                return subset
            else:
                return None

        raise StopIteration

    def map(self, func, processes=1, **kwargs):
        '''
        Maps `func` to each partition, then returns the combined result.

        Parameters
        ----------
        func : function
            Function to apply to partitions.
        processes : int
            Number of parallel processes. If less than 2, a serial mapping is
            applied.
        kwargs
            Keyword arguments passed to `func`.

        Returns
        -------
        :obj:`~pandas.DataFrame`
            Combined result of `func` applied to partitions.

        '''

        # Serial
        if processes < 2:
            result = [func(x, **kwargs) for x in self]

        # Parallel
        else:
            with mp.Pool(processes=processes) as p:
                result = list(p.imap(partial(func, **kwargs), self))

        # Add partition index
        for i in range(len(result)):
            if result[i] is not None:
                result[i]['partition_idx'] = i

        # Combine partitions
        return pd.concat(result, ignore_index=True)

def multi_sample_partition(features, split_on='mz', size=500, tol=25E-6, relative=True):
    '''
    Partitions data along a given dimension. For use with features across
    multiple samples, e.g. in alignment.

    Parameters
    ----------
    features : :obj:`~pandas.DataFrame`
        Input feature coordinates and intensities.
    split_on : str
        Dimension to partition the data.
    size : int
        Target partition size.
    tol : float
        Largest allowed distance between unique `split_on` observations.
    relative : bool
        If `True`, the `tol` parameter is interpreted as a relative tolerance.

    Returns
    -------
    :obj:`~deimos.subset.Partitions`
        A generator object that will lazily build and return each partition.

    '''

    return MultiSamplePartitions(features, split_on, size, tol, relative)
