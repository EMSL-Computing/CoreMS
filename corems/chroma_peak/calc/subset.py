# This file contains functions for subsetting dataframes that contain mass feature data.
# This is based on the deimos package, found here: https://github.com/pnnl/deimos/blob/master/deimos/subset.py with some modifications.

import multiprocessing as mp
from functools import partial

import dask.dataframe as dd
import numpy as np
import pandas as pd


class Partitions:
    '''
    Generator object that will lazily build and return each partition.

    Attributes
    ----------
    features : :obj:`~pandas.DataFrame`
        Input feature coordinates and intensities.
    split_on : str
        Dimension to partition the data.
    size : int
        Target partition size.
    overlap : float
        Amount of overlap between partitions to ameliorate edge effects.

    '''

    def __init__(self, features, split_on='mz', size=1000, overlap=0.05):
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
        overlap : float
            Amount of overlap between partitions to ameliorate edge effects.

        '''

        self.features = features
        self.split_on = split_on
        self.size = size
        self.overlap = overlap

        self._compute_splits()

    def _compute_splits(self):
        '''
        Determines data splits for partitioning.

        '''

        # Unique to split on
        idx = np.unique(self.features[self.split_on].values)

        # Number of partitions
        partitions = np.ceil(len(idx) / self.size)

        # Determine partition bounds
        bounds = [[x.min(), x.max()] for x in np.array_split(idx, partitions)]
        for i in range(1, len(bounds)):
            bounds[i][0] = bounds[i - 1][1] - self.overlap

        if (self.overlap > 0) & (len(bounds) > 1):
            # Functional bounds
            fbounds = []
            for i in range(len(bounds)):
                a, b = bounds[i]

                # First partition
                if i < 1:
                    b = b - self.overlap / 2

                # Middle partitions
                elif i < len(bounds) - 1:
                    a = a + self.overlap / 2
                    b = b - self.overlap / 2

                # Last partition
                else:
                    a = a + self.overlap / 2

                fbounds.append([a, b])
        else:
            fbounds = bounds

        self.bounds = bounds
        self.fbounds = fbounds

    def __iter__(self):
        '''
        Yields each partition.

        Yields
        ------
        :obj:`~pandas.DataFrame`
            Partition of feature coordinates and intensities.

        '''

        for a, b in self.bounds:
            yield slice(self.features, by=self.split_on, low=a, high=b)

    def map(self, func, processes=1, **kwargs):
        '''
        Maps `func` to each partition, then returns the combined result,
        accounting for overlap regions.

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

        # Reconcile overlap
        result = [slice(result[i], by=self.split_on, low=a, high=b)
                  for i, (a, b) in enumerate(self.fbounds)]

        # Combine partitions
        return pd.concat(result).reset_index(drop=True)

    def zipmap(self, func, b, processes=1, **kwargs):
        '''
        Maps `func` to each partition pair resulting from the zip operation of
        `self` and `b`, then returns the combined result, accounting for
        overlap regions.

        Parameters
        ----------
        func : function
            Function to apply to zipped partitions. Must accept and return two
            :obj:`~pandas.DataFrame` instances.
        b : :obj:`~pandas.DataFrame`
            Input feature coordinates and intensities.
        processes : int
            Number of parallel processes. If less than 2, a serial mapping is
            applied.
        kwargs
            Keyword arguments passed to `func`.

        Returns
        -------
        a, b : :obj:`~pandas.DataFrame`
            Result of `func` applied to paired partitions.

        '''

        # Partition other dataset
        partitions = (slice(b, by=self.split_on, low=a, high=b_)
                      for a, b_ in self.bounds)

        # Serial
        if processes < 2:
            result = [func(a, b_, **kwargs) for a, b_ in zip(self, partitions)]

        # Parallel
        else:
            with mp.Pool(processes=processes) as p:
                result = list(p.starmap(partial(func, **kwargs),
                                        zip(self, partitions)))

        result = {'a': [x[0] for x in result], 'b': [x[1] for x in result]}

        # Reconcile overlap
        tmp = [slice(result['a'][i], by=self.split_on, low=a, high=b_,
                     return_index=True)
               for i, (a, b_) in enumerate(self.fbounds)]

        result['a'] = [x[0] for x in tmp]
        idx = [x[1] for x in tmp]
        result['b'] = [p.iloc[i, :] if i is not None else None for p,
                       i in zip(result['b'], idx)]

        # Combine partitions
        result['a'] = pd.concat(result['a'])
        result['b'] = pd.concat(result['b'])

        return result['a'], result['b']


class MultiSamplePartitions:
    '''
    Generator object that will lazily build and return each partition constructed
    from multiple samples.

    Attributes
    ----------
    features : :obj:`~pandas.DataFrame` or :obj:`~dask.dataframe.DataFrame`
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

    def __init__(self, features, split_on='mz', size=500, tol=25E-6):
        '''
        Initialize :obj:`~deimos.subset.Partitions` instance.

        Parameters
        ----------
        features : :obj:`~pandas.DataFrame` or :obj:`~dask.dataframe.DataFrame`
            Input feature coordinates and intensities.
        split_on : str
            Dimension to partition the data.
        size : int
            Target partition size.
        tol : float
            Largest allowed distance between unique `split_on` observations.

        '''

        self.features = features
        self.split_on = split_on
        self.size = size
        self.tol = tol

        if isinstance(features, dd.DataFrame):
            self.dask = True
        else:
            self.dask = False

        self._compute_splits()

    def _compute_splits(self):
        '''
        Determines data splits for partitioning.

        '''

        self.counter = 0

        if self.dask:
            idx = self.features.groupby(
                by=self.split_on).size().compute().sort_index()
        else:
            idx = self.features.groupby(by=self.split_on).size().sort_index()

        counts = idx.values
        idx = idx.index

        dxs = np.diff(idx) / idx[:-1]

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

            if self.dask:
                subset = self.features.query(q).compute()
            else:
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


def partition(features, split_on='mz', size=1000, overlap=0.05):
    '''
    Partitions data along a given dimension.

    Parameters
    ----------
    features : :obj:`~pandas.DataFrame`
        Input feature coordinates and intensities.
    split_on : str
        Dimension to partition the data.
    size : int
        Target partition size.
    overlap : float
        Amount of overlap between partitions to ameliorate edge effects.

    Returns
    -------
    :obj:`~deimos.subset.Partitions`
        A generator object that will lazily build and return each partition.

    '''

    return Partitions(features, split_on, size, overlap)


def multi_sample_partition(features, split_on='mz', size=500, tol=25E-6):
    '''
    Partitions data along a given dimension. For use with features across
    multiple samples, e.g. in alignment.

    Parameters
    ----------
    features : :obj:`~pandas.DataFrame` or :obj:`~dask.dataframe.DataFrame`
        Input feature coordinates and intensities.
    split_on : str
        Dimension to partition the data.
    size : int
        Target partition size.
    tol : float
        Largest allowed distance between unique `split_on` observations.

    Returns
    -------
    :obj:`~deimos.subset.Partitions`
        A generator object that will lazily build and return each partition.

    '''

    return MultiSamplePartitions(features, split_on, size, tol)
