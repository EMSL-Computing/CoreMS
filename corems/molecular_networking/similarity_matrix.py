"""
SimilarityMatrix
================

Sparse matrix storage for pairwise spectral similarity scores.
One SimilarityMatrix instance is maintained per similarity metric
(e.g., "entropy_similarity", "cosine").

Uses scipy.sparse.lil_matrix for incremental construction, then converts
to csr_matrix for efficient row-slicing queries.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.sparse import lil_matrix, csr_matrix, save_npz, load_npz


class SimilarityMatrix:
    """Sparse, symmetric similarity matrix for one spectral similarity metric.

    Parameters
    ----------
    metric_name : str
        Name of the similarity metric stored in this matrix
        (e.g., ``"entropy_similarity"``, ``"cosine"``).

    Attributes
    ----------
    metric_name : str
        Name of the similarity metric.
    _id_to_idx : dict
        Mapping from spectrum ID (str) to integer matrix index.
    _idx_to_id : list
        Mapping from integer matrix index to spectrum ID.
    _matrix : scipy.sparse.lil_matrix or None
        Underlying sparse matrix (lil for construction, csr after finalise).
    _is_csr : bool
        True when the matrix has been converted to CSR format.
    """

    def __init__(self, metric_name: str = "entropy_similarity"):
        self.metric_name = metric_name
        self._id_to_idx: dict[str, int] = {}
        self._idx_to_id: list[str] = []
        self._matrix: lil_matrix | csr_matrix | None = None
        self._is_csr: bool = False

    # ── Internal helpers ──────────────────────────────────────────────────────

    def _ensure_lil(self):
        """Convert to LIL format if currently CSR (needed before adding entries)."""
        if self._is_csr and self._matrix is not None:
            self._matrix = self._matrix.tolil()
            self._is_csr = False

    def _grow_matrix(self, new_size: int):
        """Grow the sparse matrix to accommodate *new_size* spectra."""
        if self._matrix is None:
            self._matrix = lil_matrix((new_size, new_size), dtype=np.float32)
            self._is_csr = False
        else:
            self._ensure_lil()
            old_size = self._matrix.shape[0]
            if new_size > old_size:
                # Resize by creating a new larger matrix and copying
                new_mat = lil_matrix((new_size, new_size), dtype=np.float32)
                cx = self._matrix.tocsr()
                new_mat[:old_size, :old_size] = cx
                self._matrix = new_mat

    # ── Public API ────────────────────────────────────────────────────────────

    def register_spectra(self, spectrum_ids: list[str]):
        """Register new spectrum IDs, growing the matrix as needed.

        Parameters
        ----------
        spectrum_ids : list of str
            IDs to register.  Already-registered IDs are silently skipped.
        """
        new_ids = [sid for sid in spectrum_ids if sid not in self._id_to_idx]
        if not new_ids:
            return
        start = len(self._idx_to_id)
        for sid in new_ids:
            self._id_to_idx[sid] = len(self._idx_to_id)
            self._idx_to_id.append(sid)
        self._grow_matrix(len(self._idx_to_id))

    def set_similarity(self, id1: str, id2: str, score: float):
        """Store a similarity score for a pair of spectra.

        The matrix is symmetric: both (i, j) and (j, i) are set.

        Parameters
        ----------
        id1, id2 : str
            Spectrum IDs (must already be registered).
        score : float
            Similarity score to store.

        Raises
        ------
        KeyError
            If either ID has not been registered.
        """
        i = self._id_to_idx[id1]
        j = self._id_to_idx[id2]
        self._ensure_lil()
        self._matrix[i, j] = score
        self._matrix[j, i] = score

    def get_similarity(self, id1: str, id2: str) -> float:
        """Retrieve the similarity score for a pair of spectra.

        Parameters
        ----------
        id1, id2 : str
            Spectrum IDs.

        Returns
        -------
        float
            Similarity score, or 0.0 if the pair was never stored.
        """
        if id1 not in self._id_to_idx or id2 not in self._id_to_idx:
            return 0.0
        i = self._id_to_idx[id1]
        j = self._id_to_idx[id2]
        if self._matrix is None:
            return 0.0
        return float(self._matrix[i, j])

    def finalise(self):
        """Convert internal LIL matrix to CSR for efficient querying.

        Call this after all ``set_similarity`` calls for a batch are done.
        """
        if self._matrix is not None and not self._is_csr:
            self._matrix = self._matrix.tocsr()
            self._is_csr = True

    def get_pairs_above_threshold(
        self, threshold: float
    ) -> list[tuple[str, str, float]]:
        """Return all stored pairs with score >= *threshold*.

        Parameters
        ----------
        threshold : float
            Minimum score to include.

        Returns
        -------
        list of (id1, id2, score)
            Only the upper-triangle pairs are returned (no duplicates).
        """
        if self._matrix is None:
            return []
        mat = self._matrix.tocsr() if not self._is_csr else self._matrix
        cx = mat.tocoo()
        results = []
        for i, j, v in zip(cx.row, cx.col, cx.data):
            if i < j and v >= threshold:  # upper triangle only
                results.append((self._idx_to_id[i], self._idx_to_id[j], float(v)))
        return results

    def to_dense(self) -> np.ndarray:
        """Return the full similarity matrix as a dense numpy array.

        Returns
        -------
        numpy.ndarray of shape (n, n)
        """
        if self._matrix is None:
            return np.array([])
        return self._matrix.toarray().astype(np.float32)

    def to_dataframe(self, threshold: float = 0.0) -> pd.DataFrame:
        """Return pairs above *threshold* as a pandas DataFrame.

        Parameters
        ----------
        threshold : float
            Minimum score to include.  Default 0.0 (all stored pairs).

        Returns
        -------
        pandas.DataFrame
            Columns: ``["id1", "id2", "score"]``.
        """
        pairs = self.get_pairs_above_threshold(threshold)
        return pd.DataFrame(pairs, columns=["id1", "id2", "score"])

    @property
    def spectrum_ids(self) -> list[str]:
        """List of registered spectrum IDs in index order."""
        return list(self._idx_to_id)

    @property
    def n_spectra(self) -> int:
        """Number of registered spectra."""
        return len(self._idx_to_id)

    # ── Persistence ───────────────────────────────────────────────────────────

    def save(self, path: str):
        """Save the matrix and ID mapping to a ``.npz`` file.

        Parameters
        ----------
        path : str
            Output file path (will be given ``.npz`` extension if absent).
        """
        self.finalise()
        mat = self._matrix if self._matrix is not None else csr_matrix((0, 0))
        save_npz(path, mat)
        # Save ID list alongside
        ids_path = str(path).replace(".npz", "") + "_ids.npy"
        np.save(ids_path, np.array(self._idx_to_id, dtype=object))

    @classmethod
    def load(cls, path: str) -> "SimilarityMatrix":
        """Load a SimilarityMatrix from a ``.npz`` file.

        Parameters
        ----------
        path : str
            Path to the ``.npz`` file saved by :meth:`save`.

        Returns
        -------
        SimilarityMatrix
        """
        mat = load_npz(path)
        ids_path = str(path).replace(".npz", "") + "_ids.npy"
        ids = list(np.load(ids_path, allow_pickle=True))
        obj = cls()
        obj._matrix = mat.tocsr()
        obj._is_csr = True
        obj._idx_to_id = ids
        obj._id_to_idx = {sid: i for i, sid in enumerate(ids)}
        return obj

    def __repr__(self) -> str:
        n = self.n_spectra
        nnz = self._matrix.nnz // 2 if self._matrix is not None else 0
        return (
            f"SimilarityMatrix(metric='{self.metric_name}', "
            f"n_spectra={n}, n_pairs_stored={nnz})"
        )
