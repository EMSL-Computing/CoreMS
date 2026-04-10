"""
corems.molecular_networking
===========================

Molecular networking module for CoreMS.

Classes
-------
MolecularNetwork
    Main interface for building and querying molecular networks.
SimilarityMatrix
    Sparse symmetric matrix for one similarity metric.
SimilarityEngine
    Computes pairwise similarities using FlashEntropy + optional extras.
"""

from corems.molecular_networking.similarity_matrix import SimilarityMatrix
from corems.molecular_networking.similarity_engine import SimilarityEngine
from corems.molecular_networking.network_builder import MolecularNetwork

__all__ = ["MolecularNetwork", "SimilarityMatrix", "SimilarityEngine"]
