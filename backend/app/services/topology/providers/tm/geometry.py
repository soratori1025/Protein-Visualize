"""
3D slab geometry TM block.

The implementation lives in ``app.services.topology.topology_predictor`` (one copy
only). This module re-exports it, so ``providers.tm.geometry.GeometryTMProvider`` and
``topology_predictor.GeometryTMProvider`` are the same class.

(The copy that used to live here was an older version: it re-parsed the file with its
own residue filter instead of using the shared ResidueFrame - so its positions did not
match the SS block - and returned ``TMPrediction(boundaries=[])``, which raises a
pydantic ValidationError for chains of fewer than 10 residues.)
"""
from app.services.topology.topology_predictor import GeometryTMProvider

__all__ = ["GeometryTMProvider"]
