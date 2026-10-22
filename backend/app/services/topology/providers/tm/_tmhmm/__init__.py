"""Vendored TMHMM 2.0 engine (pure Python + numpy, no compiler needed).

    model.py   verbatim from pyTMHMM 1.3.6 (MIT) - parser for the model file
    hmm.py     numpy port of pyTMHMM's Cython hmm.pyx (MIT), two upstream bugs fixed
    api.py     predict() -> TMHMMResult (path, posterior, TMHMM summary numbers)

The TMHMM2.0.model file is NOT part of this package: the parameters belong to DTU
Health Tech and their academic licence does not permit redistribution.  Download TMHMM
2.0 from DTU and install the model with backend/scripts/install_tmhmm_model.py, or
point the TMHMM_MODEL environment variable at it (see LICENSE).
"""
from .api import DEFAULT_MODEL, MODEL_ENV, TMHMMResult, load_model, model_path, predict

__all__ = ["predict", "load_model", "model_path", "TMHMMResult", "DEFAULT_MODEL", "MODEL_ENV"]
