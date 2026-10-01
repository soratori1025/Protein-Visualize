"""Vendored TMHMM 2.0 engine (pure Python + numpy, no compiler needed).

    model.py   verbatim from pyTMHMM 1.3.6 (MIT) - parser for the model file
    hmm.py     numpy port of pyTMHMM's Cython hmm.pyx (MIT), two upstream bugs fixed
    api.py     predict() -> TMHMMResult (path, posterior, TMHMM summary numbers)
    TMHMM2.0.model   model parameters of TMHMM 2.0 (Krogh et al. 2001) - see LICENSE:
               the parameters come from DTU and are under DTU's academic licence
"""
from .api import DEFAULT_MODEL, MODEL_ENV, TMHMMResult, load_model, predict

__all__ = ["predict", "load_model", "TMHMMResult", "DEFAULT_MODEL", "MODEL_ENV"]
