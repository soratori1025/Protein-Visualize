"""
TMHMM 2.0 prediction (Sonnhammer, von Heijne & Krogh, ISMB 1998; Krogh, Larsson,
von Heijne & Sonnhammer, J Mol Biol 305:567, 2001).

Interface adapted from pyTMHMM 1.3.6 ``pyTMHMM/api.py`` (MIT, see LICENSE); the
summary numbers are the ones the original TMHMM 2.0 prints for every sequence.
Decoding is Viterbi, as in tmhmm.py / pyTMHMM. (The DTU program uses its own
"1-best" decoder, so a residue or two at a helix end can differ from the web server.)
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Optional

import numpy as np

from .hmm import emission_table, posterior as state_posterior, viterbi
from .model import parse

GROUP_NAMES = ("i", "M", "o")          # inside / membrane / outside
DEFAULT_MODEL = Path(__file__).with_name("TMHMM2.0.model")
MODEL_ENV = "TMHMM_MODEL"               # optional override: path to a TMHMM2.0.model


@dataclass
class TMHMMResult:
    path: str                     # one of i / M / o per residue ('O' folded into 'o')
    posterior: Optional[np.ndarray]   # (L, 3): P(inside), P(membrane), P(outside)
    exp_aa_tmh: float = 0.0       # "Exp number of AAs in TMHs"
    exp_first60: float = 0.0      # "Exp number, first 60 AAs"
    prob_n_in: float = 0.0        # "Total prob of N-in"

    @property
    def helices(self) -> list[tuple[int, int]]:
        """0-based inclusive (start, end) of every predicted TM helix."""
        out, i, n = [], 0, len(self.path)
        while i < n:
            if self.path[i] == "M":
                j = i
                while j + 1 < n and self.path[j + 1] == "M":
                    j += 1
                out.append((i, j))
                i = j + 1
            else:
                i += 1
        return out

    @property
    def possible_signal_peptide(self) -> bool:
        # TMHMM 2.0 prints "POSSIBLE N-term signal sequence" when this exceeds 10
        return self.exp_first60 > 10.0


def model_path() -> Path:
    env = os.environ.get(MODEL_ENV, "").strip()
    return Path(env) if env else DEFAULT_MODEL


@lru_cache(maxsize=4)
def load_model(path: Optional[str] = None):
    """(initial, transitions, emissions, char_map, label_map, name_map) - parsed once."""
    p = Path(path) if path else model_path()
    if not p.is_file():
        raise FileNotFoundError(
            f"TMHMM 2.0 model file not found at {p}. Put TMHMM2.0.model there or point "
            f"the {MODEL_ENV} environment variable at it.")
    _, model = parse(str(p))
    return model


def predict(sequence: str, model=None, compute_posterior: bool = True,
            reference_quirks: bool = False) -> TMHMMResult:
    initial, transitions, emissions, char_map, label_map, _ = model or load_model()
    if not sequence:
        return TMHMMResult("", np.zeros((0, 3)) if compute_posterior else None)
    emit = emission_table(sequence, emissions, char_map)
    states = viterbi(emit, initial, transitions, reference_quirks)
    labels = np.array([label_map[j] for j in range(len(initial))])
    path = "".join(labels[states]).replace("O", "o")
    if not compute_posterior:
        return TMHMMResult(path, None)

    post = state_posterior(emit, initial, transitions, reference_quirks)
    group_of = np.array([GROUP_NAMES.index("M" if lab == "M" else lab.lower())
                         for lab in labels])
    table = np.zeros((len(sequence), 3))
    for g in range(3):
        table[:, g] = post[:, group_of == g].sum(axis=1)
    return TMHMMResult(path, table,
                       exp_aa_tmh=float(table[:, 1].sum()),
                       exp_first60=float(table[:60, 1].sum()),
                       prob_n_in=float(table[0, 0]))
